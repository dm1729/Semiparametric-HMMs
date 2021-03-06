priorset <- function(R,M,alpha=rep(1,R),beta=rep(1,M)){
  #set prior
  #need to specify K Dirichlet-K vectors for rows of Q
  #need to specify K Dirichlet-M vectors for Omega
  alpha=alpha
  A=t(replicate(R,alpha)) # KxK matrix, each row corresponds to dirichlet weights for row of Q
  B=t(replicate(R,beta)) # KxM matrix, each row corresponds to dirichlet weights for omega|X=row index
  #return(list(A,B))
  return(list(A,B))
}

QGibbs <- function(X,A){#X vector of states, A prior parameters, R states
  R <- dim(A)[1] #recover the R states value
  P <- matrix(0,nrow=R,ncol=R) #reserves space. P will count transitions
  n <- length(X)
  for ( i in c(1:(n-1)) ){
    P[X[i],X[i+1]] <- P[X[i],X[i+1]] + 1 # Counts according to transition
  }
  A <- A+P #New Dirichlet weights
  Q <- matrix(0,nrow=R,ncol=R) #Allocate sapce for new draw of Q
  for ( i in c(1:R) ){
    Q[i,] <- rdirichlet(1,A[i,]) #draws Q from newly updated Dirichlet weights
  }
  return(Q)
}

WGibbs <- function(X,Y,B){ #States X (Count) Data Y Prior parameters B
  R <- dim(B)[1]
  M <- dim(B)[2]
  P <- matrix(0,nrow=R,ncol=M) #P here counts number of X in i, Y in m (to update Dir)
  n <- length(X)
  for (i in c(1:R) ){
    for (j in c(1:M) ){
      P[i,j]=(X==i)%*%(Y==j) #adds one every time both X_t=i and Y_t=j
    }
  }
  B <- B + P #Updates dirichlet weights when observation for bin j, state i comes in
  W <- matrix(0,nrow=R,ncol=M) #pre-allocate matrix
  for ( i in c(1:R) ){
    W[i,] <- rdirichlet(1,B[i,]) #draws W from newly updated Dirichlet weights
  }
  return(W)
}

XSample <- function(Y,Q,W){ #Samples X vector using forward-backward algo
  #Want to sample from the joint distribution of the latent variables, given parameters
  #library(RHmm) #gives access to forward-backward algo, https://rdrr.io/rforge/RHmm/
  M <- ncol(W) #Sets the number of bins based on number of columns of W
  pi <- abs(eigen(t(Q))$vectors[,1])/sum(abs(eigen(t(Q))$vectors[,1])) #stationary dist
  ListW <- split(t(W),rep(1:nrow(W),each=ncol(W))) #puts into list
  dist <- distributionSet(dis="DISCRETE",proba=ListW,labels=paste(c(1:M))) #paste puts class labels in as strings
  HMM <- HMMSet(pi,Q,dist)
  F <- forwardBackward(HMM,Y) #stores forwardbackward variables
  G <- F$Gamma # Marginal probabilities of X_t=i given data, params
  #Debugging: G seems to well reflect pi in the makeup of hidden states
  Xsi <- F$Xsi #List of T, Xsi[[t]]_{rs}=Xsi_{rs}(t)=P(X_t=r,X_{t+1}=s | Y,Param)
  #Chib method better? Unsure
  #Draw X_1,...,X_T sequentially
  LLH <- F$LLH #Gives log likelihood for input parameters, used in label swapping
  R <- dim(Q)[1]
  X <- rep(0,length(Y))
  X[1] <- sample(R,1,prob = G[1,])
  for (i in (2:length(Y))){
    P <- Xsi[[(i-1)]][X[(i-1)],]/G[(i-1),X[(i-1)]] #Proposal for drawing X_i | X_{i-1}
    X[i] <- c(1:R)%*%rmultinom(1,1,P) #dot product extras 1,..,R from indicator
  }
  return(list("X"=X,"LLH"=LLH))
}

SimulateHMM <- function(Q,W,N){ #Generates N samples from param Q,W
  M <- dim(W)[2] #automatically exracts #bins
  pi <- abs(eigen(t(Q))$vectors[,1])/sum(abs(eigen(t(Q))$vectors[,1])) #stationary dist
  ListW <- split(t(W),rep(1:nrow(W),each=ncol(W))) #puts into list
  dist <- distributionSet(dis="DISCRETE",proba=ListW,labels=paste(c(1:M))) #paste puts class labels in as strings
  HM <- HMMSet(pi,Q,dist)#Sets models
  return(HMMSim(N,HM))
}

SimulateHMMNorm <- function(Q,mu,sig2,N){ #Generates N samples from param Q, mu, sigma^2
  pi <- abs(eigen(t(Q))$vectors[,1])/sum(abs(eigen(t(Q))$vectors[,1])) #stationary dist
  dist <- distributionSet(dis="NORMAL",mean=mu,var=sig2) #paste puts class labels in as strings
  HM <- HMMSet(pi,Q,dist)#Sets models
  return(HMMSim(N,HM))
}

QWPosterior <- function(Y,R,M,b,I,Adir=1,Bdir=1,X=NULL){ # Y data R states M bins b burn-in I iterations
  #Initilisation of Prior
  C <- priorset(R,M,rep(Adir,R),rep(Bdir,M) ) #Adir, Bdir prior precision
  A <- C[[1]] #Initial Q Dirichlet weights
  B <- C[[2]] #Initial W Dirichlet weights
  #library(rje)
  #library(RHmm)
  #Initialisation on X
  n <- length(Y)
  if (is.null(X)){
    X <- c(t(rmultinom(n,1,rep(1,R)))%*%c(1:R)) #Initial state vector, drawn randomly
  }
  #Initialisation on Q,W not even required?
  for (i in c(1:b) ){ #Burn in phase (to get good initialisation on X, others don't matter?)
    Q <- QGibbs(X,A) #Draw Q from conditional dist (update dirichlet weights)
    W <- WGibbs(X,Y,B) #Same as Q but with relevant weight update
    X <- XSample(Y,Q,W)$X #Sample states X given Q,W,Y using Forward/Backward
  }
  LQ <- vector("list",I) #Gets a list ready to store the draws from Q, L[[i]] is draw i of Q
  LW <- vector("list",I)
  LLLH <- vector("list",I) #for storing log likelihood
  LX <- vector("list",I+1)
  LX[[1]] <- X
  for (i in c(1:I)){ #Here we will store draws on Q and W
    LQ[[i]] <- QGibbs(LX[[i]],A)
    LW[[i]] <- WGibbs(LX[[i]],Y,B)
    SamplesLLH <- XSample( Y,LQ[[i]],LW[[i]] ) #also computes log likelihood for Q[i] and W[i]
    LLLH[[i]] <- SamplesLLH$LLH
    LX[[i+1]] <- SamplesLLH$X
  }
  return(list("QList"=LQ,"WList"=LW,"XList"=LX,"LLHList"=LLLH))
}

QWPosteriorNoLatent <- function(Y,R,M,b=0,I,Adir=1,Bdir=1,X=NULL){ # Y data R states M bins b burn-in I iterations this code will not output latent states
  #Initilisation of Prior
  C <- priorset(R,M,rep(Adir,R),rep(Bdir,M) ) #Adir, Bdir prior precision
  A <- C[[1]] #Initial Q Dirichlet weights
  B <- C[[2]] #Initial W Dirichlet weights
  #library(gtools)
  #library(RHmm)
  #Initialisation on X
  n <- length(Y)
  if (is.null(X)){
    X <- c(t(rmultinom(n,1,rep(1,R)))%*%c(1:R)) #Initial state vector, drawn randomly
  }
  #Initialisation on Q,W not even required?
  for (i in c(1:b) ){ #Burn in phase (to get good initialisation on X, others don't matter?)
    Q <- QGibbs(X,A) #Draw Q from conditional dist (update dirichlet weights)
    W <- WGibbs(X,Y,B) #Same as Q but with relevant weight update
    X <- XSample(Y,Q,W)$X #Sample states X given Q,W,Y using Forward/Backward
  }
  LQ <- vector("list",I) #Gets a list ready to store the draws from Q, L[[i]] is draw i of Q
  LW <- vector("list",I)
  LLLH <- vector("list",I) #for storing log likelihood
  for (i in c(1:I)){ #Here we will store draws on Q and W
    LQ[[i]] <- QGibbs(X,A)
    LW[[i]] <- WGibbs(X,Y,B)
    SamplesLLH <- XSample( Y,LQ[[i]],LW[[i]] ) #also computes log likelihood for Q[i] and W[i]
    LLLH[[i]] <- SamplesLLH$LLH
    X <- SamplesLLH$X
  }
  return(list("QList"=LQ,"WList"=LW,"LLHList"=LLLH))
}

QWPosteriorSomeLatent <- function(Y,R,M,b=0,I,Adir=1,Bdir=1,X=NULL){ # Y data R states M bins b burn-in I iterations this code will not output latent states
  #Initilisation of Prior
  C <- priorset(R,M,rep(Adir,R),rep(Bdir,M) ) #Adir, Bdir prior precision
  A <- C[[1]] #Initial Q Dirichlet weights
  B <- C[[2]] #Initial W Dirichlet weights
  #library(gtools)
  #library(RHmm)
  #Initialisation on X
  n <- length(Y)
  if (is.null(X)){
    X <- c(t(rmultinom(n,1,rep(1,R)))%*%c(1:R)) #Initial state vector, drawn randomly
  }
  #Initialisation on Q,W not even required?
  for (i in c(1:b) ){ #Burn in phase (to get good initialisation on X, others don't matter?)
    Q <- QGibbs(X,A) #Draw Q from conditional dist (update dirichlet weights)
    W <- WGibbs(X,Y,B) #Same as Q but with relevant weight update
    X <- XSample(Y,Q,W)$X #Sample states X given Q,W,Y using Forward/Backward
  }
  LQ <- vector("list",I) #Gets a list ready to store the draws from Q, L[[i]] is draw i of Q
  LW <- vector("list",I)
  LLLH <- vector("list",I) #for storing log likelihood
  for (i in c(1:I)){ #Here we will store draws on Q and W
    LQ[[i]] <- QGibbs(X,A)
    LW[[i]] <- WGibbs(X,Y,B)
    SamplesLLH <- XSample( Y,LQ[[i]],LW[[i]] ) #also computes log likelihood for Q[i] and W[i]
    LLLH[[i]] <- SamplesLLH$LLH
    X <- SamplesLLH$X
    if (LLLH[[i]]==max(unlist(LLLH)[1:i]) ){ #store latents when encountering new maximum for likelihood
      XMLE <- list("Iteration Number"=i,"states"=X)
    }
  }
  return(list("QList"=LQ,"WList"=LW,"LLHList"=LLLH,"XMLE"=XMLE))
}

QWPosteriorFixState <- function(Y,R,M,b,I,X=NULL){ # Y data R states M bins b burn-in I iterations
  #Initilisation of Prior
  A <- priorset(R,M)[[1]] #Initial Q Dirichlet weights, can be customized
  B <- priorset(R,M)[[2]] #Initial W Dirichlet weights
  #library(rje)
  #library(RHmm)
  #Initialisation on X
  n <- length(Y)
  if (is.null(X)){
    X <- c(t(rmultinom(n,1,rep(1,R)))%*%c(1:R)) #Initial state vector, drawn randomly
  }
  #Initialisation on Q,W not even required?
  for (i in c(1:b) ){ #Burn in phase (to get good initialisation on X, others don't matter?)
    Q <- QGibbs(X,A) #Draw Q from conditional dist (update dirichlet weights)
    W <- WGibbs(X,Y,B) #Same as Q but with relevant weight update
    #X <- XSample(Y,Q,W) #Sample states X given Q,W,Y using Forward/Backward
  }
  LQ <- vector("list",I) #Gets a list ready to store the draws from Q, L[[i]] is draw i of Q
  LW <- vector("list",I)
  for (i in c(1:I)){ #Here we will store draws on Q
    LQ[[i]] <- QGibbs(X,A)
    LW[[i]] <- WGibbs(X,Y,B)
    #X <- XSample(Y,LQ[[i]],LW[[i]])
  }
  return(list("QList"=LQ,"WList"=LW))
}

Distance <- function(x,y){
  D <- sqrt(sum((x-y)^2))
  return(D)
}

LabelSwap <- function(QList,WList,Y,s,A=NULL,B=NULL){#Thins by factor s, swaps labels
  R <- dim(QList[[1]])[1] #Recovers number of hidden states
  M <- dim(WList[[1]])[2] #Recoves number of bins
  if ( is.null(A) ){
    A <- priorset(R,M)[[1]]
  }
  if ( is.null(B) ){
    B <- priorset(R,M)[[2]]
  }
  S <- length(QList)/s #s must divide length(QList)
  QThin <- vector("list",S) #Thinned list, will place permuted Q in
  WThin <- vector("list",S) #Thinned list
  PostDensity <- rep(0,S)
  for (i in c(1:S) ){ #Finds the  posterior mode
    pi <- abs(eigen(t(Q))$vectors[,1])/sum(abs(eigen(t(Q))$vectors[,1])) #same as XSample
    ListW <- split(t(WList[[s*(i-1)+1]]),rep(1:nrow(WList[[s*(i-1)+1]]),each=ncol(WList[[s*(i-1)+1]]))) 
    dist <- distributionSet(dis="DISCRETE",proba=ListW,labels=paste(c(1:M))) 
    HM <- HMMSet(pi,QList[[s*(i-1)+1]],dist)
    LogLike <- forwardBackward(HM,Y)$LLH #Sets log likelihood for these parameters
    Prior <- sum((A-1)*log(QList[[s*(i-1)+1]])) + sum((B-1)*log(WList[[s*(i-1)+1]])) #Sets Prior prob up to const(A,B)
    PostDensity[i] <- LogLike + Prior
  }
  PivotIndex <- which(PostDensity==max(PostDensity))[1] #Finds argmax of posterior
  PivotQ <- QList[[PivotIndex]]
  PivotW <- WList[[PivotIndex]]
  Perms <- permutations(R,R) #Matrix of permutations
  for (i in c(1:S)){
    D <- rep(0,dim(Perms)[1]) #Holds distances for each i to choose best one
    Q <- QList[[s*(i-1)+1]] #Initialise
    W <- WList[[s*(i-1)+1]]
    for (j in c( 1:dim(Perms)[1] ) ){#perms from gtools
      for (k in c(1:R)){
        W[k,] <- WList[[s*(i-1)+1]][Perms[j,k],] #Permute each vector of W[[i]]
        for (l in c(1:R)){
          Q[k,l] <- QList[[s*(i-1)+1]][ Perms[j,k] , Perms[j,l] ] #Permute Q[[i]] accoridng to perm j
        }
      } #End loop over matrix entries for particular perm
      D[j] <- Distance(c(Q,W),c(PivotQ,PivotW)) # Finds distance between Q,W(Q[[i], applied perm j) and Q,W(pivot)
    } #End loop over particular perm: Next apply best perm
    PermIndex <- which(D==min(D)) #Tells us which permutation was optimal for this Q
    WThin[[i]] <- matrix(0, nrow = R, ncol = M) #Create matrix of correct size in list
    QThin[[i]] <- matrix(0, nrow=R, ncol=R) # As above
    for (k in c(1:R)){ #Now we set entry of the list of outputs according to this perm
      WThin[[i]][k,] <- WList[[s*(i-1)+1]][Perms[PermIndex,k],] #Permute each vector of W[[i]]
      for (l in c(1:R)){
        QThin[[i]][k,l] <- QList[[s*(i-1)+1]][ Perms[PermIndex,k] , Perms[PermIndex,l] ] #Permute Q[[i]] according to perm j
      }
    } #End loop over matrix entries
  } #End loop over S, so now all elements of QThin are entered with appropriate permutation
  return(list("QThin"=QThin,"WThin"=WThin))
}

LabelSwapTruth <- function(QList,WList,Q,W,s,A=NULL,B=NULL){#Thins by factor s, swaps labels
  R <- dim(QList[[1]])[1] #Recovers number of hidden states
  M <- dim(WList[[1]])[2] #Recoves number of bins
  if ( is.null(A) ){
    A <- priorset(R,M)[[1]]
  }
  if ( is.null(B) ){
    B <- priorset(R,M)[[2]]
  }
  S <- length(QList)/s #s must divide length(QList)
  QThin <- vector("list",S) #Thinned list, will place permuted Q in
  WThin <- vector("list",S) #Thinned list
  PivotQ <- Q
  PivotW <- W
  Perms <- permutations(R,R) #Matrix of permutations
  for (i in c(1:S)){
    D <- rep(0,dim(Perms)[1]) #Holds distances for each i to choose best one
    Q <- QList[[s*(i-1)+1]] #Initialise
    W <- WList[[s*(i-1)+1]]
    for (j in c( 1:dim(Perms)[1] ) ){#perms from gtools
      for (k in c(1:R)){
        W[k,] <- WList[[s*(i-1)+1]][Perms[j,k],] #Permute each vector of W[[i]]
        for (l in c(1:R)){
          Q[k,l] <- QList[[s*(i-1)+1]][ Perms[j,k] , Perms[j,l] ] #Permute Q[[i]] accoridng to perm j
        }
      } #End loop over matrix entries for particular perm
      D[j] <- Distance(c(Q,W),c(PivotQ,PivotW)) # Finds distance between Q,W(Q[[i], applied perm j) and Q,W(pivot)
    } #End loop over particular perm: Next apply best perm
    PermIndex <- which(D==min(D)) #Tells us which permutation was optimal for this Q
    WThin[[i]] <- matrix(0, nrow = R, ncol = M) #Create matrix of correct size in list
    QThin[[i]] <- matrix(0, nrow=R, ncol=R) # As above
    for (k in c(1:R)){ #Now we set entry of the list of outputs according to this perm
      WThin[[i]][k,] <- WList[[s*(i-1)+1]][Perms[PermIndex,k],] #Permute each vector of W[[i]]
      for (l in c(1:R)){
        QThin[[i]][k,l] <- QList[[s*(i-1)+1]][ Perms[PermIndex,k] , Perms[PermIndex,l] ] #Permute Q[[i]] accoridng to perm j
      }
    } #End loop over matrix entries
  } #End loop over S, so now all elements of QThin are entered with appropriate permutation
  return(list("QThin"=QThin,"WThin"=WThin))
}

LabelSwapLLH <- function(QList,WList,LLHList,s,A=NULL,B=NULL){#Thins by factor s, swaps labels
  R <- dim(QList[[1]])[1] #Recovers number of hidden states
  M <- dim(WList[[1]])[2] #Recoves number of bins
  if ( is.null(A) ){
    A <- priorset(R,M)[[1]]
  }
  if ( is.null(B) ){
    B <- priorset(R,M)[[2]]
  }
  S <- length(QList)/s #s must divide length(QList)
  QThin <- vector("list",S) #Thinned list, will place permuted Q in
  WThin <- vector("list",S) #Thinned list
  PostDensity <- rep(0,S)
  for (i in c(1:S) ){ #Finds the  posterior mode
    Prior <- sum((A-1)*log(QList[[s*(i-1)+1]])) + sum((B-1)*log(WList[[s*(i-1)+1]])) #Sets Prior prob up to const(A,B)
    PostDensity[i] <- LLHList[[s*(i-1)+1]] + Prior
  }
  PivotIndex <- which(PostDensity==max(PostDensity))[1] #Finds argmax of posterior
  PivotQ <- QList[[PivotIndex]]
  PivotW <- WList[[PivotIndex]]
  Perms <- permutations(R,R) #Matrix of permutations
  for (i in c(1:S)){
    D <- rep(0,dim(Perms)[1]) #Holds distances for each i to choose best one
    Q <- QList[[s*(i-1)+1]] #Initialise
    W <- WList[[s*(i-1)+1]]
    for (j in c( 1:dim(Perms)[1] ) ){#perms from gtools
      for (k in c(1:R)){
        W[k,] <- WList[[s*(i-1)+1]][Perms[j,k],] #Permute each vector of W[[i]]
        for (l in c(1:R)){
          Q[k,l] <- QList[[s*(i-1)+1]][ Perms[j,k] , Perms[j,l] ] #Permute Q[[i]] accoridng to perm j
        }
      } #End loop over matrix entries for particular perm
      D[j] <- Distance(c(Q,W),c(PivotQ,PivotW)) # Finds distance between Q,W(Q[[i], applied perm j) and Q,W(pivot)
    } #End loop over particular perm: Next apply best perm
    PermIndex <- which(D==min(D)) #Tells us which permutation was optimal for this Q
    WThin[[i]] <- matrix(0, nrow = R, ncol = M) #Create matrix of correct size in list
    QThin[[i]] <- matrix(0, nrow=R, ncol=R) # As above
    for (k in c(1:R)){ #Now we set entry of the list of outputs according to this perm
      WThin[[i]][k,] <- WList[[s*(i-1)+1]][Perms[PermIndex,k],] #Permute each vector of W[[i]]
      for (l in c(1:R)){
        QThin[[i]][k,l] <- QList[[s*(i-1)+1]][ Perms[PermIndex,k] , Perms[PermIndex,l] ] #Permute Q[[i]] accoridng to perm j
      }
    } #End loop over matrix entries
  } #End loop over S, so now all elements of QThin are entered with appropriate permutation
  return(list("QThin"=QThin,"WThin"=WThin))
}

#LabelSwapLLHPi2 <- function(QList,WList,LLHList,s,A=NULL,B=NULL){#Thins by factor s, swaps labels
#  R <- dim(QList[[1]])[1] #Recovers number of hidden states
#  M <- dim(WList[[1]])[2] #Recoves number of bins
#  if ( is.null(A) ){
#    A <- priorset(R,M)[[1]]
#  }
#  if ( is.null(B) ){
#    B <- priorset(R,M)[[2]]
#  }
#  S <- length(QList)/s #s must divide length(QList)
#  QThin <- vector("list",S) #Thinned list, will place permuted Q in
#  WThin <- vector("list",S) #Thinned list
#  PostDensity <- rep(0,S)
# for (i in c(1:S) ){ #Finds the  posterior mode
#    Prior <- sum((A-1)*log(QList[[s*(i-1)+1]])) + sum((B-1)*log(WList[[s*(i-1)+1]])) #Sets Prior prob up to const(A,B)
#    PostDensity[i] <- LLHList[[s*(i-1)+1]] + Prior
#  }
#  PivotIndex <- which(PostDensity==max(PostDensity))[1] #Finds argmax of posterior
#  PivotQ <- QList[[PivotIndex]]
#  PivotW <- WList[[PivotIndex]]
#  Perms <- permutations(R,R) #Matrix of permutations
#  for (i in c(1:S)){
  #  D <- rep(0,dim(Perms)[1]) #Holds distances for each i to choose best one
  #  Q <- QList[[s*(i-1)+1]] #Initialise
  #  W <- WList[[s*(i-1)+1]]
  #  for (j in c( 1:dim(Perms)[1] ) ){#perms from gtools
    #  for (k in c(1:R)){
       # W[k,] <- WList[[s*(i-1)+1]][Perms[j,k],] #Permute each vector of W[[i]]
       # for (l in c(1:R)){
      #    Q[k,l] <- QList[[s*(i-1)+1]][ Perms[j,k] , Perms[j,l] ] #Permute Q[[i]] accoridng to perm j
     #   }
    #  } #End loop over matrix entries for particular perm
   #   D[j] <- Distance(c(Q,W),c(PivotQ,PivotW)) # Finds distance between Q,W(Q[[i], applied perm j) and Q,W(pivot)
  #  } #End loop over particular perm: Next apply best perm
  #  PermIndex <- which(D==min(D)) #Tells us which permutation was optimal for this Q
  #  WThin[[i]] <- matrix(0, nrow = R, ncol = M) #Create matrix of correct size in list
  #  QThin[[i]] <- matrix(0, nrow=R, ncol=R) # As above
  #  for (k in c(1:R)){ #Now we set entry of the list of outputs according to this perm
    #  WThin[[i]][k,] <- WList[[s*(i-1)+1]][Perms[PermIndex,k],] #Permute each vector of W[[i]]
    #  for (l in c(1:R)){
     #   QThin[[i]][k,l] <- QList[[s*(i-1)+1]][ Perms[PermIndex,k] , Perms[PermIndex,l] ] #Permute Q[[i]] accoridng to perm j
    #  }
   # } #End loop over matrix entries
  #} #End loop over S, so now all elements of QThin are entered with appropriate permutation
 # return(list("QThin"=QThin,"WThin"=WThin))
#}

#LogPriorProb <- function(Q,W,A,B){ #Input transition matrix and weights, output log dirichlet density
  #*UP TO CONSTANTS DEPENDING ON A,B SINCE A,B ARE FIXED IN MAXIMISATION STEP*
  #P <- sum((A-1)*log(Q)) + sum((B-1)*log(W))
#}

PostMean <- function(QList){#Input MCMC sample list QList (or thinned list QThin)
  QMean <- Reduce("+",QList)/length(QList)
  return(QMean)
}

EntryDraws <- function(QList,i,j,s=1){ #Put in QList of post draws. Extract vector of (i,j) entries
  S <- length(QList)/s
  q<-rep( 0 , S )
  for (k in c( 1:S ) ){
    q[k] <- QList[[(k-1)*s+1]][i,j]
  }
  return(q)
}

Bin <- function(Y,M,link=NULL){
  if (is.null(link)==FALSE){
    Y <- link(Y)
  } #if link has been specified, apply it
  BinY <- floor(M*Y)+1 #Gives bin label in ascending order
  BinY[BinY==(M+1)] <- M #Makes sure {1} is not a singleton bin (makes last interval closed both ends)
  return(BinY)
}

MyLink<- function(y,A=0,B=0){
M <- function(yM,AM,BM){ #Define function which we will vectorize
  if (A==B){
    M <- inv.logit(yM)
  }else if (yM < A){
    M <- inv.logit(yM)
  } else if (yM >B){
    M <- inv.logit(yM)
  } else{
    M <- inv.logit(AM)+( (yM-AM)/(BM-AM) )*(inv.logit(BM)-inv.logit(AM))
  }
  return(M)
}
MV <- Vectorize(M) #Vectorize
MyLink <- MV(y,A,B) #Return
return(MyLink)
}

MyLinkAB <- function(A,B){
  Func <- function(y){
    M <- MyLink(y,A,B)
    return(M)
  }
  return(Func)
}


EmissionPosterior <- function(Y,R,QList,Bdir=1,Mpois=exp(1),X=NULL){ # Y data R states M bins b burn-in I iterations QList Q draws
  #Initilisation of Prior
  I <- length(QList)
  #library(gtools)
  #library(RHmm)
  #Initialisation on X
  n <- length(Y)
  if (is.null(X)){
    X <- c(t(rmultinom(n,1,rep(1,R)))%*%c(1:R)) #Initial state vector, drawn randomly
  }
  #LinkA <- Y[order(Y)][floor(length(Y)/20)] #Takes roughly the 5% lower quantile
  #LinkB <- Y[order(Y,decreasing = TRUE)][floor(length(Y)/20)] #Takes approx 5% upper quantile
  #LinkY <- MyLink(Y,-3,3) #Transforms data (uses same link as used before) *input already transformed*
  LM <- vector("list",I) #for storing selected number of bins
  LW <- vector("list",I) #for storting selected weights. LW[[i]] is matrix of dim R x LM[[i]]
  LLLH <- vector("list",I) #for storing log likelihood
  for (i in c(1:I) ){
    Q <- QList[[i]] #Uses previous draw from Q. At this point, could enter "mini chain" loop
    LM[[i]] <- BinCountSample(X,Y,R,Bdir,Link=NULL,Mpois) #Samples number of bins
    B <- priorset(R,M,rep(0,R),rep(Bdir,M) )[[2]] #Gets the prior matrix for that many bins
    YBin <- factor( Bin(Y,M) , c(1:M) ) #Puts into count data for bins selected previously
    LW[[i]] <- WGibbs(X,YBin,B) # Samples weights
    SamplesLLH <- XSample( YBin,Q,LW[[i]] ) #also computes log likelihood for Q[i] and W[i]
    LLLH[[i]] <- SamplesLLH$LLH
    X <- SamplesLLH$X
  }
  return(list("QList"=QList,"MList"=LM,"WList"=LW,"LLHList"=LLLH)) #(Also returns inputs Q to have on hand)
}

FullPi2 <- function(Y,R,I,Adir=1,Bdir=1,Mpois=exp(1),X=NULL){ # Y data in [0,1] R states M bins b burn-in I iterations QList Q draws
  #Initilisation of Prior
  A <- priorset(R,2,rep(Adir,R),rep(0,2) )[[1]] #Sets prior weights for Q
  #library(gtools)
  #library(RHmm)
  #Initialisation on X
  n <- length(Y)
  if (is.null(X)){
    X <- c(t(rmultinom(n,1,rep(1,R)))%*%c(1:R)) #Initial state vector, drawn randomly
  }
  #LinkA <- Y[order(Y)][floor(length(Y)/20)] #Takes roughly the 5% lower quantile
  #LinkB <- Y[order(Y,decreasing = TRUE)][floor(length(Y)/20)] #Takes approx 5% upper quantile
  #LinkY <- MyLink(Y,-3,3) #Transforms data (uses same link as used before) *input already transformed*
  LQ <- vector("list",I) #for storing draws of Q
  LM <- vector("list",I) #for storing selected number of bins
  LW <- vector("list",I) #for storting selected weights. LW[[i]] is matrix of dim R x LM[[i]]
  LLLH <- vector("list",I) #for storing log likelihood
  for (i in c(1:I) ){
    LQ[[i]] <- QGibbs(X,A) #Uses previous draw from Q. At this point, could enter "mini chain" loop
    LM[[i]] <- BinCountSample(X,Y,R,Bdir,Link=NULL,Mpois) #Samples number of bins
    B <- priorset(R,LM[[i]],rep(0,R),rep(Bdir,LM[[i]]) )[[2]] #Gets the prior matrix for that many bins
    YBin <- factor( Bin(Y,LM[[i]]) , c(1:LM[[i]]) ) #Puts into count data for bins selected previously
    LW[[i]] <- WGibbs(X,YBin,B) # Samples weights
    SamplesLLH <- XSample( YBin,LQ[[i]],LW[[i]] ) #also computes log likelihood for Q[i] and W[i]
    LLLH[[i]] <- SamplesLLH$LLH
    X <- SamplesLLH$X
  }
  return(list("QList"=LQ,"MList"=LM,"WList"=LW,"LLHList"=LLLH)) #(Also returns inputs Q to have on hand)
}

BinCountSample <- function(X,Y,R,Bdir=1,Link=NULL,Mpois=exp(1),TruncM=30){#usual input Y in [0,1]
  #First define function for use later
  if (is.null(TruncM) ){
    TruncM <- length(Y)
  }
  
  if ( is.null(Link)==FALSE ){ #(if not already transformed)
    Y <- Link(Y)
  }
  x <- rep(0,TruncM)
  for (i in c((R+1):TruncM)){
    x[i] <- mLogPost(i,X,Y,R,Bdir,Mpois)
  }
  x <- x-max(x[(R+1):TruncM]) #Rescales, doesn't affect sample but stops overflow with exp
  M <- sample( c((R+1):TruncM) , 1 , prob = exp(x)[(R+1):TruncM] )
  #return(list(M,x[(R+1):TruncM]))
  #return(list(x[(R+1):TruncM]))
  return(M)
}

mLogPost <- function(m,X,Y,R,Bdir,Mpois=exp(1)){ #posterior proba *up to normalising constants*
  LogPriorProb <- -Mpois+(m-R-1)*log(Mpois) - lgamma(m-R) #last term for log (M-R-1)!
  YBin <- factor( Bin(Y,m) , c(1:m) )
  P <- matrix(0,R,m)
  for (i in c(1:R) ){
    for (j in c(1:m) ){
      P[i,j]=(X==i)%*%(YBin==j) #adds one every time both X_t=i and YBin_t=j
    }
  }
  LogLike1 <- sum( lgamma(P+Bdir)-lgamma(Bdir) ) #Contribution from the double product
  LogLike2 <- sum( lgamma( Bdir*m )  - lgamma( rowSums(P) + Bdir*m  ) ) #rowsums of P gives me total count for X_t=i
  LogLike3 <- length(X)*log(m)
  return(LogPriorProb+LogLike1+LogLike2+LogLike3)
}

BinCountMH <- function(X,Y,R,M,Bdir,Link,Mpois=exp(1)){
  #First define function for use later
  mLogPost <- function(m){ #log posterior proba *up to normalising/log additive constants* #uses Rfast
    LogPriorProb <- -Mpois+(m-R-1)*log(Mpois) - log(factorial(m-R-1))
    YBin <- factor( Bin(Y,m,Link) , c(1:m) )
    P <- matrix(0,R,m)
    for (i in c(1:R) ){
      for (j in c(1:m) ){
        P[i,j]=(X==i)%*%(YBin==j) #adds one every time both X_t=i and YBin_t=j
      }
    }
    LogLike1 <- sum( Lgamma(P+Bdir)-Lgamma(Bdir) ) #Contribution from the double product in likelihood
    LogLike2 <- sum( Lgamma( Bdir*m )  - Lgamma( rowSums(P) + Bdir*m  ) ) #rowsums of P gives me total count for X_t=i
    return(LogPriorProb+LogLike1+LogLike2)
  }
  Mcand <- R+1+rpois(1,Mpois) #using q(x,y)=Pois(y) no idea what proposal to use!!
  LogPoisOld <- -Mpois+(M-R-1)*log(Mpois) - log(factorial(M-R-1))
  LogPoisCand <- -Mpois+(Mcand-R-1)*log(Mpois) - log(factorial(Mcand-R-1))
  logalpha <- LogPoisOld-LogPoisCand+mLogPost(Mcand)-mLogPost(M)
  alpha <- min(1,exp(logalpha))
  U <- runif(1)
  if (U <= alpha){
    M <- Mcand
  }
  return(list(logalpha,M))
}
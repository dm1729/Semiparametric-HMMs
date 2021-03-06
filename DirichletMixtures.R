MDPTrunc <- function(Y,M,cmu,cvar,igshape,igrate,QList,Xinit=NULL,SMax=NULL,C=10){ #M precision cmucvar params for centre measure
  #Q list of draws from pi1 posterior eps tolerance SMax largest number of Dirichlet components allowed
  
  start_time=Sys.time()
  
  #RECOVERING PARAMETERS OF MODEL + ITERATIONS
  N <- length(Y)
  R <- nrow(QList[[1]]) #Number of states
  L <- length(QList) #QList is from Step 1 of cut posterior algo (implemented previously)
  
  #SETTING TRUNCATION LEVEL
  if (is.null(SMax)){
    SMax <- max( 20,floor(sqrt(N)) )
  }
  #DEFINING LISTS
  ThList <- vector("list",L) #Initialize with prior draws
  WList <- vector("list",L) #Initialize with prior draws
  LLHList <- vector("list",L)
  PresList <- vector("list",L)
  
  #INITIALISATION (l=0)
  
  #Initialise inverse variance from prior
  pres <- rgamma(R,igshape,igrate) 
  
  #Initialize stick breaks from prior
  W <- t(gtools::rdirichlet(R,rep((M/SMax),SMax)))
  
  if (is.null(Xinit)){ #random initialisation of X if none specified
    Xinit <- c(t(rmultinom(N,1,rep(1,R)))%*%c(1:R)) #Random init of X (Posterior MAP of pi1?)
  }
  
  #initialize pointers from prior
  
  
  S <- sample(c(1:SMax),N,replace=TRUE,prob = W[,1]) #Initialise pointers
  S[Xinit==2] <- sample(c(1:SMax),sum(Xinit==2),replace = TRUE,prob=W[,2]) #Change the ones for state 2
  
  #Initial allocation of theta array
  
  Th <- matrix(0,nrow=SMax,ncol=R)
  
  #SAMPLER
  for ( l in c(1:L) ){
    #LLH = -Inf #Initialize log likelihood at -Inf
    X <- Xinit #for l=1 returned from input. for l>1 returned from previous loop
    
    for (c in c(1:C) ){
      for (r in c(1:R) ){ #Theta, Precisions vary state-by-state
        
        #DEBUGGING
        #print(Th)
        #print(X)
        #print(r)
        #print(C)
        #print(unique(S[X==r]))
        
        
        #UPDATE THETA (COMPONENT MEANS)
        for (j in unique(S[X==r]) ){ #go through pairs for which product in notes is non-empty
          #(product over i s.t. X_i=r and S_i=j)
          Th[j,r] <- MDPThSample(j,r,S,X,Y,cmu,cvar,pres)
        }
        Th[-unique(S[X==r]),r] <- rnorm(SMax-length(unique(S[X==r])),cmu,sqrt(cvar)) #PRIOR DRAWS
        
        #UPDATE INVERSE VARIANCE
        
        pres[r] <- MDPPresSample(r,S,X,Y,Th,igshape,igrate)
        
      }#END LOOP OVER R
      
      
      #UPDATE STICK BREAKS
      W <- MDPWSample(SMax,S,X,M,R)
      
      
      #UPDATE LATENT CHAIN X AND LATENT POINTERS S
      Latents <- MDPLatentSample(Y,QList[[l]],W,Th,pres)
      X <- Latents$X #Comment if want oracle version
      S <- Latents$S
      
      
      #if (Latents$LLH>LLH){ # Likelihood of new point is higher
        #STORE MODEL PARAMETERS FOR MLE (MAP) ALONG MINI CHAIN
        #Xinit <- Latents$X
        #LLH <- Latents$LLH #update LLH
      #}
      
      if (c==C){ #stores final draw of chain
        LLHList[[l]] <- Latents$LLH
        WList[[l]] <- W
        ThList[[l]] <- Th
        PresList[[l]] <- pres
      }
      
    }#End loop over minichain
    
  }
  
  end_time=Sys.time()
  
  return( list("Thetas"=ThList,"StickBreaks"=WList,"LogLikes"=LLHList,"Precisions"=PresList,"Elapsed"=end_time-start_time) ) #Then draw is sum( W_i*f(.|theta_i) )
}

MDPTruncStore <- function(Y,M,cmu,cvar,igshape,igrate,QList,Xinit=NULL,SMax=NULL,C=10){ #M precision cmucvar params for centre measure
  #Q list of draws from pi1 posterior eps tolerance SMax largest number of Dirichlet components allowed
  #STORES MINICHAIN DRAWS
  start_time=Sys.time()
  #RECOVERING PARAMETERS OF MODEL + ITERATIONS
  N <- length(Y)
  R <- nrow(QList[[1]]) #Number of states
  L <- length(QList) #QList is from Step 1 of cut posterior algo (implemented previously)
  
  #SETTING TRUNCATION LEVEL
  if (is.null(SMax)){
    SMax <- max( 20,floor(sqrt(N)) )
  }
  #DEFINING LISTS
  ThList <- vector("list",L*C) #Initialize with prior draws
  WList <- vector("list",L*C) #Initialize with prior draws
  LLHList <- vector("list",L*C)
  PresList <- vector("list",L*C)
  #XList <- vector("list",L*C)
  
  #INITIALISATION (l=0)
  
  #Initialise inverse variance from prior
  pres <- rgamma(R,igshape,igrate) 
  
  #Initialize stick breaks from prior
  W <- t(gtools::rdirichlet(R,rep((M/SMax),SMax)))
  
  if (is.null(Xinit)){ #random initialisation of X if none specified
    Xinit <- c(t(rmultinom(N,1,rep(1,R)))%*%c(1:R)) #Random init of X (Posterior MAP of pi1?)
  }
  
  #initialize pointers from prior
  
  S <- sample(c(1:SMax),N,replace=TRUE,prob = W[,1]) #Initialise pointers
  S[X==2] <- sample(c(1:SMax),sum(X==2),replace = TRUE,prob=W[,2]) #Change the ones for state 2
  
  #Initial allocation of theta array
  
  Th <- matrix(0,nrow=SMax,ncol=R)
  
  #SAMPLER
  
  e <- 1
  for ( l in c(1:L) ){
    LLH = -Inf #Initialize log likelihood at -Inf
    X <- Xinit #for l=1 returned from input. for l>1 returned from previous loop
    
    for (c in c(1:C) ){
      for (r in c(1:R) ){ #Theta, Precisions vary state-by-state
        
        #DEBUGGING
        #print(Th)
        #print(X)
        #print(r)
        #print(C)
        #print(unique(S[X==r]))
        
        
        #UPDATE THETA (COMPONENT MEANS)
        for (j in unique(S[X==r]) ){ #go through pairs for which product in notes is non-empty
          #(product over i s.t. X_i=r and S_i=j)
          Th[j,r] <- MDPThSample(j,r,S,X,Y,cmu,cvar,pres)
        }
        Th[-unique(S[X==r]),r] <- rnorm(SMax-length(unique(S[X==r])),cmu,sqrt(cvar)) #PRIOR DRAWS
        
        #UPDATE INVERSE VARIANCE
        
        pres[r] <- MDPPresSample(r,S,X,Y,Th,igshape,igrate)
        
      }#END LOOP OVER R
      
      
      #UPDATE STICK BREAKS
      W <- MDPWSample(SMax,S,X,M,R)
      
      
      #UPDATE LATENT CHAIN X AND LATENT POINTERS S
      Latents <- MDPLatentSample(Y,QList[[l]],W,Th,pres)
      X <- Latents$X #Comment if want oracle version
      S <- Latents$S
      
      
      #if (Latents$LLH>LLH){ # Likelihood of new point is higher
      #STORE MODEL PARAMETERS FOR MLE (MAP) ALONG MINI CHAIN
      #Xinit <- Latents$X
      #LLH <- Latents$LLH #update LLH
      #}
      LLHList[[e]] <- Latents$LLH
      WList[[e]] <- W
      ThList[[e]] <- Th
      PresList[[e]] <- pres
      #XList[[e]] <- X
      
      e <- e+1
      
    }#End loop over minichain
    
  }
  end_time=Sys.time()
  return( list("Thetas"=ThList,"StickBreaks"=WList,"LogLikes"=LLHList,"Precisions"=PresList,"XList"=X,"Elapsed"=end_time-start_time) ) #Then draw is sum( W_i*f(.|theta_i) )
}

#FULL Pi2 (no cut, using Dirichlet prior on Q)

FullPi2 <- function(Y,M,L=5000,cmu,cvar,igshape,igrate,Xinit=NULL,SMax=NULL,C=10){ #M precision cmucvar params for centre measure
  #Q list of draws from pi1 posterior eps tolerance SMax largest number of Dirichlet components allowed
  #STORES MINICHAIN DRAWS
  start_time=Sys.time()
  #RECOVERING PARAMETERS OF MODEL + ITERATIONS
  N <- length(Y)
  R <- 2 #Number of states
  #L <- length(QList) #QList is from Step 1 of cut posterior algo (implemented previously)
  
  #SETTING TRUNCATION LEVEL
  if (is.null(SMax)){
    SMax <- max( 20,floor(sqrt(N)) )
  }
  #DEFINING LISTS
  ThList <- vector("list",L) #Initialize with prior draws
  WList <- vector("list",L) #Initialize with prior draws
  LLHList <- vector("list",L)
  PresList <- vector("list",L)
  QList <- vector("list",L)
  #XList <- vector("list",L*C)
  
  #INITIALISATION (l=0)
  
  #Initialise inverse variance from prior
  pres <- rgamma(R,igshape,igrate) 
  
  #Initialize stick breaks from prior
  W <- t(gtools::rdirichlet(R,rep((M/SMax),SMax)))
  
  if (is.null(Xinit)){ #random initialisation of X if none specified
    Xinit <- c(t(rmultinom(N,1,rep(1,R)))%*%c(1:R)) #Random init of X (Posterior MAP of pi1?)
  }
  
  #initialize pointers from prior
  
  X <- Xinit #for l=1 returned from input. for l>1 returned from previous loop
  S <- sample(c(1:SMax),N,replace=TRUE,prob = W[,1]) #Initialise pointers
  S[X==2] <- sample(c(1:SMax),sum(X==2),replace = TRUE,prob=W[,2]) #Change the ones for state 2
  
  #Initial allocation of theta array
  
  Th <- matrix(0,nrow=SMax,ncol=R)
  
  #SAMPLER
  A <- matrix(1,nrow=R,ncol=R)
  
  e <- 1
  
  for ( l in c(1:L) ){
    Q <- QGibbs(X,A)
    if (e%%2500==0){
      print(Sys.time()-start_time) #speed check
    }
    
      for (r in c(1:R) ){ #Theta, Precisions vary state-by-state
        
        #DEBUGGING
        #print(Th)
        #print(X)
        #print(r)
        #print(C)
        #print(unique(S[X==r]))
        
        
        #UPDATE THETA (COMPONENT MEANS)
        for (j in unique(S[X==r]) ){ #go through pairs for which product in notes is non-empty
          #(product over i s.t. X_i=r and S_i=j)
          Th[j,r] <- MDPThSample(j,r,S,X,Y,cmu,cvar,pres)
        }
        Th[-unique(S[X==r]),r] <- rnorm(SMax-length(unique(S[X==r])),cmu,sqrt(cvar)) #PRIOR DRAWS
        
        #UPDATE INVERSE VARIANCE
        
        pres[r] <- MDPPresSample(r,S,X,Y,Th,igshape,igrate)
        
      }#END LOOP OVER R
      
      
      #UPDATE STICK BREAKS
      W <- MDPWSample(SMax,S,X,M,R)
      
      
      #UPDATE LATENT CHAIN X AND LATENT POINTERS S
      Latents <- MDPLatentSample(Y,Q,W,Th,pres)
      X <- Latents$X #Comment if want oracle version
      S <- Latents$S
      
      
      #if (Latents$LLH>LLH){ # Likelihood of new point is higher
      #STORE MODEL PARAMETERS FOR MLE (MAP) ALONG MINI CHAIN
      #Xinit <- Latents$X
      #LLH <- Latents$LLH #update LLH
      #}
      LLHList[[e]] <- Latents$LLH
      WList[[e]] <- W
      ThList[[e]] <- Th
      PresList[[e]] <- pres
      QList[[e]] <- Q
      #XList[[e]] <- X
      
      e <- e+1
    
  }
  end_time=Sys.time()
  return( list("Thetas"=ThList,"StickBreaks"=WList,"LogLikes"=LLHList,"Precisions"=PresList,"QList"=QList,"Elapsed"=end_time-start_time) ) #Then draw is sum( W_i*f(.|theta_i) )
}

#THETA SAMPLER
MDPThSample <- function(j,r,S,X,Y,cmu,cvar,pres){ #Update thetas ( Step 2b(ii) )
  #Use conjugacy of base measure alpha
  suff <- sum(Y[(S==j)&(X==r)]) #The sufficient stat sum(relevant obs) appearing in the update formula
  sampsize <- sum((S==j)&(X==r)) #The 'effective sample size' for updating this theta
  postvar <- ( (1/cvar^2)+sampsize*pres[r] )^(-1) #Using wikipedia https://en.wikipedia.org/wiki/Conjugate_prior
  postmean <- postvar*( cmu/cvar + suff*pres[r] ) #see also Walker paper
  theta <- rnorm(1,postmean,sqrt(postvar))
  return(theta)
}

#PRECISION SAMPLER
MDPPresSample <- function(r,S,X,Y,Th,igshape,igrate){
  sampsize <- sum(X==r) #effective sample size
  postshape <- igshape + (sampsize/2) #using conjugacy
  postrate <- igrate + 0.5*sum( (Y[X==r]-Th[S[X==r],r])^2 )
  tau <- rgamma(1,postshape,rate=postrate)
  return(tau)
}

#STICK BREAK SAMPLER (DIRICHLET MULTINOMIAL)

MDPWSample <- function(SMax,S,X,M,R){
  W <- matrix(0,nrow=SMax,ncol=R)
  for (r in c(1:R)){
    alpha <- rep(SMax/M,SMax)
    #increases dirichlet weight by one per unit assignment
    #e.g. S[X==r] = (1,1,1,3,3,2) then alpha[1] <- alpha[1] +3, alpha[3] <- alpha[3] + 2 etc.
    alpha[as.numeric(names(table(S[X==r])))]<- alpha[as.numeric(names(table(S[X==r])))]+table(S[X==r])
    W[,r] <- gtools::rdirichlet(1,alpha) # Samples weights
  }
  return(W)
}

#LATENT SAMPLER (STATES AND MIXTURE ALLOCATIONS)

MDPLatentSample <- function(Y,Q,W,Th,Pres){ #eps tol, data, Qmat, V betas, S pointers, Th locations
  C <- nrow(W) #(Recovers SMax)
  R <- ncol(W) #Recovers R
  Q1 <- Q1set(Q,W) #sets transition matrix for bivariate state space
  pi <- abs(eigen(t(Q1))$vectors[,1])/sum(abs(eigen(t(Q1))$vectors[,1])) #stationary dist
  #first state X=1 S=1 second state X=2 S=1 third state X=1 S=2 etc.
  dist <- distributionSet(dis="NORMAL",mean = as.vector(t(Th)),var = rep(Pres^(-1),C) )
  HMM <- HMMSet(pi,Q1,dist)
  F <- forwardBackward(HMM,Y) #stores forwardbackward variables
  G <- F$Gamma # Marginal probabilities of X_t=i given data, params
  #Debugging: G seems to well reflect pi in the makeup of hidden states
  Xsi <- F$Xsi #List of T, Xsi[[t]]_{rs}=Xsi_{rs}(t)=P(X_t=r,X_{t+1}=s | Y,Param)
  #Draw X_1,...,X_T sequentially
  LLH <- F$LLH #Gives log likelihood for input parameters
  
  X <- rep(0,length(Y))
  X[1] <- sample(R*C,1,prob = G[1,])
  for (i in (2:length(Y))){
    P <- Xsi[[(i-1)]][X[(i-1)],]/G[(i-1),X[(i-1)]] #Proposal for drawing (X_i,s_i) | (X_{i-1},s_{i-1})
    X[i] <- sample(R*C,1,prob = P)
    #X[i] <- c(1:R)%*%rmultinom(1,1,P) #dot product extras 1,..,R from indicator
  }
  Xdraws <- Vectorize(mod)(X,R) #gives corresponding latent chain states
  Sdraws <- floor((X-1)/R)+1 #gives corresponding mixture allocation
  return( list("X"=Xdraws,"S"=Sdraws,"LLH"=LLH) )
  #return(X)
}

#AUXILIARY FUNCTIONS FOR LATENT SAMPLER

Q1set <- function(Q,W){ #W is SMax by R
  R <- ncol(W)
  SMax <- nrow(W)
  Q1 <- matrix(0,nrow=R*SMax,ncol=R*SMax) #allocates enlarged transition matrix
  for (i in c(1:(R*SMax)) ){ #find a better way of doing this
    for (j in c(1:(R*SMax)) ){
      Q1[i,j] <- Q[mod(i,R),mod(j,R)]*W[floor((j-1)/R)+1,mod(j,R)] # = P(X_{i+1}|X_i)*P(S_{i+1}|X_{i+1})
    }
  }
  return(Q1)
}

mod <- function(i,R){ #returns i mod R but R=R rather than R=0
  return( (i+(R-1))%%R+1 )
}

#PLOTTER FOR OUTPUT OF MDPTRUNC

MDPFullPlot <- function(Data,IndexSet=NULL){ #Change to preprocess Data
  if (is.null(IndexSet)){
    IndexSet <- c(1:length(Data$Thetas))
  }
  par(mfrow=c(1,1))
  for (j in c(1:2)){
    x <- seq(-5,5,0.01)
    #print(x[1])
    #Unlist Data matrix
    f <- matrix(0,nrow=length(x),ncol=length(IndexSet))
    fup <- rep(0,length(x))
    flow <- rep(0,length(x))
    fmean <- rep(0,length(x))
    for (l in c(1:length(x)) ){
      n <- 0
      for (k in IndexSet ){
        n <- n + 1
        t <- Data$Thetas[[k]][,j]
        w <- Data$StickBreaks[[k]][,j]
        tau <- Data$Precisions[[k]][j]
        for (i in c(1:length(Data$StickBreaks[[k]][,j]))) {
          f[l,n] = f[l,n] + w[i]*( tau^(0.5)*(1/sqrt(2*pi))*exp(-0.5*tau*(x[l]-t[i])^2) )
        }
      }
      #print(f)
      #print(length(x))
      #print(length(Data$Thetas))
      fmean[l] <- mean(f[l,])
      #fup[l] <- quantile(f[l,],0.975)
      #flow[l] <- quantile(f[l,],0.025)
      
    }
    #print(x[which(fmean==max(fmean))])
    #print(max(fup))
    ftrue <- (1/sqrt(2*pi)*1)*exp(-0.5*1*(x-sign(x[which(fmean==max(fmean))]))^2)
    if (j==1){
      if (sign(x[which(fmean==max(fmean))])<0){
        color = "red"
      }
      else{
        color = "blue"
      }
    plot(x,ftrue,col=color,"l",xlab="x",ylab="f(x)",lty=2,lwd=2)
    lines(x,fmean,col=color,lwd=2)
    }
    else if (j==2){
      if (sign(x[which(fmean==max(fmean))])<0){
        color = "red"
      }
      else{
        color = "blue"
      }
      lines(x,ftrue,col=color,"l",xlab="x",ylab="f(x)",lty=2,lwd=2)
      lines(x,fmean,col=color,lwd=2)
    }
    #lines(x,flow,col="green")
    #lines(x,fup, col="green")
    #return(list("fmean"=fmean,"fup"=fup,"flow"=flow))
  }
}

SmoothingPi2 <- function(Pi2In,Pi2Out){ #Pi2 in contains data and Q, Pi2Out contains emissions
  
  I <- length(Pi2Out$Thetas) #should also be length of QThin
  Y <- Pi2In$Y
  M <- floor(sqrt(length(Y)))
  SmoothMat <- matrix(0,nrow = I, ncol = length(Y))
  #print(SmoothMat)
  for (i in c(1:I)){
    #print(i)
    Q <- Pi2In$QList[[i]]
    pi <- abs(eigen(t(Q))$vectors[,1])/sum(abs(eigen(t(Q))$vectors[,1]))
    distr <- distributionSet(dis="MIXTURE",mean=list(Pi2Out$Thetas[[i]][,1],Pi2Out$Thetas[[i]][,2]),var = list(rep(Pi2Out$Precisions[[i]][1],M),rep(Pi2Out$Precisions[[i]][2],M)), proportion = list(Pi2Out$StickBreaks[[i]][,1],Pi2Out$StickBreaks[[i]][,2]) )
    HMM <- HMMSet(pi,Q,distr)
    a <- forwardBackward(HMM,Y)$Gamma[,1]
    #print(a)
    #print(SmoothMat[i,])
    SmoothMat[i,] <- a #P(X=0|data)
  }
  return(SmoothMat) #each row gives the vector P(X_k=1|data) for k=1,...,N=length(Y)
}

SmoothingPi1 <- function(Pi1Data,Index){ #Pi2 in contains data and Q, Pi2Out contains emissions
  I <- length(Pi2Out$Thetas) #should also be length of QThin
  Y <- Pi1Data$Inputs[[Index]]$Y
  SmoothMat <- matrix(0,nrow = I, ncol = length(Y))
  #print(SmoothMat)
  for (i in c(1:I)){
    #print(i)
    Q <- Pi1$QList[[i]]
    pi <- abs(eigen(t(Q))$vectors[,1])/sum(abs(eigen(t(Q))$vectors[,1]))
    #distr <- distributionSet
    HMM <- HMMSet(pi,Q,distr)
    a <- forwardBackward(HMM,Y)$Gamma[,1]
    #print(a)
    #print(SmoothMat[i,])
    SmoothMat[i,] <- a #P(X=0|data)
  }
  return(SmoothMat) #each row gives the vector P(X_k=1|data) for k=1,...,N=length(Y)
}

TrueSmoothing <- function(Y,Q=NULL,means=NULL,vars=NULL){ #takes about 0.002 seconds
  a <- Sys.time()
  if (is.null(Q) ){
    Q <- t(matrix(c(0.7,0.3,0.2,0.8),2,2))
  }
  if (is.null(means)){
    means <- c(-1,1)
  }
  if (is.null(vars)){
    vars <- c(1,1)
  }
  pi <- abs(eigen(t(Q))$vectors[,1])/sum(abs(eigen(t(Q))$vectors[,1]))
  distr <- distributionSet(dis="NORMAL",means,vars)
  HMM <- HMMSet(pi,Q,distr)
  Smooth <- forwardBackward(HMM,Y)$Gamma[,1] #P(X=0|data)
  elapsed <- Sys.time() - a
  return(Smooth)
}

SmoothingPlot <- function(kvec,SmoothMat,TrueSmooth,swap=FALSE){ #index k, returns P(X_k=1)
  if (swap){
    TrueSmooth <- 1-TrueSmooth
  }
  for (k in kvec){
  SmoothDraws <- SmoothMat[,k] #as many entries as iterations of Pi2
  Upper <- quantile(SmoothDraws,0.95)
  Lower <- quantile(SmoothDraws,0.05)
  hist(SmoothDraws,breaks=seq(0,1,0.025),main=paste("Conditional probability of state ",paste(k),"being equal 1, given N=",paste(length(TrueSmooth)),"observations"))
  lines(c(TrueSmooth[k],TrueSmooth[k]),c(0,10000),col="green")
  print(TrueSmooth[k])
  lines(c(Upper,Upper),c(0,10000),col="red")
  print(Upper)
  lines(c(Lower,Lower),c(0,10000),col="red")
  print(Lower)
  }
}

ConfCoverage <- function(SmoothMat,TrueSmooth,swap = FALSE){
  Count <- 0
  if (swap){
    TrueSmooth <- 1-TrueSmooth
  }
  for (k in c(1:length(TrueSmooth))){
    SmoothDraws <- SmoothMat[,k] #as many entries as iterations of Pi2
    Upper <- quantile(SmoothDraws,0.95)
    Lower <- quantile(SmoothDraws,0.05)
    if (Lower < TrueSmooth[k] & TrueSmooth[k] < Upper){
      Count <- Count+1
    }
  }
  Coverage <- Count/length(TrueSmooth)
  return(Coverage)
}

ErrorPlot <- function(SmoothMat,TrueSmooth,Input,swap = FALSE){ #index k, returns P(X_k=1)
  Err = 0*TrueSmooth
  if (swap){
    TrueSmooth <- 1-TrueSmooth
  }
  for (k in c(1:length(TrueSmooth))){
    SmoothDraws <- SmoothMat[,k] #as many entries as iterations of Pi2
    Err[k] = abs(TrueSmooth[k]-mean(SmoothDraws))
    
  }
  plot(c(1:k),Err,xlab="Hidden state index",ylab="Absolute error in estimating hidden state with given index")
  #hist(Err,breaks=seq(0,1,0.025),main=paste("Absolute difference between posterior mean prob. and true prob for N=",length(TrueSmooth)))
}

ClassificationLoss <- function(SmoothMat,TrueSmooth,swap=TRUE){
  C <- 0
  if (swap){
    TrueSmooth <- 1-TrueSmooth
  }
  for (k in c(1:length(TrueSmooth))){
    SmoothDraws <- SmoothMat[,k] #as many entries as iterations of Pi2
    if (mean(SmoothDraws)<0.5 & TrueSmooth[k]<0.5){
      C = C+1
    }
    if (mean(SmoothDraws)>=0.5 & TrueSmooth[k]>=0.5){
      C = C+1
    }
  }
  return(C/length(TrueSmooth))
}
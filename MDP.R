MDPPost <- function(Y,M,cmu,cvar,igpar,QList,SMax=NULL){ #M precision cmucvar params for centre measure
  #Q list of draws from pi1 posterior eps tolerance SMax largest number of Dirichlet components allowed
  N <- length(Y)
  if (is.null(SMax)){
    SMax <- max( 20,sqrt(N) )
  }
  var <- rep(1,R) #variance of mixture comps: need to change so that this has inv gamma prior
  R <- nrow(QList[[1]]) #Number of states
  L <- length(QList) #QList is from Step 1 of cut posterior algo (implemented previously)
  ThList <- vector("list",L) #Initialize with prior draws
  WList <- vector("list",L) #Initialize with prior draws
  V <- matrix(0,nrow=SMax,ncol=R)
  W <- matrix(0,nrow=SMax,ncol=R)
  for (r in c(1:R) ){
    V[,r] <- c(rbeta((SMax-1),1,(M)),1) #Prior draws for V. Make last one =1 to make a unit stick for W
    W[,r] <- V[,r]*c(1,cumprod(1-V[,r]))[1:SMax]
  }
  X <- c(t(rmultinom(n,1,rep(1,R)))%*%c(1:R)) #Random init of X (Posterior MAP of pi1?)
  S <- sample(c(1:SMax),N,replace=TRUE,prob = W[,1]) #Initialise pointers
  S[X==2] <- sample(c(1:SMax),sum(X==2),replace = TRUE,prob=W[,2]) #Change the ones for state 2
  #Edit initialisation to be for general # of states!
  for ( l in c(1:L) ){
    #Step 2a (done outside for loop for l=1 in order to intialise S)
    if (l>1){
      X <- MDPXSample(R,Y,QList[[l]],WList[[l-1]],S,ThList[[l-1]])  
    }
    
    #Step 2bi
    if (l==1){
      U <- MDPUSample(W,X,S)
    }else{
      U <- MDPUSample(W[[(l-1)]],X,S)  
    }
    
    #Step 2bii
    ThList[[l]] <- ThList[[l-1]] #Preallocates for theta update (keeps most entries the same)
    #Do I really want to do this? Or do I want to draw a whole new lot each time?
    #See comment at end. Probably need prior draws to update latents
    for (j in unique(S) ){ #Theta update S can go between 1 and SMax but won't contain most of them
      #Only need to update those j for which s_i=j for some i
      for (r in c(1:R) ){
        ThList[[l]][j,r] <- MDPThSample(j,r,S,X,Y,cmu,cvar,var)
      }
    }
    
    #Step 2biii
    if (l=1){
      for (j in unique(S) ){ #When j is one of the elements of S, j=s_i for some i in product so impacts dist
        #Maybe I want j in c(1:max(S)) instead?
        for (r in c(1:R) ){
          V[j,r] <- MDPVSample(j,r,U,V,W,X,M)
        }
      }
    } else {
    for (j in unique(S) ){ #When j is one of the elements of S, j=s_i for some i in product so impacts dist
      #Maybe I want j in c(1:max(S)) instead?
      for (r in c(1:R) ){
        V[j,r] <- MDPVSample(j,r,U,V,W[[l-1]],X,M)
      }
    }
    }
    WList[[l]]  <- V*c(1,cumprod(1-V))[1:SMax] #Stores the stick associated to most recent V update
    
    #Step 2biv
    for (i in c(1:N) ){ #For updating slices
      S[i] <- MDPSSample(i,U,WList[[l]],X,Y,ThList[[l]]) #need to vectorize!
    }
  }
  
  return( list("Thetas"=ThList,"StickBreaks"=WList) ) #Then draw is sum( W_i*f(.|theta_i) )
  #The 'empty states' can be then filled with prior draws for getting proper posterior draws
  #actually might need to do these on the fly to be able to update X!
}



MDPUSample <- function(W,X,S){ #Update slicing ( Step 2b(i) ) This is just a unif(0,Trunc) for specified trunc point
  N <- length(X)
  if (S[i]==1){
    Trunc <- V[1]
  } else {
  Trunc <- prod( (1-V)[1:(S[i]-1),X[i] ] )*V[S[i],X[i]]
  }
  U <- runif(1,0,Trunc) #Draws from uniform on 0 to W[s_i,X_i]
  return(U)
}

MDPThSample <- function(j,r,S,X,Y,cmu,cvar,var){ #Update thetas ( Step 2b(ii) )
#Use conjugacy of base measure alpha
suff <- sum(Y[S[i==j]&&X[i]==r]) #The sufficient stat sum(relevant obs) appearing in the update formula
sampsize <- sum(Y[S[i==j]&&X[i]==r]) #The 'effective sample size' for updating this theta
postvar <- ( (1/cvar[r]^2)+sampsize/(var[r]^2) )^(-1) #Using wikipedia https://en.wikipedia.org/wiki/Conjugate_prior
postmean <- postvar*( cmu/cvar + suff/var )
theta <- rnorm(1,postmean,postvar)
return(theta)
}

MDPVSample <- function(j,r,U,V,W,X,M){#Update relative stick weights ( 2b(iii) )
  #Use quantile method for sampling from continuous distributions
  #Set endpoints as per pg 109 of vdV (but only looking at the relevant for each state)
  a <- max( U[ S==j&&X==r ] )/(cumprod(1-V[,r])[j-1]) #truncation lower bound 
  b <- min( 1- U[X==r&&S>j]/( V[S[X==r&&S>j,r]]*cumprod(1-V[,r])[S[X==r&&S>j]-1] )  )
  Vdraw <- truncdist::rtrunc(1,"beta",a,b,shape1=1,shape2=M) #library(truncdist)
  return(Vdraw)
}

MDPSSample <- function(i,U,W,X,Y,Th){#Update pointer variables ( 2b(iv) )
  Trunc <- max(which(W[,X[i]]>U[i])) #Tells us the max W s_i val we can accept given that W_{s_i}>U_i
  prob <- exp( ( Y[i]-Th[1:Trunc,X[i]] ) ) #We can point anywhere compatible with the indicator condition
  s <- sample(c(SMin:SMax),1,replace=TRUE,prob=prob) #Samples SMin to truncation SMax according to density
  return(s)
}

#JUST TWO STATES FOR NOW (just need to change the distributionset line but don't want inefficient code)
MDPXSample <- function(R,Y,Q,W,S,Th){ #eps tol, data, Qmat, V betas, S pointers, Th locations
  sig2 <- 1 #CHANGE TO UPDATE WITH INVERSE GAMMA
  C <- nrow(W)
  pi <- abs(eigen(t(Q))$vectors[,1])/sum(abs(eigen(t(Q))$vectors[,1])) #stationary dist
  dist <- distributionSet(dis="MIXTURE",mean=list(Th[,1],Th[,2]),var = list(rep(sig2,C),rep(sig2,C), W )
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
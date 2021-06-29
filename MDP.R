MDPPost <- function(Y,M,cmu,cvar,igshape,igrate,QList,SMax=NULL){ #M precision cmucvar params for centre measure
  #Q list of draws from pi1 posterior eps tolerance SMax largest number of Dirichlet components allowed
  N <- length(Y)
  if (is.null(SMax)){
    SMax <- max( 20,sqrt(N) )
  }
  pres <- rgamma(R,igshape,igrate) #inv variance of mixture comps
  R <- nrow(QList[[1]]) #Number of states
  L <- length(QList) #QList is from Step 1 of cut posterior algo (implemented previously)
  ThList <- vector("list",L) #Initialize with prior draws
  WList <- vector("list",L) #Initialize with prior draws
  LLHList <- vector("list",L)
  PresList <- vector("list",L)
  V <- matrix(0,nrow=SMax,ncol=R)
  W <- matrix(0,nrow=SMax,ncol=R)
  for (r in c(1:R) ){
    V[,r] <- c(rbeta((SMax-1),1,(M)),1) #Prior draws for V. Make last one =1 to make a unit stick for W
    W[,r] <- V[,r]*c(1,cumprod(1-V[,r]))[1:SMax]
  }
  X <- c(t(rmultinom(N,1,rep(1,R)))%*%c(1:R)) #Random init of X (Posterior MAP of pi1?)
  S <- sample(c(1:SMax),N,replace=TRUE,prob = W[,1]) #Initialise pointers
  S[X==2] <- sample(c(1:SMax),sum(X==2),replace = TRUE,prob=W[,2]) #Change the ones for state 2
  #Edit initialisation to be for general # of states!
  for ( l in c(1:L) ){
    #Step 2a (done outside for loop for l=1 in order to intialise S)
    if (l>1){
      StatesLLH <- MDPXSample(R,Y,QList[[l-1]],WList[[l-1]],S,ThList[[l-1]],PresList[[l-1]])
      X <- StatesLLH$X
      LLHList[[l-1]] <- StatesLLH$LLH
    }
    
    #Step 2bi
    if (l==1){
      U <- MDPUSample(W,X,S)
    }else{
      U <- MDPUSample(W[[(l-1)]],X,S)  
    }
    
    #Step 2bii
    
    if (l==1){
      for (j in unique(S) ){ #Theta update S can go between 1 and SMax but won't contain most of them
        #Only need to update those j for which s_i=j for some i
        for (r in c(1:R) ){
          ThList[[l]][j,r] <- MDPThSample(j,r,S,X,Y,cmu,cvar,pres)
        }
      }
    }else{
      ThList[[l]] <- ThList[[l-1]] #Preallocates for theta update (keeps most entries the same)
      #Update empty entries with new prior draws (FILL IN)
      #INSERT PRIOR DRAWS HERE
      for (j in unique(S) ){ #Theta update S can go between 1 and SMax but won't contain most of them
        #Only need to update those j for which s_i=j for some i
        for (r in c(1:R) ){
          ThList[[l]][j,r] <- MDPThSample(j,r,S,X,Y,cmu,cvar,PresList[[l-1]])
        }
      }
    }
    
    #UPDATE TAU
    PresList[[l]] <- rep(0,R)
    for (r in c(1:R) ){
      PresList[[l]][r] <- MDPPresSample(r,S,X,Y,ThList[[l]],igshape,igrate)
    } 
    
    #Step 2biii
    if (l==1){
      for (j in unique(S) ){ #When j is one of the elements of S, j=s_i for some i in product so impacts dist
        #Maybe I want j in c(1:max(S)) instead?
        for (r in unique(X[S==j]) ){ #only loop over those r for which U[S == j && X == r] not empty
          V[j,r] <- MDPVSample(j,r,U,V,W,S,X,M)
        }
      }
    } else {
    for (j in unique(S) ){ #When j is one of the elements of S, j=s_i for some i in product so impacts dist
      #Maybe I want j in c(1:max(S)) instead?
      for (r in unique(X[S==j]) ){
        V[j,r] <- MDPVSample(j,r,U,V,W[[l-1]],S,X,M)
      }
    }
    }
    #ALSO WANT SOME PRIOR DRAWS. DO V[-UNIQUE(S),r] (need state by state too)
    WList[[l]]  <- V*c(1,cumprod(1-V))[1:SMax] #Stores the stick associated to most recent V update
    
    #Step 2biv
    for (i in c(1:N) ){ #For updating slices
      S[i] <- MDPSSample(i,U,WList[[l]],X,Y,ThList[[l]],PresList[[l]]) #need to vectorize!
    }
  }
  
  return( list("Thetas"=ThList,"StickBreaks"=WList,"LogLikes"-LLHList) ) #Then draw is sum( W_i*f(.|theta_i) )
  #The 'empty states' can be then filled with prior draws for getting proper posterior draws
  #actually might need to do these on the fly to be able to update X!
}



MDPUSample <- function(W,X,S){ #Update slicing ( Step 2b(i) ) This is just a unif(0,Trunc) for specified trunc point
  N <- length(X)
  U <- runif(N) #Draws from uniform on 0 to W[s_i,X_i]
  for (r in c(1:2) ){ #make R-dep. Not sure how to do without loop over r?
  U[X==r] <- U[X==r]*W[S[X==r],r] #Want to make it at most W[S_i,X_i]
  }
  return(U)
}

MDPThSample <- function(j,r,S,X,Y,cmu,cvar,pres){ #Update thetas ( Step 2b(ii) )
#Use conjugacy of base measure alpha
suff <- sum(Y[(S==j)&(X==r)]) #The sufficient stat sum(relevant obs) appearing in the update formula
sampsize <- sum((S==j)&(X==r)) #The 'effective sample size' for updating this theta
postvar <- ( (1/cvar[r]^2)+sampsize*pres(r) )^(-1) #Using wikipedia https://en.wikipedia.org/wiki/Conjugate_prior
postmean <- postvar*( cmu/cvar + suff*pres )
theta <- rnorm(1,postmean,postvar)
return(theta)
}

MDPPresSample <- function(r,S,X,Y,Th,igshape,igrate){
  sampsize <- sum(X==r) #effective sample size
  postshape <- igshape + (sampsize/2) #using conjugacy
  postrate <- igrate + 0.5*sum(Y[X==r]-Theta[S[X==r],r])
  tau <- rgamma(1,postshape,rate=postrate)
}

MDPVSample <- function(j,r,U,V,W,S,X,M){#Update relative stick weights ( 2b(iii) )
  #Use quantile method for sampling from continuous distributions
  #Set endpoints as per pg 109 of vdV (but only looking at the relevant for each state)
  a <- max( U[ S==j&X==r ] )/((cumprod(1-V[,r])[j])/(1-V[j,r])) #truncation lower bound
  b <- min( 1- U[X==r&S>j]/( V[S[X==r&S>j],r]*(cumprod(1-V[,r])[S[X==r&S>j]-1]/(1-V[j,r])) )  )
  #check implementation of b
  Vdraw <- truncdist::rtrunc(1,"beta",a,b,shape1=1,shape2=M) #library(truncdist)
  return(Vdraw)
}

MDPSSample <- function(i,U,W,X,Y,Th,Pres){#Update pointer variables ( 2b(iv) )
  Trunc <- max(which(W[,X[i]]>U[i])) #Tells us the max W s_i val we can accept given that W_{s_i}>U_i
  prob <- sqrt(Pres[X[i]])*exp( -0.5*Pres[X[i]]*( Y[i]-Th[1:Trunc,X[i]] )^2 ) #We can point anywhere compatible with the indicator condition
  #makes a vector for each index of Theta (possible vals of s_i)
  s <- sample(c(1:Trunc),1,replace=TRUE,prob=prob) #Samples SMin to truncation SMax according to density
  return(s)
}

#JUST TWO STATES FOR NOW (just need to change the distributionset line but don't want inefficient code)
MDPXSample <- function(R,Y,Q,W,S,Th,Pres){ #eps tol, data, Qmat, V betas, S pointers, Th locations
  C <- nrow(W)
  pi <- abs(eigen(t(Q))$vectors[,1])/sum(abs(eigen(t(Q))$vectors[,1])) #stationary dist
  dist <- distributionSet(dis="MIXTURE",mean=list(Th[,1],Th[,2]),var = list(rep(Pres[1]^(-1),C),rep(Pres[2]^(-1),C), W ) )
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
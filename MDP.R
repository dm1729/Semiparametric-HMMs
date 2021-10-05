MDPTrunc <- function(Y,M,cmu,cvar,igshape,igrate,QList,X=NULL,SMax=NULL,C=10){ #M precision cmucvar params for centre measure
  #Q list of draws from pi1 posterior eps tolerance SMax largest number of Dirichlet components allowed
  
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
  
  if (is.null(X)){ #random initialisation of X if none specified
    X <- c(t(rmultinom(N,1,rep(1,R)))%*%c(1:R)) #Random init of X (Posterior MAP of pi1?)
  }
  #initialize pointers from prior
  
  S <- sample(c(1:SMax),N,replace=TRUE,prob = W[,1]) #Initialise pointers
  S[X==2] <- sample(c(1:SMax),sum(X==2),replace = TRUE,prob=W[,2]) #Change the ones for state 2
  
  #Initial allocation of theta array
  
  Th <- matrix(0,nrow=SMax,ncol=R)
  
  #SAMPLER
  
  for ( l in c(1:L) ){
    LLHList[[l]] = -Inf #Initialize log likelihood at -Inf
    for (c in c(1:C) ){
      for (r in c(1:R) ){ #Theta, Precisions vary state-by-state
        
        #UPDATE THETA
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
      
      
      #UPDATE POINTERS
      for (i in c(1:N) ){ #need to vectorize!
        S[i] <- MDPSSample(i,SMax,W,X,Y,Th,pres)
      }
      
      
      #UPDATE LATENTS
      StatesLLH <- MDPXSample(R,Y,QList[[l]],W,S,Th,pres)
      #X <- StatesLLH$X #Comment if want oracle version
      if (StatesLLH$LLH>LLHList[[l]]){ # Likelihood of new point is higher
        #STORE MODEL PARAMETERS FOR MLE (MAP) ALONG MINI CHAIN
        LLHList[[l]] <- StatesLLH$LLH
        WList[[l]] <- W
        ThList[[l]] <- Th
        PresList[[l]] <- pres
      }
      
    }#End loop over minichain
    
  }
  
  return( list("Thetas"=ThList,"StickBreaks"=WList,"LogLikes"=LLHList,"Precisions"=PresList) ) #Then draw is sum( W_i*f(.|theta_i) )
}

MDPSlice <- function(Y,M,cmu,cvar,igshape,igrate,QList,X=NULL,SMax=NULL,C=10){ #M precision cmucvar params for centre measure
  #Q list of draws from pi1 posterior eps tolerance SMax largest number of Dirichlet components allowed
  
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
  V <- matrix(0,nrow=SMax,ncol=R)
  W <- matrix(0,nrow=SMax,ncol=R)
  
  for (r in c(1:R) ){
    V[,r] <- c(rbeta((SMax-1),1,(M)),1) #Prior draws for V. Make last one =1 to make a unit stick for W
    W[,r] <- V[,r]*c(1,cumprod(1-V[,r]))[1:SMax]
  }
  
  if (is.null(X)){ #random initialisation of X if none specified
    X <- c(t(rmultinom(N,1,rep(1,R)))%*%c(1:R)) #Random init of X (Posterior MAP of pi1?)
  }
  
  #initialize pointers from prior
  
  S <- sample(c(1:SMax),N,replace=TRUE,prob = W[,1]) #Initialise pointers
  S[X==2] <- sample(c(1:SMax),sum(X==2),replace = TRUE,prob=W[,2]) #Change the ones for state 2
  
  #Initial allocation of theta array
  
  Th <- matrix(0,nrow=SMax,ncol=R)
  
  #SAMPLER
  
  for ( l in c(1:L) ){
    LLHList[[l]] = -Inf #Initialize log likelihood at -Inf
    
    for (c in c(1:C) ){
      #UPDATE SLICES
      U <- MDPUSample(W,X,S)
      
      
      for (r in c(1:R) ){ #Theta, Precisions and Stick-breaks vary state-by-state
        
        #UPDATE THETA
        for (j in unique(S[X==r]) ){ #go through pairs for which product in notes is non-empty
          #(product over i s.t. X_i=r and S_i=j)
          Th[j,r] <- MDPThSample(j,r,S,X,Y,cmu,cvar,pres)
        }
        Th[-unique(S[X==r]),r] <- rnorm(SMax-length(unique(S[X==r])),cmu,sqrt(cvar)) #PRIOR DRAWS

      #UPDATE TAU

        pres[r] <- MDPPresSample(r,S,X,Y,Th,igshape,igrate)
      
      #UPDATE STICK BREAKS
        for ( j in c(1:max(S[X==r])) ){ #the distinct levels which are occupied for that state
          V[j,r] <- MDPVSample(j,r,U,V,S,X,M)
        }
        V[-c(1:max(S[X==r])),r] <- rbeta(SMax-max(S[X==r]),1,M) #prior draws for rest
        V[SMax,r] <- 1
        W[,r]  <- V[,r]*c(1,cumprod(1-V[1:(SMax-1),r]))[1:SMax] 
      
      }#END LOOP OVER R
      
      #UPDATE POINTERS
      for (i in c(1:N) ){ #need to vectorize!
        S[i] <- MDPSSampleSlice(i,U,W,X,Y,Th,pres) 
      }
      
      #UPDATE LATENTS
      StatesLLH <- MDPXSample(R,Y,QList[[l]],W,S,Th,pres)
      #X <- StatesLLH$X #Comment if want oracle version
      if (StatesLLH$LLH>LLHList[[l]]){ # Likelihood of new point is higher
        #STORE MODEL PARAMETERS FOR MLE (MAP) ALONG MINI CHAIN
        LLHList[[l]] <- StatesLLH$LLH
        WList[[l]] <- W
        ThList[[l]] <- Th
        PresList[[l]] <- pres
      }
      
    }#End loop over minichain

  }
  
  return( list("Thetas"=ThList,"StickBreaks"=WList,"LogLikes"=LLHList,"Precisions"=PresList) ) #Then draw is sum( W_i*f(.|theta_i) )
}

#POINTER SAMPLE (DIRICHLET MULTINOMIAL)


MDPSSample <- function(i,SMax,W,X,Y,Th,pres){
  
    #prob of each s proportional to associated p(s|W)*p(Y|s,W)
    prob <- rep(0,SMax)
    
    for (s in c(1:SMax) ){
    prob[s] <- W[s,X[i]]*(sqrt(pres[X[i]])/sqrt(2*pi))*exp(-0.5*pres[X[i]]*(Y[i]-Th[s,X[i]])^2)
    }
    #Update S value
    SUpdate <- sample(c(1:SMax),1,replace=TRUE,prob=prob) #Need to update components for which X==r
    return(SUpdate)
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

#SLICE SAMPLER
MDPUSample <- function(W,X,S){ #Update slicing ( Step 2b(i) ) This is just a unif(0,Trunc) for specified trunc point
  N <- length(X)
  U <- runif(N) #Draws from uniform on 0 to W[s_i,X_i]
  for (r in c(1:ncol(W)) ){ #make R-dep. Not sure how to do without loop over r?
  U[X==r] <- U[X==r]*W[S[X==r],r] #Want to make it at most W[S_i,X_i]
  #Problem where we draw an exact zero if we take W[22,r]?
  }
  return(U)
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
}


#V SAMPLES WITH SLICING
MDPVSample <- function(j,r,U,V,S,X,M){#Update relative stick weights ( 2b(iii) )
  #Use quantile method for sampling from continuous distributions
  #Set endpoints as per pg 109 of vdV (but only looking at the relevant for each state)
  if (sum( (S==j)&(X==r) )>0){ #if there exist some relevant points
  #a <- ((1-V[j,r])*max( (U[ S==j&X==r ] )))/(cumprod(1-V[,r])[j]) #truncation lower bound
  a <- (max( (U[ S==j&X==r ] )))/(min(1,cumprod(1-V[,r])[j-1])) #truncation lower bound. when j=1 takes 1
  }else{
  a <- 0
  }
  #Check a in the case for j<max but not one of the s_i (gap)
  #I suppose just take 0 (hence the max with 0)
  if (sum( (S>j)&(X==r) )>0 ){
  b <- min(1- ( (U[X==r&S>j])*( (1-V[j,r]) ) )/(V[S[X==r&S>j],r]*(cumprod(1-V[,r])[S[X==r&S>j]-1]) ) )
 
  }else{
  b <- 1  
  }
  #check b in the case for j=max (so min is empty)
  #I suppose just take 1 (hence the min with 1)
  
  if ( abs(a-b) <10^(-5) ) { #truncdist doesnt work if R thinks a=b
    #issues here with getting a>b? weird. will check out
    #issues with a and b not even defined??
    Vdraw <- a
  } else{
  Vdraw <- truncdist::rtrunc(1,"beta",a,b,shape1=1,shape2=M) #library(truncdist)
  }
  return(Vdraw)
}


#POINTER SAMPLE WITH SLICES
MDPSSampleSlice <- function(i,U,W,X,Y,Th,Pres){#Update pointer variables ( 2b(iv) )
  Trunc <- max(which(W[,X[i]]>=U[i])) #Tells us the max W s_i val we can accept given that W_{s_i}>U_i
  #changed to geq because of zero bug
  if ( max(W[,X[i]]-U[i]) <=0 ){ #debug loop
    print(max(W[,X[i]]-U[i]) )
  }
   #there should be one of these over zero so max should remain positive.
  #print(Trunc)
  #getting no non-missing arguments error, but how?
  prob <- sqrt(Pres[X[i]])*exp( -0.5*Pres[X[i]]*( Y[i]-Th[1:Trunc,X[i]] )^2 ) #We can point anywhere compatible with the indicator condition
  #makes a vector for each index of Theta (possible vals of s_i)
  s <- sample(c(1:Trunc),1,replace=TRUE,prob=prob) #Samples SMin to truncation SMax according to density
  return(s)
}

#JUST TWO STATES FOR NOW (just need to change the distributionset line but don't want inefficient code)
MDPXSample <- function(R,Y,Q,W,S,Th,Pres){ #eps tol, data, Qmat, V betas, S pointers, Th locations
  C <- nrow(W) #(Recovers SMax)
  pi <- abs(eigen(t(Q))$vectors[,1])/sum(abs(eigen(t(Q))$vectors[,1])) #stationary dist
  dist <- distributionSet(dis="MIXTURE",mean=list(Th[,1],Th[,2]),var = list(rep(Pres[1]^(-1),C),rep(Pres[2]^(-1),C)), list(W[,1],W[,2]) )
  HMM <- HMMSet(pi,Q,dist)
  F <- forwardBackward(HMM,Y) #stores forwardbackward variables
  G <- F$Gamma # Marginal probabilities of X_t=i given data, params
  #Debugging: G seems to well reflect pi in the makeup of hidden states
  Xsi <- F$Xsi #List of T, Xsi[[t]]_{rs}=Xsi_{rs}(t)=P(X_t=r,X_{t+1}=s | Y,Param)
  #Draw X_1,...,X_T sequentially
  LLH <- F$LLH #Gives log likelihood for input parameters
  R <- dim(Q)[1]
  X <- rep(0,length(Y))
  X[1] <- sample(R,1,prob = G[1,])
  for (i in (2:length(Y))){
    P <- Xsi[[(i-1)]][X[(i-1)],]/G[(i-1),X[(i-1)]] #Proposal for drawing X_i | X_{i-1}
    X[i] <- c(1:R)%*%rmultinom(1,1,P) #dot product extras 1,..,R from indicator
  }
  return( list("X"=X,"LLH"=LLH) )
}



####PLOTS###

MDPMLEPlot <- function(Data){#MLE plot
  for (j in c(1:2)){
  m <- which(unlist(Data$LogLikes)==max(unlist(Data$LogLikes)))[1]
  t <- Data$Thetas[[m]][,j]
  w <- Data$StickBreaks[[m]][,j]
  tau <- Data$Precisions[[m]][j]
  x <- seq(-5,5,0.001)
  y <- x*0
  for (i in c(1:length(w) )) {
  y <- y+ w[i]*( tau^(0.5)*(1/sqrt(2*pi))*exp(-0.5*tau*(x-t[i])^2) )
  }
  plot(x,y)
  print(tau)
  }
}

MDPFullPlot <- function(Data){ #Change to preprocess Data
  for (j in c(1:2)){
    x <- seq(-5,5,0.01)
    #print(x[1])
    #Unlist Data matrix
    f <- matrix(0,nrow=length(x),ncol=(length(Data$Thetas)-100))
    fup <- rep(0,length(x))
    flow <- rep(0,length(x))
    fmean <- rep(0,length(x))
    for (l in c(1:length(x)) ){
      for (k in c(101:length(Data$Thetas))){
        t <- Data$Thetas[[k]][,j]
        w <- Data$StickBreaks[[k]][,j]
        tau <- Data$Precisions[[k]][j]
        for (i in c(1:length(Data$StickBreaks[[k]][,j]))) {
          f[l,(k-100)] = f[l,(k-100)] + w[i]*( tau^(0.5)*(1/sqrt(2*pi))*exp(-0.5*tau*(x[l]-t[i])^2) )
          }
      }
      #print(f)
      #print(length(x))
      #print(length(Data$Thetas))
      fmean[l] <- mean(f[l,])
      fup[l] <- quantile(f[l,],0.95)
      flow[l] <- quantile(f[l,],0.05)
      
    }
    print(x[which(fmean==max(fmean))])
    print(max(fup))
    ftrue <- (1/sqrt(2*pi))*exp(-0.5*(x-sign(x[which(fmean==max(fmean))]))^2)
    plot(x,fup,col="red","l")
    lines(x,fmean,col="blue")
    lines(x,flow,col="red")
    lines(x,ftrue, col="green")
    #return(list("fmean"=fmean,"fup"=fup,"flow"=flow))
  }
}
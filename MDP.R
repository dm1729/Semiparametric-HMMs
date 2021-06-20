MDPPost <- function(Y,M,cmu,cvar,igpar,R,QList,eps,SMax){ #M precision cmucvar params for centre measure
  #R states Q list of draws from pi1 posterior eps tolerance SMax largest number of Dirichlet components allowed
  L <- length(QList)
  ThList <- vector("list",L) #Initialize with prior draws
  VList <- vector("list",L) #Initialize with prior draws
  ThList[[]]
  N <- length(Y)
  S <- rep(1,n)
  V <- c( rep(0.5,VMax-1),1 )
  U <- rep(0,N) #Initialise slicing variables
  for ( l in c(1:L) ){
    #Step 2a
    if (l==1){
    X <- MDPXSample(eps,Y,QList[[l]],V,S,Th)
    }else{
      X <- MDPXSample(eps,Y,QList[[l]],VList[[l]],S,ThList[[l]])  
    }
    #Step 2bi
    for (i in c(1:N) ){ #For updating slices
      if (l==1){
        U[i] <- MDPUSample(i,V,X,S)
      }else{
        U[i] <- MDPUSample(i,VList[[l]],X,S)
      }
    }
    #Step 2bii
    ThList[[l]] <- ThList[[l-1]] #Preallocates for theta update (keeps most entries the same)
    for (j in unique(S) ){ #Theta update S can go between 1 and SMax but won't contain most of them
      #Only need to update those j for which s_i=j for some i
      for (r in c(1:R) ){
        ThList[[l]][j,r] <- MDPThSample(j,r,S,X,Y,cmu,cvar)
      }
    }
    #Step 2biii
    VList[[l]] <- matrix(0,nrow=SMax,ncol=R) #Preallocation for S update
    for (j in unique(S) ){ #When j is one of the elements of S, j=s_i for some i in product so impacts dist
      #Maybe I want j in c(1:max(S)) instead?
      for (r in c(1:R) ){
        VList[[l]][j,r] <- MDPVSample(j,r,U,V,X,M)
      }
    }
    #Step 2biv
    for (i in c(1:N) ){ #For updating slices
      S[i] <- MDPSSample(i,U,VList[[l]],X,Y,ThList[[l]],SMax)
    }
  }
  
  return( list("Thetas"=ThList,"StickBreaks"=VList) ) #Then draw is sum( W_i*f(.|theta_i) )
}

MDPUSample <- function(i,V,X,S){ #Update slicing ( Step 2b(i) ) This is just a unif(0,Trunc) for specified trunc point
  if (S[i]==1){
    Trunc <- V[1]
  } else {
  Trunc <- prod( (1-V)[1:(S[i]-1),X[i] ] )*V[S[i],X[i]]
  }
  U <- runif(1,0,Trunc) #Draws from uniform on 0 to W[s_i,X_i]
  return(U)
}

MDPThSample <- function(j,r,S,X,Y){ #Update thetas ( Step 2b(ii) )
#Use conjugacy of base measure alpha
  
}

MDPVSample <- function(j,r,U,V,X,M){#Update relative stick weights ( 2b(iii) )
  #Use quantile method for sampling from continuous distributions
  #Set endpoints as per pg 109 of vdV
  a <- max( U[ S==j&&X==r ] ) #truncation lower bound 
  b <- min(1 - U[X==r]/(prod( (1-V)[1:(S[X==r]-1),r] )*V[S[X==r],r]) ) #truncation upper bound
  Vdraw <- truncdist::rtrunc(1,"beta",a,b,shape1=1,shape2=M)
  return(Vdraw)
}

MDPSSample <- function(i,U,V,X,Y,Th,SMax){#Update pointer variables ( 2b(iv) )
  #These updates are from discrete distribution of "likelihood mass" supported on SMin to infty
  #Where SMin is the smallest index s of W for which W_{s} > U_i. Then we take s_i in SMin to infty (in practice SMax)
  l <- 1
  W <- V[l,X[i]]
  while (W <= U[i] ){ #Dependence on r surpressed
    l <- l+1
    W <- (W/V[(l-1),X[i]])*(1-V[(l-1),X[i]])*V[ l,X[i] ]
  }
  SMin <- l #Tells us the lowest W entry we can accept given that W_{s_i}>U_i
  prob <- exp( ( Y[i]-Th[SMin:SMax,X[i]] ) ) #We can point anywhere compatible with the indicator condition
  s <- sample(c(SMin:SMax),1,replace=TRUE,prob=prob) #Samples SMin to truncation SMax according to density
  return(s)
}

#NEED TO CHANGE THIS TO MAKE IT A DIFFERENT MIXTURE FOR EACH HMM STATE
MDPXSample <- function(eps,Y,Q,V,S,Th){ #tol, data, Qmat, V betas, S pointers, Th locations
  sig2 <- 1 #CHANGE TO UPDATE WITH INVERSE GAMMA
  W <- rep(0,length(V))
  i<-1
  W[1]=V[1]
  while ( sum(W)<=(1-eps) ){
    i <- i+1
    W[i] <- (W[i-1]/V[i-1])*(1-V[i-1])*V[i]
  }
  W <- W[W>0]/sum(W) #Takes the first 1-eps and renormalises
  M <- length(W) #Number of mixture components
  pi <- abs(eigen(t(Q))$vectors[,1])/sum(abs(eigen(t(Q))$vectors[,1])) #stationary dist
  dist <- distributionSet(dis="MIXTURE",mean=Th[1:M],var = rep(sig2,M), W ) #Argh all wrong! Correct! Only specifying one states emissions!
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
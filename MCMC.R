HmmMCMC <- function(M,N,mu,sig2,b,I){ #vector of bins, data, emissions (mu, sig2 length R*length(M))
R = 2 #states
Q <- t(matrix(c(0.7,0.3,0.2,0.8),2,2))
E <- length(M) #How many experiments
ListInputs <- vector("list",E)
ListData <- vector("list",E)
ListOutputs <- vector("list",E)
for (e in c(1:E) ){
  means <- mu[(R*e-1):(R*e)]
  vars <- sig2[(R*e-1):(R*e)]
  Link <- MyLinkAB( min(means) - 2*sqrt(sig2 [ which(means==min(means))[1] ] ) , max(means) + 2*sqrt( sig2[ which( means==max(means) )[1] ] ) )
  ListInputs[[e]] <- list("means"=means,"vars"=vars, "SampleSize" = N[e], "Link"=Link, "BinCount"=M[e], "BurnIn" = b[e], "Iterations" = I[e])
  ListData[[e]] <- SimulateHMMNorm( Q,means,vars,N[e] )
  Y <- ListData[[e]]$obs
  #X0 <- ListData[[e]]$states #Commented because can't see why this needs to be there.
  YBin <- factor( Bin(Y,M[e],Link) , c(1:M[e]) )
  ListOutputs[[e]] <- QWPosterior( YBin,R,M[e],b[e],I[e] ) #Gives simulated Q, W, states and loglikelihoods
} #can change last line to use QWPosteriorNoLatent if memory issues
return(list("Inputs"=ListInputs,"Data"=ListData,"Outputs"=ListOutputs))
}

HmmMCMC2 <- function(M,N,D,mu,sig2,I,Y=NULL,X0=NULL){ #D distinct data sets, N number of obs. No burn in. Uniform prior.
  R = 2 #states
  Q <- t(matrix(c(0.7,0.3,0.2,0.8),2,2))
  B <- length(M) #How many different bin sizes
  S <- length(N) #How many different sample sizes
  E <- S*B*D
  ListInputs <- vector("list",E)
  ListData <- vector("list",E)
  ListOutputs <- vector("list",E)
  e <- 1 #Initialise
  for (d in c(1:D) ){
    means <- mu[(R*d-1):(R*d)]
    vars <- sig2[(R*d-1):(R*d)]
    Link <- MyLinkAB( min(means) - 2*sqrt(sig2 [ which(means==min(means))[1] ] ) , max(means) + 2*sqrt( sig2[ which( means==max(means) )[1] ] ) )
    if (is.null(Y)){#can instead specify data set if don't want to generate new one
      Sims <- SimulateHMMNorm( Q,means,vars,max(N) ) #Simulates a full data set for the highest N valu
      Y <- Sims$obs
    }
    if (is.null(X0)){
      X0 <- Sims$states
    }
    for (j in c(1:S)){
      for (i in c(1:B)){
        ListInputs[[e]] <- list("means"=means,"vars"=vars, "Link"=Link, "SampleSize" = N[j], "BinCount"=M[i], "Iterations" = I)
        ListData[[e]] <- list("obs"=Y[1:N[j]],"states"=X0[1:N[j]]) #Observations and states for this experiment
        YBin <- factor( Bin(Y[1:N[j]],M[i],Link) , c(1:M[i]) ) #Bins first N[j] observations into M[i] bins
        ListOutputs[[e]] <- QWPosteriorSomeLatent( YBin,R,M[i],0,I) #Gives simulated Q, W, states and loglikelihoods
        e <- e+1 #Move to next entry of lists for following iteration of loop.
      }
    }
  }
  return(list("Inputs"=ListInputs,"Data"=ListData,"Outputs"=ListOutputs))
}

MCMCEmission <- function(Data,b=500,s=10){ #Does the pi2 posterior for the MCMC draws
  R <- 2
  #L <- length(Data$Outputs)
  L <- 1
  A <- vector("list",L)
  for (l in c(1:L)){
    Y <- Data$Inputs[[l]]$Link( Sims$Data[[l]]$obs ) #Transform data to [0,1] (using same link as before)
    QList <- Data$Outputs[[l]]$QList[seq(b,length(Data$Outputs[[l]]$QList),s)] #Provides thinned list
    A[[l]] <- EmissionPosterior(Y,R,QList)
  }
  return(A)
}

MCMCpi2 <- function(Data,I){ #Does the pi2 posterior for the MCMC draws (could make funct of Y)
  R <- 2
  #L <- length(Data$Outputs)
  L <- 1
  A <- vector("list",L)
  for (l in c(1:L)){
    Y <- Data$Inputs[[l]]$Link( Sims$Data[[l]]$obs ) #Transform data to [0,1] (using same link as before)
    A[[l]] <- FullPi2(Y,R,I)
  }
  return(A)
}

UnstoreLatent <- function(Data){
  L <- length(Data$Outputs)
  for (E in c(1:L) ){
    Data$Outputs[[E]]$XList <- NULL
  }
  return(Data)
}
  
MCMCPi1PlotsQ <- function(Data,b,s,Q=NULL,plotrow=4,plotcol=2,indexSET=NULL,both=FALSE,truth=TRUE,conf=TRUE,postmu=TRUE,basemu=TRUE){ #Data frame e.g. ExperimentsN500 , N1000 etc. b burn in vector. Q true
  if (is.null(indexSET)){
    indexSET = c(1:length(Data$Outputs))
    L = length(indexSET)
  }
  if (is.null(Q)){
    Q <- 0*Data$Outputs[[L]]$QList[[1]] #Gives 0 matrix of correct size
  }
  b <- rep(b,L/length(b))
  ExperimentsQThin <- vector("list",L)
  M <- vector("list",L)
  V <- vector("list",L)
  SD <- vector("list",L)
  par(mfrow=c(plotrow,plotcol))
  baseE=2
  for (E in indexSET ){
    QList <- Data$Outputs[[E]]$QList
    WList <- Data$Outputs[[E]]$WList
    LLHList <- Data$Outputs[[E]]$LLHList
    QList <- QList[(b[E]+1):30000]
    WList <- WList[(b[E]+1):30000]
    LLHList <- LLHList[(b[E]+1):30000]
    ExperimentsQThin[[E]] <- LabelSwapLLH(QList,WList,LLHList,s)$QThin
    R <- nrow(QList[[1]])
    M[[E]] <- PostMean(ExperimentsQThin[[E]])
    V[[E]] <- matrix(0,R,R)
    for (i in c(1:R)){
      for (j in c(1:R) ){
        V[[E]][i,j] <- var( EntryDraws( ExperimentsQThin[[E]],i,j ) )
      }
    }
    SD[[E]] = sqrt(V[[E]])
    Perms <- permutations(R,R)
    PermutedTrueQ <- Q
    if (sum(Q)>0){ #Only do if the Q was not NULL
    D <- rep(0,dim(Perms)[1]) #allocate
    for (j in c(1:dim(Perms)[1]) ){
      for (r in c(1:R)){
        for (s in c(1:R)){
          PermutedTrueQ[r,s] <- Q[Perms[j,r],Perms[j,s]] #permuting truth equivalent to permuting mean
        }
      }
      D[j] <- Distance(M[[E]],PermutedTrueQ)
    }
    PermIndex <- which(D==min(D))[1]
    for (r in c(1:R)){
      for (s in c(1:R)){
        PermutedTrueQ[r,s] <- Q[ Perms[PermIndex,r],Perms[PermIndex,s] ] #see which version of truth it fits best
      }
    }
    }
    j=1
    if (PermutedTrueQ[2,2]==0.8){
      j=2
    }
    if (E%%7==1){
      baseE = E+1
      
      frame()
      mtext(paste("Posterior draws for Q_{22} for ",Data$Inputs[[E]]$SampleSize," samples"),side = 3, line = -1.5)
    }
    if (E%%7==2){
      basej = j
    }
    #hist(EntryDraws(ExperimentsQThin[[E]],1,1),breaks=seq(0,1,0.0125),main = paste("Q11 Histogram in Experiment #",E,"of", L,"N=",Data$Inputs[[E]]$SampleSize, "Bins=",Data$Inputs[[E]]$BinCount),xlab = "Q(1,1)")
    #hist(EntryDraws(ExperimentsQThin[[E]],j,j),breaks=seq(0,1,0.005),main = paste("Sample size of ",Data$Inputs[[E]]$SampleSize, " with ",Data$Inputs[[E]]$BinCount," bins"),xlab = "First diagonal element")
    hist(EntryDraws(ExperimentsQThin[[E]],j,j),breaks=seq(0,1,0.005),main = paste(Data$Inputs[[E]]$BinCount," bins"),xlab = "Q_{22} value")
    if (postmu){ #Posterior mean for this experiment
      lines(c(M[[E]][j,j],M[[E]][j,j]),c(0,10000),col="red")
    }
    if (truth){ #True Q
      lines(c(PermutedTrueQ[j,j],PermutedTrueQ[j,j]),c(0,10000),col="red",lty=2)
    }
    if (conf){ #95% confidence bounds
      #print(quantile(EntryDraws(ExperimentsQThin[[E]],1,1),0.975))
      lines(c(quantile(EntryDraws(ExperimentsQThin[[E]],j,j),0.05),quantile(EntryDraws(ExperimentsQThin[[E]],j,j),0.05)),c(0,10000),col="blue")
      lines(c(quantile(EntryDraws(ExperimentsQThin[[E]],j,j),0.95),quantile(EntryDraws(ExperimentsQThin[[E]],j,j),0.95)),c(0,10000),col="blue")
    }
    if (basemu){ #Baseline posterior mean
      if (E>baseE){
        print(c(M[[baseE]][basej,basej],M[[baseE]][basej,basej]))
        lines(c(M[[baseE]][basej,basej],M[[baseE]][basej,basej]),c(0,10000),col="green")
      }
      
    }
    
    if (both){
      j = 3 - j
      #hist(EntryDraws(ExperimentsQThin[[E]],1,1),breaks=seq(0,1,0.0125),main = paste("Q11 Histogram in Experiment #",E,"of", L,"N=",Data$Inputs[[E]]$SampleSize, "Bins=",Data$Inputs[[E]]$BinCount),xlab = "Q(1,1)")
      hist(EntryDraws(ExperimentsQThin[[E]],j,j),breaks=seq(0,1,0.005),main = paste("Sample size of ",Data$Inputs[[E]]$SampleSize, " with ",Data$Inputs[[E]]$BinCount," bins"),xlab = "Second diagonal element")
      if (postmu){ #Posterior mean for this experiment
        lines(c(M[[E]][j,j],M[[E]][j,j]),c(0,10000),col="blue")
      }
      if (truth){ #True Q
        lines(c(PermutedTrueQ[j,j],PermutedTrueQ[j,j]),c(0,10000),col="red")
      }
      if (conf){ #95% confidence bounds
        #print(quantile(EntryDraws(ExperimentsQThin[[E]],1,1),0.975))
        lines(c(quantile(EntryDraws(ExperimentsQThin[[E]],j,j),0.025),quantile(EntryDraws(ExperimentsQThin[[E]],j,j),0.025)),c(0,10000),col="aquamarine")
        lines(c(quantile(EntryDraws(ExperimentsQThin[[E]],j,j),0.975),quantile(EntryDraws(ExperimentsQThin[[E]],j,j),0.975)),c(0,10000),col="aquamarine")
      }
      if (basemu){ #Baseline posterior mean
        if (E>baseE){
          print(c(M[[baseE]][3-basej,3-basej],M[[baseE]][3-basej,3-basej]))
          lines(c(M[[baseE]][3-basej,3-basej],M[[baseE]][3-basej,3-basej]),c(0,10000),col="green")
        }
        
      }
    }
  }
  return(list("QThin"=ExperimentsQThin,"PostMean"=M,"PostVariance"=V, "PostSD"=SD))
}

MCMCPi1PlotsQW <- function(Data,b,s,Q=NULL){ #Data frame e.g. ExperimentsN500 , N1000 etc. b burn in vector. Q true
  L <- length(Data$Outputs)
  if (is.null(Q)){
    Q <- 0*Data$Outputs[[L]]$QList[[1]] #Gives 0 matrix of correct size
  }
  b <- rep(b,L/length(b))
  ExperimentsQThin <- vector("list",L)
  ExperimentsWThin <- vector("list",L)
  MQ <- vector("list",L)
  VQ <- vector("list",L)
  SDQ <- vector("list",L)
  MW <- vector("list",L)
  VW <- vector("list",L)
  SDW <- vector("list",L)
  par(mfrow=c(2,4))
  for (E in c(1:L) ){
    QList <- Data$Outputs[[E]]$QList
    WList <- Data$Outputs[[E]]$WList
    LLHList <- Data$Outputs[[E]]$LLHList
    QList <- QList[(b[E]+1):30000]
    WList <- WList[(b[E]+1):30000]
    LLHList <- LLHList[(b[E]+1):30000]
    ExperimentsQThin[[E]] <- LabelSwapLLH(QList,WList,LLHList,s)$QThin
    ExperimentsWThin[[E]] <- LabelSwapLLH(QList,WList,LLHList,s)$WThin
    R <- nrow(QList[[1]])
    M <- ncol(WList[[1]])
    MQ[[E]] <- PostMean(ExperimentsQThin[[E]])
    VQ[[E]] <- matrix(0,R,R)
    for (i in c(1:R)){
      for (j in c(1:R) ){
        VQ[[E]][i,j] <- var( EntryDraws( ExperimentsQThin[[E]],i,j ) )
      }
    }
    SDQ[[E]] = sqrt(VQ[[E]])
    
    MW[[E]] <- PostMean(ExperimentsWThin[[E]])
    #VW[[E]] <- matrix(0,R,M)
    #for (i in c(1:R)){
     # for (j in c(1:M) ){
     #   VW[[E]][i,j] <- var( EntryDraws( ExperimentsWThin[[E]],i,j ) )
    #  }
    #}
    #SDW[[E]] = sqrt(VW[[E]])
    
    Perms <- permutations(R,R)
    PermutedTrueQ <- Q
    if (sum(Q)>0){ #Only do if the Q was not NULL
      D <- rep(0,dim(Perms)[1]) #allocate
      for (j in c(1:dim(Perms)[1]) ){
        for (r in c(1:R)){
          for (s in c(1:R)){
            PermutedTrueQ[r,s] <- Q[Perms[j,r],Perms[j,s]] #permuting truth equivalent to permuting mean
          }
        }
        D[j] <- Distance(MQ[[E]],PermutedTrueQ) #SHOULD BRING IN DISTANCE FOR W ALSO
      }
      PermIndex <- which(D==min(D))[1]
      for (r in c(1:R)){
        for (s in c(1:R)){
          PermutedTrueQ[r,s] <- Q[ Perms[PermIndex,r],Perms[PermIndex,s] ] #see which version of truth it fits best
        }
      }
    }
    #hist(EntryDraws(ExperimentsQThin[[E]],1,1),breaks=seq(0,1,0.0125),main = paste("Q11 Histogram in Experiment #",E,"of", L,"N=",Data$Inputs[[E]]$SampleSize, "Bins=",Data$Inputs[[E]]$BinCount),xlab = "Q(1,1)")
    hist(EntryDraws(ExperimentsQThin[[E]],1,1),breaks=seq(0,1,0.0125),main = paste("Sample size of ",Data$Inputs[[E]]$SampleSize, " with ",Data$Inputs[[E]]$BinCount," bins"),xlab = "First diagonal element")
    lines(c(MQ[[E]][1,1],MQ[[E]][1,1]),c(0,10000),col="blue")
    lines(c(PermutedTrueQ[1,1],PermutedTrueQ[1,1]),c(0,10000),col="red")
    #hist(EntryDraws(ExperimentsQThin[[E]],2,2),breaks=seq(0,1,0.0125),main = paste("Q22 Histogram in Experiment #",E,"of", L,"N=",Data$Inputs[[E]]$SampleSize, "Bins=",Data$Inputs[[E]]$BinCount),xlab = "Q(2,2)")
    hist(EntryDraws(ExperimentsQThin[[E]],2,2),breaks=seq(0,1,0.0125),xlab = "Second diagonal element")
    lines(c(MQ[[E]][2,2],MQ[[E]][2,2]),c(0,10000),col="blue")
    lines(c(PermutedTrueQ[2,2],PermutedTrueQ[2,2]),c(0,10000),col="red")
    Link <- Data$Inputs[[E]]$Link #Retrieves the link function
    Line <- seq(-5,5,0.01) #Defines the line space
    Hist <- WHist(MW[[E]][1,],Link)(Line) #Defines the (transformed) hist to be plotted
    mu <- 2*(Perms[PermIndex,1]-1.5) #Decides which is the correct mu
    plot( Line,(1/sqrt(2*pi))*exp(-(Line-mu)^2/2),"l",col="red",main = paste("Emission 1, with N=",Data$Inputs[[E]]$SampleSize, "Bins=",Data$Inputs[[E]]$BinCount),xlab = "x",ylab="Density True/PostMean" )
    lines(Line,Hist,"l",col="blue")
    Hist <- WHist(MW[[E]][2,],Link)(Line)
    mu <- 2*(Perms[PermIndex,2]-1.5)
    plot( Line,(1/sqrt(2*pi))*exp(-(Line-mu)^2/2) ,"l",col="red",main = paste("Emission 2, with N=",Data$Inputs[[E]]$SampleSize, "Bins=",Data$Inputs[[E]]$BinCount),xlab = "x",ylab="Density True/PostMean" )
    lines(Line,Hist,"l",col="blue")
    }
  return( list("QThin"=ExperimentsQThin,"PostMeanQ"=MQ,"PostVarianceQ"=VQ, "PostSDQ"=SDQ,"WThin"=ExperimentsWThin,"PostMeanW"=MW,"PostVarianceW"=VW, "PostSDW"=SDW) )
}

MCMCPi1PlotsW <- function(Data,b,s,Q=NULL){ #Data frame e.g. ExperimentsN500 , N1000 etc. b burn in vector. Q true
  L <- length(Data$Outputs)
  if (is.null(Q)){
    Q <- 0*Data$Outputs[[L]]$QList[[1]] #Gives 0 matrix of correct size
  }
  b <- rep(b,L/length(b))
  ExperimentsQThin <- vector("list",L)
  ExperimentsWThin <- vector("list",L)
  MQ <- vector("list",L)
  VQ <- vector("list",L)
  SDQ <- vector("list",L)
  MW <- vector("list",L)
  VW <- vector("list",L)
  SDW <- vector("list",L)
  par(mfrow=c(1,2))
  for (E in c(1:L) ){
    QList <- Data$Outputs[[E]]$QList
    WList <- Data$Outputs[[E]]$WList
    LLHList <- Data$Outputs[[E]]$LLHList
    QList <- QList[(b[E]+1):30000]
    WList <- WList[(b[E]+1):30000]
    LLHList <- LLHList[(b[E]+1):30000]
    ExperimentsQThin[[E]] <- LabelSwapLLH(QList,WList,LLHList,s)$QThin
    ExperimentsWThin[[E]] <- LabelSwapLLH(QList,WList,LLHList,s)$WThin
    R <- nrow(QList[[1]])
    M <- ncol(WList[[1]])
    MQ[[E]] <- PostMean(ExperimentsQThin[[E]])
    VQ[[E]] <- matrix(0,R,R)
    for (i in c(1:R)){
      for (j in c(1:R) ){
        VQ[[E]][i,j] <- var( EntryDraws( ExperimentsQThin[[E]],i,j ) )
      }
    }
    SDQ[[E]] = sqrt(VQ[[E]])
    
    MW[[E]] <- PostMean(ExperimentsWThin[[E]])
    #VW[[E]] <- matrix(0,R,M)
    #for (i in c(1:R)){
    # for (j in c(1:M) ){
    #   VW[[E]][i,j] <- var( EntryDraws( ExperimentsWThin[[E]],i,j ) )
    #  }
    #}
    #SDW[[E]] = sqrt(VW[[E]])
    
    Perms <- permutations(R,R)
    PermutedTrueQ <- Q
    if (sum(Q)>0){ #Only do if the Q was not NULL
      D <- rep(0,dim(Perms)[1]) #allocate
      for (j in c(1:dim(Perms)[1]) ){
        for (r in c(1:R)){
          for (s in c(1:R)){
            PermutedTrueQ[r,s] <- Q[Perms[j,r],Perms[j,s]] #permuting truth equivalent to permuting mean
          }
        }
        D[j] <- Distance(MQ[[E]],PermutedTrueQ) #SHOULD BRING IN DISTANCE FOR W ALSO
      }
      PermIndex <- which(D==min(D))[1]
      for (r in c(1:R)){
        for (s in c(1:R)){
          PermutedTrueQ[r,s] <- Q[ Perms[PermIndex,r],Perms[PermIndex,s] ] #see which version of truth it fits best
        }
      }
    }
    #hist(EntryDraws(ExperimentsQThin[[E]],1,1),breaks=seq(0,1,0.0125),main = paste("Q11 Histogram in Experiment #",E,"of", L,"N=",Data$Inputs[[E]]$SampleSize, "Bins=",Data$Inputs[[E]]$BinCount),xlab = "Q(1,1)")
    #hist(EntryDraws(ExperimentsQThin[[E]],1,1),breaks=seq(0,1,0.0125),main = paste("Q11 Histogram for N=",Data$Inputs[[E]]$SampleSize, "Bins=",Data$Inputs[[E]]$BinCount),xlab = "Q(1,1)")
    #lines(c(MQ[[E]][1,1],MQ[[E]][1,1]),c(0,10000),col="blue")
    #lines(c(PermutedTrueQ[1,1],PermutedTrueQ[1,1]),c(0,10000),col="red")
    #hist(EntryDraws(ExperimentsQThin[[E]],2,2),breaks=seq(0,1,0.0125),main = paste("Q22 Histogram in Experiment #",E,"of", L,"N=",Data$Inputs[[E]]$SampleSize, "Bins=",Data$Inputs[[E]]$BinCount),xlab = "Q(2,2)")
    #hist(EntryDraws(ExperimentsQThin[[E]],2,2),breaks=seq(0,1,0.0125),main = paste("Q22 Histogram for N=",Data$Inputs[[E]]$SampleSize, "Bins=",Data$Inputs[[E]]$BinCount),xlab = "Q(2,2)")
    #lines(c(MQ[[E]][2,2],MQ[[E]][2,2]),c(0,10000),col="blue")
    #lines(c(PermutedTrueQ[2,2],PermutedTrueQ[2,2]),c(0,10000),col="red")
    Link <- Data$Inputs[[E]]$Link #Retrieves the link function
    Line <- seq(-5,5,0.01) #Defines the line space
    Hist <- WHist(MW[[E]][1,],Link)(Line) #Defines the (transformed) hist to be plotted
    mu <- 2*(Perms[PermIndex,1]-1.5) #Decides which is the correct mu
    plot( Line,(1/sqrt(2*pi))*exp(-(Line-mu)^2/2),"l",col="red",main = paste("Emission 1, with N=",Data$Inputs[[E]]$SampleSize, "Bins=",Data$Inputs[[E]]$BinCount),xlab = "x",ylab="Density True/PostMean" )
    lines(Line,Hist,"l",col="blue")
    Hist <- WHist(MW[[E]][2,],Link)(Line)
    mu <- 2*(Perms[PermIndex,2]-1.5)
    plot( Line,(1/sqrt(2*pi))*exp(-(Line-mu)^2/2) ,"l",col="red",main = paste("Emission 2, with N=",Data$Inputs[[E]]$SampleSize, "Bins=",Data$Inputs[[E]]$BinCount),xlab = "x",ylab="Density True/PostMean" )
    lines(Line,Hist,"l",col="blue")
  }
  return( list("QThin"=ExperimentsQThin,"PostMeanQ"=MQ,"PostVarianceQ"=VQ, "PostSDQ"=SDQ,"WThin"=ExperimentsWThin,"PostMeanW"=MW,"PostVarianceW"=VW, "PostSDW"=SDW) )
}

WHist <- function(W,Link=NULL){#Plots histogram with weights W (W is 1 by M)
  
  LinkDeriv <- function(x){ #For evluating the derivative of the link function we use
    if (x >= -3 && x<=3){
      Out <- ( inv.logit(c(-3,3))[2]-inv.logit(c(-3,3))[1] )/6 #In the piecewise linear part
      #THIS IS ONLY THE CORRECT DERIVATIVE FOR MY CASE (linear -3 to +3)
    }else{
      Out <- (exp(x))/(1+exp(x)^2)
    }
    return(Out)
  }
  
  A <- function(x){ #Defines the function of interest to return (outputs the density corresponding to W)
    M <- length(W)
    if (is.null(Link)==FALSE){
      xLink <- Link(x)
      xLinkderiv <- (Vectorize(LinkDeriv))(x) #computes the derivative at each point on the x axis
    }
    for (m in c(1:M)){
      if ( (m-1)/M <= xLink && xLink < (m/M) ){ #so x is in bin m
        C <- W[m]*M*xLinkderiv
      }
    }
    if (xLink==1){
      C <- 0
    }
    return(C)
  }
  B <- Vectorize(A)
  return(B)
}

MCMCPi2PlotsQW <- function(Data,b,s,Q=NULL){ #Data frame e.g. ExperimentsN500 , N1000 etc. b burn in vector. Q true
  L <- length(Data)
  if (is.null(Q)){
    Q <- 0*Data[[L]]$QList[[1]] #Gives 0 matrix of correct size
  }
  b <- rep(b,L/length(b))
  ExperimentsQThin <- vector("list",L)
  ExperimentsWThin <- vector("list",L)
  MQ <- vector("list",L)
  VQ <- vector("list",L)
  SDQ <- vector("list",L)
  MW <- vector("list",L)
  VW <- vector("list",L)
  SDW <- vector("list",L)
  par(mfrow=c(2,5))
  for (E in c(1:L) ){
    QList <- Data[[E]]$QList
    WList <- Data[[E]]$WList
    LLHList <- Data[[E]]$LLHList
    QList <- QList[(b[E]+1):30000]
    WList <- WList[(b[E]+1):30000]
    LLHList <- LLHList[(b[E]+1):30000]
    ExperimentsQThin[[E]] <- LabelSwapLLH(QList,WList,LLHList,s)$QThin
    ExperimentsWThin[[E]] <- LabelSwapLLH(QList,WList,LLHList,s)$WThin
    R <- nrow(QList[[1]])
    M <- ncol(WList[[1]])
    MQ[[E]] <- PostMean(ExperimentsQThin[[E]])
    VQ[[E]] <- matrix(0,R,R)
    for (i in c(1:R)){
      for (j in c(1:R) ){
        VQ[[E]][i,j] <- var( EntryDraws( ExperimentsQThin[[E]],i,j ) )
      }
    }
    SDQ[[E]] = sqrt(VQ[[E]])
    
    MW[[E]] <- PostMean(ExperimentsWThin[[E]])
    VW[[E]] <- matrix(0,R,M)
    for (i in c(1:R)){
      for (j in c(1:M) ){
        VW[[E]][i,j] <- var( EntryDraws( ExperimentsWThin[[E]],i,j ) )
      }
    }
    SDW[[E]] = sqrt(VW[[E]])
    
    Perms <- permutations(R,R)
    PermutedTrueQ <- Q
    if (sum(Q)>0){ #Only do if the Q was not NULL
      D <- rep(0,dim(Perms)[1]) #allocate
      for (j in c(1:dim(Perms)[1]) ){
        for (r in c(1:R)){
          for (s in c(1:R)){
            PermutedTrueQ[r,s] <- Q[Perms[j,r],Perms[j,s]] #permuting truth equivalent to permuting mean
          }
        }
        D[j] <- Distance(MQ[[E]],PermutedTrueQ)
      }
      PermIndex <- which(D==min(D))[1]
      for (r in c(1:R)){
        for (s in c(1:R)){
          PermutedTrueQ[r,s] <- Q[ Perms[PermIndex,r],Perms[PermIndex,s] ] #see which version of truth it fits best
        }
      }
    }
    #hist(EntryDraws(ExperimentsQThin[[E]],1,1),breaks=seq(0,1,0.0125),main = paste("Q11 Histogram in Experiment #",E,"of", L,"N=",Data$Inputs[[E]]$SampleSize, "Bins=",Data$Inputs[[E]]$BinCount),xlab = "Q(1,1)")
    hist(EntryDraws(ExperimentsQThin[[E]],1,1),breaks=seq(0,1,0.0125),main = paste("Q11 Histogram for N=",Data[[E]]$SampleSize), xlab = "Q(1,1)")
    lines(c(MQ[[E]][1,1],MQ[[E]][1,1]),c(0,10000),col="blue")
    lines(c(PermutedTrueQ[1,1],PermutedTrueQ[1,1]),c(0,10000),col="red")
    #hist(EntryDraws(ExperimentsQThin[[E]],2,2),breaks=seq(0,1,0.0125),main = paste("Q22 Histogram in Experiment #",E,"of", L,"N=",Data$Inputs[[E]]$SampleSize, "Bins=",Data$Inputs[[E]]$BinCount),xlab = "Q(2,2)")
    hist(EntryDraws(ExperimentsQThin[[E]],2,2),breaks=seq(0,1,0.0125),main = paste("Q22 Histogram for N=",Data[[E]]$SampleSize), xlab = "Q(2,2)")
    lines(c(MQ[[E]][2,2],MQ[[E]][2,2]),c(0,10000),col="blue")
    lines(c(PermutedTrueQ[2,2],PermutedTrueQ[2,2]),c(0,10000),col="red")
    Link <- Data$Inputs[[E]]$Link
    Line <- seq(-5,5,0.01)
    Hist <- WHist(MW[[E]][1,],Link)(Line)
    plot(Line,Hist,"l",col="green",main = paste("Emission 1, with N=",Data[[E]]$SampleSize),xlab = "x")
    Hist <- WHist(MW[[E]][2,],Link)(Line)
    plot(Line,Hist,"l",col="green",main = paste("Emission 2, with N=",Data[[E]]$SampleSize),xlab = "x")
  }
  return( list("QThin"=ExperimentsQThin,"PostMeanQ"=MQ,"PostVarianceQ"=VQ, "PostSDQ"=SDQ,"WThin"=ExperimentsWThin,"PostMeanW"=MW,"PostVarianceW"=VW, "PostSDW"=SDW) )
}
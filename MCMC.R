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


HmmMCMC2 <- function(M,N,D,mu,sig2,I){ #D distinct data sets, N number of obs. No burn in. Uniform prior.
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
    Sims <- SimulateHMMNorm( Q,means,vars,max(N) ) #Simulates a full data set for the highest N value
    Y <- Sims$obs
    X0 <- Sims$states
    for (j in c(1:S)){
      for (i in c(1:B)){
        ListInputs[[e]] <- list("means"=means,"vars"=vars, "Link"=Link, "SampleSize" = N[j], "BinCount"=M[i], "Iterations" = I)
        ListData[[e]] <- list("obs"=Y[1:N[j]],"states"=X0[1:N[j]]) #Observations and states for this experiment
        YBin <- factor( Bin(Y[1:N[j]],M[i],Link) , c(1:M[i]) ) #Bins first N[j] observations into M[i] bins
        ListOutputs[[e]] <- QWPosteriorNoLatent( YBin,R,M[i],0,I) #Gives simulated Q, W, states and loglikelihoods
        e <- e+1 #Move to next entry of lists for following iteration of loop.
      }
    }
  }
  return(list("Inputs"=ListInputs,"Data"=ListData,"Outputs"=ListOutputs))
}



UnstoreLatent <- function(Data){
  L <- length(Data$Outputs)
  for (E in c(1:L) ){
    Data$Outputs[[E]]$XList <- NULL
  }
  return(Data)
}
  


MCMCPlots <- function(Data,b,s,Q=NULL){ #Data frame e.g. ExperimentsN500 , N1000 etc. b burn in vector. Q true
  L <- length(Data$Outputs)
  if (is.null(Q)){
    Q <- 0*Data$Outputs[[L]]$QList[[1]] #Gives 0 matrix of correct size
  }
  b <- rep(b,L/length(b))
  ExperimentsQThin <- vector("list",L)
  M <- vector("list",L)
  V <- vector("list",L)
  SD <- vector("list",L)
  par(mfrow=c(2,5))
  for (E in c(1:L) ){
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
    #hist(EntryDraws(ExperimentsQThin[[E]],1,1),breaks=seq(0,1,0.0125),main = paste("Q11 Histogram in Experiment #",E,"of", L,"N=",Data$Inputs[[E]]$SampleSize, "Bins=",Data$Inputs[[E]]$BinCount),xlab = "Q(1,1)")
    hist(EntryDraws(ExperimentsQThin[[E]],1,1),breaks=seq(0,1,0.0125),main = paste("Q11 Histogram for N=",Data$Inputs[[E]]$SampleSize, "Bins=",Data$Inputs[[E]]$BinCount),xlab = "Q(1,1)")
    lines(c(M[[E]][1,1],M[[E]][1,1]),c(0,10000),col="blue")
    lines(c(PermutedTrueQ[1,1],PermutedTrueQ[1,1]),c(0,10000),col="red")
    #hist(EntryDraws(ExperimentsQThin[[E]],2,2),breaks=seq(0,1,0.0125),main = paste("Q22 Histogram in Experiment #",E,"of", L,"N=",Data$Inputs[[E]]$SampleSize, "Bins=",Data$Inputs[[E]]$BinCount),xlab = "Q(2,2)")
    hist(EntryDraws(ExperimentsQThin[[E]],2,2),breaks=seq(0,1,0.0125),main = paste("Q22 Histogram for N=",Data$Inputs[[E]]$SampleSize, "Bins=",Data$Inputs[[E]]$BinCount),xlab = "Q(2,2)")
    lines(c(M[[E]][2,2],M[[E]][2,2]),c(0,10000),col="blue")
    lines(c(PermutedTrueQ[2,2],PermutedTrueQ[2,2]),c(0,10000),col="red")
  }
  return(list("QThin"=ExperimentsQThin,"PostMean"=M,"PostVariance"=V, "PostSD"=SD))
}

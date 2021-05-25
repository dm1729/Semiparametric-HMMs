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
  Link <- MyLinkAB( min(means) - 2*sig2 [ which(means==min(means))[1] ] , max(means) + 2*sig2[ which( means==max(means) )[1] ] )
  ListInputs[[e]] <- list("means"=means,"vars"=vars, "SampleSize" = N[e], "Link"=Link, "BinCount"=M[e], "BurnIn" = b[e], "Iterations" = I[e])
  ListData[[e]] <- SimulateHMMNorm( Q,means,vars,N[e] )
  Y <- ListData[[e]]$obs
  X0 <- ListData[[e]]$states
  YBin <- factor( Bin(Y,M[e],Link) , c(1:M[e]) )
  ListOutputs[[e]] <- QWPosterior( YBin,R,M[e],b[e],I[e] ) #Gives simulated Q, W, states and loglikelihoods
} #can change last line to use QWPosteriorNoLatent if memory issues
return(list("Inputs"=ListInputs,"Data"=ListData,"Outputs"=ListOutputs))
}

MCMCPlots <- function(Data,b,s){ #Data frame e.g. ExperimentsN500 , N1000 etc. b burn in vector
  L <- length(Data$Outputs)
  b <- rep(b,L/length(b))
  ExperimentsQThin <- vector("list",L)
  M <- vector("list",L)
  V <- vector("list",L)
  SD <- vector("list",L)
  par(mfrow=c(2,3))
  for (E in c(1:L) ){
    QList <- Data$Outputs[[E]]$QList
    WList <- Data$Outputs[[E]]$WList
    LLHList <- Data$Outputs[[E]]$LLHList
    QList <- QList[(b[E]+1):30000]
    WList <- WList[(b[E]+1):30000]
    LLHList <- LLHList[(b[E]+1):30000]
    ExperimentsQThin[[E]] <- LabelSwapLLH(QList,WList,LLHList,s)$QThin
    hist(EntryDraws(ExperimentsQThin[[E]],1,1),breaks=seq(0,1,0.0125),main = paste("Posterior MCMC for Diagonal 1, in Experiment #",E,"of 9,N=",Data$Inputs[[E]]$SampleSize),xlab = "Q(1,1)")
    hist(EntryDraws(ExperimentsQThin[[E]],2,2),breaks=seq(0,1,0.0125),main = paste("Posterior MCMC for Diagonal 2, in Experiment #",E,"of 9,N=",Data$Inputs[[E]]$SampleSize),xlab = "Q(2,2)")
    R <- nrow(QList[[1]])
    M[[E]] <- PostMean(ExperimentsQThin[[E]])
    V[[E]] <- matrix(0,R,R)
    for (i in c(1:R)){
      for (j in c(1:R) ){
        V[[E]][i,j] <- var( EntryDraws( ExperimentsQThin[[E]],i,j ) )
      }
    }
    SD[[E]] = sqrt(V[[E]])
  }
  return(list("QThin"=ExperimentsQThin,"PostMean"=M,"PostVariance"=V, "PostSD"=SD))
}

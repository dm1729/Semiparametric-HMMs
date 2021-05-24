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
}
return(list("Inputs"=ListInputs,"Data"=ListData,"Outputs"=ListOutputs))
}
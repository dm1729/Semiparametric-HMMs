install.packages("RHmm", repos="http://R-Forge.R-project.org")
install.packages("gtools")
require(RHmm)
require(gtools)
source('~/Documents/Semiparametric HMMs/DirichletMixtures.R')
Inputs <- readRDS("Inputs.Rda")
Outputs <- vector("list",length = length(Inputs) )
for ( i in c(1:Inputs) ){
  Outputs[[i]] < - MDPTrunc(Inputs[[i]]$Y,Inputs[[i]]$M,Inputs[[i]]$cmu,Inputs[[i]]$cvar,Inputs[[i]]$igshape,Inputs[[i]]$igrate,Inputs[[i]]$QList,Inputs[[i]]$Xinit,SMax=NULL,Inputs[[i]]$C)
  saveRDS(Outputs[[i]],paste("Outputs",i,".Rda",recycle0 = TRUE))
}
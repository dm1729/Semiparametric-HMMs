library(MASS)
library(nlme)
library(gtools)
library(RHmm)
source('DirichletMixtures.R')
#source('RHmm.R')
#source('RHmm-internals.R')
Inputs <- readRDS("~/Documents/R Code/Semiparametric-HMMs/Inputs24112021.Rda")
Outputs <- vector("list",length = length(Inputs) )
for ( i in c(3,4,2,1) ){
  Outputs[[i]] <- MDPTrunc(Inputs[[i]]$Y,Inputs[[i]]$M,Inputs[[i]]$cmu,Inputs[[i]]$cvar,Inputs[[i]]$igshape,Inputs[[i]]$igrate,Inputs[[i]]$QList,Inputs[[i]]$Xinit,SMax=NULL,Inputs[[i]]$C)
  saveRDS(Outputs[[i]],paste("Pi2draws",i,".Rda",recycle0 = TRUE))
}

#1 N=2500: 8 bins
#2 N=5000: 8 bins
#3 N=10000: 16 bins
#4 N=10000: 4 bins

#N=1000 4 bins already done (in Outputs 1 .Rda)
#Diagnostics - Draws samples with true paramaters (to check Gibbs sampler code works)
IT <- 10000 #iterations
Qd <- vector("list",IT)
Wd <- vector("list",IT)
Xd <- vector("list",IT)
A <- priorset(R,M)[[1]]
B <- priorset(R,M)[[2]]
for (i in c(1:IT)){
  Qd[[i]] <- QGibbs(X0,A)
}
for (i in c(1:IT)){
  Wd[[i]] <- WGibbs(X0,Y,B)
}
for (i in c(1:IT) ){
  Xd[[i]] <- XSample(Y,Q,W)
}

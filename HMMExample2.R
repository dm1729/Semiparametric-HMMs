#HMM Example 2 uses mixture of normals to generate data then map
R = 2 #states
M = 8 #bins
Q <- t(matrix(c(0.6,0.4,0.2,0.8),2,2))
N <- 500 #Number of samples to draw
mu <- c(-2,2)
v <- c(0.5,0.5)
S <- SimulateHMMNorm(Q,mu,v,N)
Y <- S$obs
X0 <- S$states
YBin <- factor(Bin(Y,M,inv.logit),c(1:M)) #Needs to be factor type to work with multinom
b <- 0 #burn in
I <- 1000 #iterations
Adir <- 1
Bdir <- 1
        Out <- QWPosterior(YBin,R,M,b,I,Adir,Bdir)
#Out <- QWPosteriorFixState(YBin,R,M,b,I,X0) #Diagnostic to see problems with X simulations
QList <- Out$QList
WList <- Out$WList
XList <- Out$XList
#Q11 <- EntryDraws(QList,1,1)
#Q22 <- EntryDraws(QList,2,2)
#Qdata<-data.frame(x=Q11,y=Q22)
#ggplot(Qdata, aes(x=x, y=y) ) +
# geom_bin2d(bins = 40) +
#scale_fill_continuous(type = "viridis") +
#theme_bw()
# Adjust binning (interpolate - can be computationally intensive for large datasets)
#k <- kde2d(Qdata$x, Qdata$y, n=200)
#image(k)
s <- 20 #thinning parameter
Thin <- LabelSwap(QList,WList,YBin,s)
QThin <- Thin$QThin
WThin <- Thin$WThin
hist(EntryDraws(QThin,1,1),breaks=seq(0,1,0.0125))
hist(EntryDraws(QThin,2,2),breaks=seq(0,1,0.0125))

#HMM Example directly generates multinomial data for testing
R = 2
M = 3
Q <- t(matrix(c(0.3,0.7,0.2,0.8),2,2))
W <- matrix(0,R,M)
W[1,] <- rdirichlet(1,2*c((M):1)) #set up emissions
W[2,] <- rdirichlet(1,2*c(1:M))
N <- 10 #Number of samples to draw
S <- SimulateHMM(Q,W,N)
Y <- S$obs
X0 <- S$states
b <- 50 #burn in
I <- 200 #iterations
Out <- QWPosterior(Y,R,M,b,I)
#Out <- QWPosteriorFixState(Y,R,M,b,I,X0) #Diagnostic to see problems with X simulations
QList <- Out$QList
WList <- Out$WList
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
Thin <- LabelSwap(QList,WList,Y,s)
QThin <- Thin$QThin
#WThin <- Thin$WThin
hist(EntryDraws(QThin,1,1),breaks=seq(0,1,0.05))
hist(EntryDraws(QThin,2,2),seq(0,1,0.05))

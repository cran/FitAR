`LBQPlot` <-
function(obj, SquaredQ=FALSE){
if (SquaredQ){
    res<-residuals(obj)
    QQ<-LjungBoxTest(res, StartLag=1, SquaredQ=TRUE, lag.max=10)
}
else 
    QQ<-obj$LjungBoxQ
    
plot(QQ[,1],QQ[,3], xlab="lag", ylab="p-value", ylim=c(0,1), main="Ljung-Box Test")
abline(h=0.05, col="red", lty=2)
}


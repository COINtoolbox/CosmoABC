library(sm)
library(fields)


pdf("weights_all.pdf",width=12,height=8)

for (i in 1:6){
  
    round = as.character( i )
    data_param<-read.table(file=paste("SMC_ABC_SZ_", round,".dat", sep=""),header=T)
    wt<-read.table(file=paste("weights_",round,".dat",sep=""),header=F, nrows=1)

    y1<-data_param[,c("Om")]
    wt1<-as.numeric(wt[-1])

    d1 = density(y1,from = -5, to = 5,n = 2048,kernel="epanechnikov")  #<- estimate density of sample without weights; should look uniform
    d2 = density(y1,weights = wt1/sum(wt1),kernel="epanechnikov") # re-weight particle system to target Gaussian

    par(mfrow = c(1,2))
    hist(y1,prob = TRUE, main = paste("No weights - time.steps = ", round,sep=""), xlim=c(0.15,0.3))
    lines(d1$x, d1$y,col="red")
    plot(x=d2$x, y=d2$y, type="l",main = paste("With weights - time.steps = ", round,sep=""), xlim=c(0.15,0.3),col="red")
}
dev.off()

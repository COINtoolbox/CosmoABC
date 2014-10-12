
library(sm)
library(fields)


#  A two-dimensional example

pdf(paste("weights_all.pdf",sep=""),width=12,height=8)

for(i in 1:12){
    round = as.character(i)

    data_param<-read.table(file=paste("SMC_ABC_SZ_", round,".dat",sep=""),header=T)
    wt<-read.table(file=paste("weights_", round, ".dat", sep=""))

    y1<-data_param[,c("Om","sigma8")]
    wt1<-as.numeric(wt[,-1])

    
    par(mfrow = c(1,2))
    out = sm.density(y1, display = "image", xlim=c(0.1,.5),ylim=c(0.7,0.9))
    contour(out$eval.points[,1],out$eval.points[,2],out$estimate,add = TRUE)
    title({paste("No weights - time.steps = ",round)})
    out.wt = sm.density(y1, display = "image", weights = wt1/sum(wt1), nbins=0,xlim=c(0.1,.5),ylim=c(0.7,0.9))
    contour(out.wt$eval.points[,1],out.wt$eval.points[,2],out.wt$estimate,add = TRUE,main = "With weights")
    title({paste("With weights - time.steps = ",round,sep="")})
    
}
dev.off()
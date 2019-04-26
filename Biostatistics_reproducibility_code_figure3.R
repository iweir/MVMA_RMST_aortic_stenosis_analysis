library(metaRMST)
library(survRM2)

# read in built-in dataset 
data(AorticStenosisTrials)

palette2 <- c("#FF57A8FF", "#C729D6FF", "#5000FFFF", "#FF9C63FF", "#FFFF60FF")

J <- 5 #number of trials
t <- seq(0,40,by=0.25)
RMSTcurveRes <- matrix(NA, length(t), J+1)

for (j in 1:J){
	dat <- AorticStenosisTrials[which(AorticStenosisTrials$trialID==j),]
	index <- 0
	FU <- min(max(dat[which(dat$Arm==1),]$Time), max(dat[which(dat$Arm==0),]$Time)) 
	print(FU) 
	for (tau in t){
		index <- index+1
		if(FU>=tau){obj<-rmst2(dat$Time, dat$Event, dat$Arm, tau=tau)}else{obj <- NA}
		if(FU>=tau){RMSTcurveRes[index,j+1] <-  obj$unadjusted.result[1]}else{RMSTcurveRes[index,j+1] <- NA}
		RMSTcurveRes[index,1] <- tau	
}
}

#colnames(RMSTcurveRes) <- c("t",names(ASdatalist))
colnames(RMSTcurveRes) <- c("t","Notion", "Partner", "Surtavi", "Partner2", "Corevalve")

# fit Royston Parmar model in each trial
t <- seq(0.25,36,by=0.25)
RPres <- matrix(NA, length(t), J+1)
for (j in 1:J){
  dat <- AorticStenosisTrials[which(AorticStenosisTrials$trialID==j),]
  index <- 0
  MC <- stpm2(Surv(Time, Event==1)~Arm, data=dat, smooth.formula=~ns(log(Time),df=3)+log(Time):Arm)
  
  for (tau in t){
    index <- index+1
  RPres[index, j+1]	<- predict(MC, newdata=data.frame(Arm=1, Time=tau), type="rmst",se.fit=TRUE)$Estimate	- predict(MC, newdata=data.frame(Arm=0, Time=tau), type="rmst",se.fit=TRUE)$Estimate
  RPres[index,1] <- tau
  }
}

# calculate meta-analytic results by all methods:

mvma_res <- metaRMSTD(AorticStenosisTrials, time_horizons=c(12,24,36), MA_method="mvma")
boot_res <- metaRMSTD(AorticStenosisTrials, time_horizons=c(12,24,36), MA_method="mvma_boot", nboot=500)
uni_res <- metaRMSTD(AorticStenosisTrials, time_horizons=c(12,24,36), MA_method="uni")
uni_flex_res <- metaRMSTD(AorticStenosisTrials, time_horizons=c(12,24,36), MA_method="uni_flex")
results <- rbind(mvma_res$REresult, boot_res$REresult, uni_flex_res$result[,1:9], uni_res$result[,1:9])

### plot of RMSTD + RP predictions + MA results ###
plot(RMSTcurveRes[,1], RMSTcurveRes[,3], type="l", col=palette2[3], lwd=3, ylim=c(-0.75,2.75), xlim=c(0,36),
		main="", xlab="Time Horizon (months)", ylab="Difference in RMST (months)", xaxt="n", yaxt="n")
axis(2, at=c(-0.5, 0,0.5,1, 1.5, 2, 2.5))
axis(1, at=c(0,12,24,36))

lines(RMSTcurveRes[,1], RMSTcurveRes[,2], col=palette2[1], lwd=3)
lines(RMSTcurveRes[,1], RMSTcurveRes[,4], col=palette2[5], lwd=3)
lines(RMSTcurveRes[,1], RMSTcurveRes[,5], col=palette2[4], lwd=3)
lines(RMSTcurveRes[,1], RMSTcurveRes[,6], col=palette2[2], lwd=3)

# add RP lines
lines(RPres[,1], RPres[,2], col=palette2[1], lty=2, lwd=3)
lines(RPres[,1], RPres[,3], col=palette2[3], lty=2, lwd=3)
lines(RPres[,1], RPres[,4], col=palette2[5], lty=2, lwd=3)
lines(RPres[,1], RPres[,5], col=palette2[4], lty=2, lwd=3)
lines(RPres[,1], RPres[,6], col=palette2[2], lty=2, lwd=3)

# add meta-analysis results
# MVMA
points(results[1,2]-0.9, results[1,4], pch=19, col="black")
arrows(results[1,2]-0.9, results[1,6], results[1,2]-0.9, results[1,7], angle=90, code=3, length=0.05, col="black")
points(results[2,2]-0.9, results[2,4], pch=19, col="black")
arrows(results[2,2]-0.9, results[2,6], results[2,2]-0.9, results[2,7], angle=90, code=3, length=0.05, col="black")
points(results[3,2]-0.9, results[3,4], pch=19, col="black")
arrows(results[3,2]-0.9, results[3,6], results[3,2]-0.9, results[3,7], angle=90, code=3, length=0.05, col="black")

# MVMA boot
points(results[4,2]-0.3, results[4,4], col="black", pch=17)
arrows(results[4,2]-0.3, results[4,6], results[4,2]-0.3, results[4,7], angle=90, code=3, length=0.05, col="black")
points(results[5,2]-0.3, results[5,4], col="black", pch=17)
arrows(results[5,2]-0.3, results[5,6], results[5,2]-0.3, results[5,7], angle=90, code=3, length=0.05, col="black")
points(results[6,2]-0.3, results[6,4], col="black", pch=17)
arrows(results[6,2]-0.3, results[6,6], results[6,2]-0.3, results[6,7], angle=90, code=3, length=0.05, col="black")

# royston parmar
points(results[7,2]+0.3, results[7,4], col="black", pch=15)
arrows(results[7,2]+0.3, results[7,6], results[7,2]+0.3, results[7,7], angle=90, code=3, length=0.05, col="black")
points(results[8,2]+0.3, results[8,4], col="black", pch=15)
arrows(results[8,2]+0.3, results[8,6], results[8,2]+0.3, results[8,7], angle=90, code=3, length=0.05, col="black")
points(results[9,2]+0.3, results[9,4], col="black", pch=15)
arrows(results[9,2]+0.3, results[9,6], results[9,2]+0.3, results[9,7], angle=90, code=3, length=0.05, col="black")

# univariate
points(results[10,2]+0.9, results[10,4], col="black", pch=18)
arrows(results[10,2]+0.9, results[10,6], results[10,2]+0.9, results[10,7], angle=90, code=3, length=0.05, col="black")
points(results[11,2]+0.9, results[11,4], col="black", pch=18)
arrows(results[11,2]+0.9, results[11,6], results[11,2]+0.9, results[11,7], angle=90, code=3, length=0.05, col="black")
points(results[12,2]+0.9, results[12,4], col="black", pch=18)
arrows(results[12,2]+0.9, results[12,6], results[12,2]+0.9, results[12,7], angle=90, code=3, length=0.05, col="black")

# trial legend
legend(x=1, y=2.75, legend=c("Partner", "Corevalve", "Notion", "Partner2", "Surtavi"), 
		col=c(palette2[c(3,2,1,4,5)]), 
		horiz=F, lty=c(1,1,1,1,1), pch=c(NA,NA,NA,NA,NA), 
		lwd=c(3,3,3,3,3), seg.len=1, bty="n")

# MA method legend
legend(x=10, y=2.75, legend=c("Multivariate, analytic covariance", "Multivariate, bootstrap covariance", "Univariate, flexible parametric model estimates", "Univariate, available data" ), 
		col=c("black", "black", "black", "black"), 
		horiz=F, lty=c(1,1,1,1), pch=c(19,17,15,18), 
		lwd=c(1,1,1,1), seg.len=1, bty="n")

text(x=2, y=2.75, "Trial")
text(x=15.5, y=2.75, "Meta-analytic Method")

abline(h=0)


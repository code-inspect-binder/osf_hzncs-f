sampleSizes <- c(40,50,60,100)
sampleSizes <- seq(10,100,by=10)
effectSizes <- rep(.3,length(sampleSizes))
numStudies <- 1000

mean1 <- .50
mySDs <- c(.1,.1)

sdPooled <- sqrt((mySDs[1]^2 + mySDs[2]^2) / 2)
saveSims <- as.data.frame(matrix(0,ncol=7, nrow=length(sampleSizes) * (numStudies)))
colnames(saveSims) <- c("groupNum","modeledD","sampleSize","actualD","t","p","BF")
sa <- 0

for (i in 1:length(sampleSizes)) {
  N <- sampleSizes[i]
  
    mean2 <- mean1 - (effectSizes[i] * sdPooled)
    
    for (k in 1:numStudies) {
      
      #simulate data
      group1 <- rnorm(N,mean1,mySDs[1])
      group2 <- rnorm(N,mean2,mySDs[2])
      
      #run t-test
      if  (mySDs[1]==mySDs[2]) {
        a <- t.test(group2, group1,var.equal=T)
      } else {
        a <- t.test(group2, group1)
      }
      
      #save outcomes
      sa <- sa+1
      saveSims$modeledD[sa] <- effectSizes[i]
      saveSims$sampleSize[sa] <- N
      saveSims$actualD[sa] <- (mean(group1) - mean(group2))/sqrt((sd(group1)^2 + sd(group2)^2)/2)
      
      saveSims$p[sa] <- a$p.value
      saveSims$t[sa] <- a$statistic
#      saveSims$BF[sa] <- exp(ttest.tstat(t=a$statistic, n1=N, n2=N, rscale = 0.707)[['bf']])
      
    } #end for k numstudies  
} #end for i sampleSizes


print(paste("mean of all:",round(mean(saveSims$actualD),2)))
print(paste("mean of sig only:",round(mean(saveSims$actualD[which(saveSims$p < .05)]),2)))

mm <- aggregate(actualD ~ sampleSize + modeledD, saveSims, mean)
mm2 <- aggregate(actualD ~ sampleSize + modeledD, saveSims[which(saveSims$p <= .05),], mean)

mm3 <- aggregate(actualD ~ sampleSize + modeledD, saveSims, length)
mm4 <- aggregate(actualD ~ sampleSize + modeledD, saveSims[which(saveSims$p <= .05),], length)
mm4$prop <- mm4$actualD / mm3$actualD
mm4$cex <- mm4$prop * 10 + 1

xLab <- mm$sampleSize
for (i in 1:length(xLab)) {
  xLab[i] <- paste(mm$modeledD[i],"(",mm$sampleSize[i],")",sep="")
}

#plot(seq(1,length(xLab)), mm$actualD,xaxt="n",bty="l",ylim=c(min(mm$modeledD)/2,max(mm2$actualD)*1.2),xlab="Effect size (Sample size)", ylab="Effect size d")
#points(seq(1,length(xLab)),mm2$actualD,pch=19)
#legend("topright",pch=c(19,1),legend = c("Sig only","All"))
#axis(side=1,at=seq(1,length(xLab)), labels = xLab)
#abline(h=unique(effectSizes),col=rainbow(length(unique(effectSizes))))

plot(mm$sampleSize, mm$actualD,bty="l",ylim=c(min(mm$modeledD)/2,max(mm2$actualD)*1.2),xlab="Sample size", ylab="Cohen's d",cex=2)
points(mm$sampleSize,mm2$actualD,pch=19,cex=2)
legend("topright",pch=c(19,1),legend = c("Published","All"))


plot(mm$sampleSize, abs(mm$actualD),bty="l",ylim=c(min(mm$modeledD)/2,max(mm2$actualD)*1.2),xlab="Sample size", ylab="Absolute value of Cohen's d",cex=2)
points(mm$sampleSize,abs(mm2$actualD),pch=19,cex=mm4$cex)
legend("topright",pch=c(19,1),legend = c("Published","All"))


require(pwr)
a <- pwr.t.test(n=NULL, d = .5, sig.level = .05,power=.8,type="two.sample",alternative="two.sided")
b <- pwr.t.test(n=a$n, d = .25, sig.level = .05,power=NULL,type="two.sample",alternative="two.sided")
print(b$power)

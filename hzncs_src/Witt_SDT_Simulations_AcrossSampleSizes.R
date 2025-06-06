##
## Jessica K. Witt - Colorado State University - Psychology
##
## November 15, 2017 (last revised 3/28/18)
##
## Uses signal detection theory techniques to compare and contrast
## across different criteria for statistical significance
##
## Criterion are: (1) p < .05; (2) Bayes Factor > 3; etc...
##
## Hit = data were modeled as a real effect (Cohen's d > 0) and 
##       criteria indicates the effect is significant (e.g. p < .05)
## False alarm = data were modeled as a null effect (Cohen's d = 0) and 
##       criteria indicates the effect is significant (e.g. p < .05)
##
## Strategy: Simulate means from two groups and compute relevant statistics
## Mean for group 1 is .5 (as in 50% on a memory test)
## Mean for group 2 is a function of the specified effect size
##
## Two outcome measures:
##  (1) Area under the curve (AUC) - measure of discriminability
##  (2) Distance to perfection - Euclidean distance between a given point on the ROC curve and perfect performance (100% hits and 0% false alarms)




########## Packages #################
require(BayesFactor)
require(pROC)
require(ggplot2)
require(reshape2)



########## Functions #################

## Assumes 2 groups, so two in mySDs
runSims <- function(sampleSizes, effectSizes, mean1, mySDs, numStudies) {
  sdPooled <- sqrt((mySDs[1]^2 + mySDs[2]^2) / 2)
  allPs <- array(0,c(length(sampleSizes),length(effectSizes),numStudies))
  allBFs <- array(0,c(length(sampleSizes),length(effectSizes),numStudies))
  allESs <- array(0,c(length(sampleSizes),length(effectSizes),numStudies))
  allTs <- array(0,c(length(sampleSizes),length(effectSizes),numStudies))
  
  for (i in 1:length(sampleSizes)) {
    N <- sampleSizes[i]

    for (j in 1:length(effectSizes)) {
      mean2 <- mean1 - (effectSizes[j] * sdPooled)

      for (k in 1:numStudies) {
        
        group1 <- rnorm(N,mean1,mySDs[1])
        group2 <- rnorm(N,mean2,mySDs[2])

        if  (mySDs[1]==mySDs[2]) {
          a <- t.test(group2, group1,var.equal=T)
        } else {
          a <- t.test(group2, group1)
        }

        allPs[i,j,k] <- a$p.value
        allTs[i,j,k] <- a$statistic
        allBFs[i,j,k] <- exp(ttest.tstat(t=a$statistic, n1=N, n2=N, rscale = 0.707)[['bf']])
        allESs[i,j,k] <- (mean(group1) - mean(group2))/sqrt((sd(group1)^2 + sd(group2)^2)/2)
        
        
      } #end for k numstudies  
    } #end for j effectSizes
  } #end for i sampleSizes
  
  return(list(allPs,allBFs,allESs))
} #end function

### Example call
#myOut <- runSims(sampleSizes,effectSizes, mean1, mySDs, numStudies)
#allPs <- myOut[[1]]
#allBFs <- myOUt[[2]]
#allESs <- myOut[[3]]



#AUC Function
#tips taken from:
#https://www.r-bloggers.com/roc-curves-in-two-lines-of-r-code/

## example call:
## myAUCs <- plotAUC(allPs,allBFs,dt,numStudies,TRUE)

plotAUC <- function(allPs,allBFs,dt,numStudies,plotIt) {

  auc <- c(0,0)  #auc for p-value and Bayes factor
  
  sg <- c(log(allPs[1,2,]),log(allPs[1,1,]))
  sf <- c(rep(1,numStudies),rep(0,numStudies))

  fit_glm <- glm(sf ~ sg, family=binomial(link="logit"))
  glm_response_scores <- predict(fit_glm, data.frame(sg,sf), type="response")

  
  sg2 <- c(log(allBFs[1,2,]),log(allBFs[1,1,]))
  sf2 <- c(rep(1,numStudies),rep(0,numStudies))

  fit_glm2 <- glm(sf2 ~ sg2, family=binomial(link="logit"))
  glm_response_scores2 <- predict(fit_glm2, data.frame(sg2,sf2), type="response")


  auc[1] <- auc(sf, glm_response_scores)
  auc[2] <- auc(sf2, glm_response_scores2)


  if(plotIt) {  
    plot(roc(sf, glm_response_scores, direction="<"), col="green", lwd=4)
  
    legend("bottomright",col=seq(1:8),legend = crits,pch=19)
    abline(v=1)
    abline(v=0)
  
    lines(roc(sf2, glm_response_scores2, direction="<"), col="blue", lwd=1)

    points(1-(dt$fa/numStudies),dt$hits/numStudies,col=dt$criterion,pch=19)
  } #end if plotIt
  
  return(auc)
}  # end function plotAUC




## sample call:
# dt <- evalSig(effectSizes,critCat,critValue,allPs,allBFs,allESs)

evalSig <- function(effectSizes,critCat,critValue,allPs,allBFs,allESs) {
  
  numEffectSizes <- length(effectSizes) - 1
  numCriteria <- length(critCat)
  allC <- length(sampleSizes) * numEffectSizes * numCriteria 
  dt <- as.data.frame(matrix(0,ncol=12,nrow=allC))
  colnames(dt) <- c("sampleSize","effectSizes","criterion","hits","fa","actEffectSize","criterionF","miss","corrRej","numAmbigEff","numAmbigNull", "ROCdist")
  da <- 0  #index for dt

  for(i in 1:length(sampleSizes)) {
    for(j in 1:numEffectSizes) {
      for (k in 1:numCriteria) {
      
        da <- da+1
        dt$sampleSize[da] <- sampleSizes[i]
        dt$effectSizes[da] <- effectSizes[j+1]
        dt$criterion[da] <- k
      
        jj <- j+1  #assumes first effect size is 0, working with next one
        
        if (critCat[k] == 1) {  #p-value
          dt$hits[da] <- length(which(allPs[i,jj,] <= critValue[k]))
          dt$fa[da] <-length(which(allPs[i,1,] <= critValue[k]))
          dt$miss[da] <- numStudies - dt$hits[da]
          dt$corrRej[da] <- numStudies - dt$fa[da]
          dt$criterionF[da] <- paste("p",critValue[k])
        } else {  #Bayes factor
          dt$hits[da] <- length(which(allBFs[i,jj,] >= critValue[k]))
          dt$fa[da] <- length(which(allBFs[i,1,] >= critValue[k]))
          dt$miss[da] <- length(which(allBFs[i,jj,] <= (1/critValue[k])))
          dt$corrRej[da] <- length(which(allBFs[i,1,] <= (1/critValue[k])))
          dt$numAmbigEff[da] <- length(which(allBFs[i,jj,] < critValue[k] & allBFs[i,jj,] > (1/critValue[k]) ))
          dt$numAmbigNull[da] <- length(which(allBFs[i,1,] < critValue[k] & allBFs[i,1,] > (1/critValue[k]) ))
          dt$criterionF[da] <- paste("BF",critValue[k])
        } 
      
      h1 <- 1 - (dt$hits[da]/numStudies)  #vertical distance
      h2 <- dt$fa[da]/numStudies  #horizDistance
      dt$ROCdist[da] <- sqrt(h1^2 + h2^2)
      
      dt$actEffectSize[da] <- mean(allESs[i,jj,])
      
      }
   }
  }
  
  return(dt)
}



############################
####### Start Code #########
############################


#Num of studies to be run and number of times that number of studies should be run
#numStudies <- 10 #for Fig 4, panel a; Fig 8
#numTimes <- 100 #for Fig 4, panel a; Fig 8


numStudies <- 1000 #for Fig 4, panels b,c and d; Fig 6
numTimes <- 1  #for Fig 4, panels b,c and d; Fig 6


#Sample sizes and Effect sizes (note: sample size = sample size per group; assumes equal samples per group)
sSizesAll <- c(50,176,290,64,105,26,45)
effSzAll <- c(.3,.3,.3,.5,.5,.8,.8)
sSizesAll <- seq(32,200,length.out=7)  #for Fig 8
sSizesAll[length(sSizesAll)+1] <- 2000 #for Fig 8
sSizesAll <- seq(20,2000,length.out=30) #for Fig 4
sSizesAll <- c(seq(5,100,length.out=20),seq(105,2000,length.out = 10)) #for Fig 6
sSizesAll <- round(sSizesAll)
effSzAll <- rep(.5,length(sSizesAll))

#data frame to save SDT analyses
allDt <- as.data.frame(matrix(0,ncol=13,nrow=1))
colnames(allDt) <- c("sampleSize","effectSizes","criterion","hits","fa","actEffectSize","criterionF","miss","corrRej","numAmbigEff","numAmbigNull", "ROCdist","runNum")

#criteria for significance to be used
critCat <- c(1,1,1,1,2,2,2,2) # 1 = p-value; 2 = BF
critValue <- c(.1,.05,.005,.001,1,2,3,10)  #critical value for significance,, aligns with critCat
crits <- c("p<.10","p<.05","p<.005","p<.001","BF>1","BF>2","BF>3","BF>10")  #label not needed for function

#Group mean and SDs
mean1 <- .5
mySDs <- c(.1, .1)

#AUC estimates
aucs <- as.data.frame(matrix(0,ncol=5,nrow=length(sSizesAll)*numTimes))
colnames (aucs) <- c("sampleSize","effectSize","runTime","aucP","aucBF")
ai <- 0

###### Loops for Data simulation and SDT calculations ####

##if want to do 1 sample size at a time  (Figs 4a, 8)
#for (ssa in 1:length(sSizesAll)) {
  #sampleSizes <- sSizesAll[ssa]  

## if want to do all sample sizes at once (Figs 4b,c,d; Fig 6)
for (ssa in 1:1) {
  sampleSizes <- sSizesAll      
    
  effectSizes <- c(0,effSzAll[ssa])  #sets first effect size to 0; necessary for SDT analyses
  
  for(i in 1:numTimes) {
    
    ## Run simulations
    myOut <- runSims(sampleSizes,effectSizes, mean1, mySDs, numStudies)
    allPs <- myOut[[1]]
    allBFs <- myOut[[2]]
    allESs <- myOut[[3]]
    

    ## Run SDT analyses
    ##critical: Assumes the first effect size is 0
    dt <- evalSig(effectSizes,critCat,critValue,allPs,allBFs,allESs)
    dt$runNum <- i
    allDt <- rbind(allDt,dt)
    
    # Plot ROC curves and caculate AUC (if desired)
    plotIt <- FALSE
    currAUC <- plotAUC(allPs,allBFs,dt,numStudies,plotIt)  #currently does not save the AUCs
    ai <- ai+1
    aucs[ai,1] <- sampleSizes[1]
    aucs[ai,2] <- effectSizes[2]
    aucs[ai,3] <- i #runNumber
    aucs[ai,4] <- currAUC[1]
    aucs[ai,5] <- currAUC[2]

  } #end for i runTimes
} #end for ssa

allDt <- allDt[-1,]  #get rid of initial dummy row that was just used to init allDt

######################
###### Analyses ######
######################


#### Compare AUCs  (Fig 4a) ####
ap <- aggregate(aucP ~ sampleSize + effectSize,aucs,mean)
ap2 <- aggregate(aucBF ~ sampleSize + effectSize,aucs,mean)
plot(seq(1:length(sSizesAll)),ap$aucP,col=rainbow(length(sSizesAll),start=0,end=.8),bty="l",ylim=c(.7,1),pch=1,cex=2,xaxt="n",xlab="N (per group)",ylab="Area Under the Curve")
points(seq(1:length(sSizesAll)),ap2$aucBF,col=rainbow(length(sSizesAll),start=0,end=.8),pch=19)
legend("bottomleft",pch=c(1,19),pt.cex=c(2,1),legend = c("p-values","Bayes factors"))
la <- c(1,length(sSizesAll)/3,length(sSizesAll)/3*2,length(sSizesAll))
axis(side=1,at=la,labels = sSizesAll[la])


#### Plot distance to perfection  (fig 8) ####
aa <- aggregate(ROCdist ~ criterion + sampleSize, data=allDt,mean)
ab <- aggregate(ROCdist ~ criterion + sampleSize, data=allDt,sd)
ab$ci <- qnorm(.975) * ab$ROCdist / sqrt(numStudies)
ac <- length(aa$criterion)
ad <- length(unique(aa$criterion))
aa$sampleSize[which(aa$sampleSize == 2000)] <- 250 
plot(aa$sampleSize+(2*aa$criterion),aa$ROCdist,bty="l",xaxt="n",xlab="Sample Size (n per group)",ylab="Distance to Perfection",pch=19,col=rainbow(length(crits)),cex=2,ylim=c(0,max(aa$ROCdist)+max(ab$ci)))
axis(side=1,at=sort(unique(aa$sampleSize)) + length(crits),labels = sSizesAll)
legend("topright",pch=19,col=rainbow(length(crits)),legend = crits)
for(i in 1:ac) {
#  segments(aa$sampleSize[i]+(2*aa$criterion[i]),aa$ROCdist[i] - ab$ci[i],aa$sampleSize[i]+(2*aa$criterion[i]),aa$ROCdist[i]+ab$ci[i])
  segments(aa$sampleSize[i]+(2*aa$criterion[i]),aa$ROCdist[i] - ab$ci[i],aa$sampleSize[i]+(2*aa$criterion[i]),aa$ROCdist[i]+ab$ci[i],col=rainbow(length(crits))[aa$criterion[i]])
}



#### Plot Fig 4b,c,d ####
#plot BF vs p-value across sample Size
if (1>0) {
  saveLM <- as.data.frame(matrix(0,ncol=8,nrow=length(sampleSizes)))
  colnames(saveLM) <- c("N","intercept","slope","rSq","p10","p05","p005","p44")
#  plot(log(allPs), log(allBFs),bty="l",col=rainbow(length(sampleSizes)),xlab="p-value",ylab="Bayes Factor",yaxt="n",xaxt="n",ylim=c(-4,15),xlim=c(-20,0),cex.lab=1.5)
#  plot(log(allPs), log(allBFs),bty="l",col=rainbow(length(sampleSizes)),xlab="p-value",ylab="Bayes Factor",yaxt="n",xaxt="n",ylim=c(-4,5),xlim=c(-10,0),cex.lab=1.5)
  plot(log(allPs), log(allBFs),bty="l",pch=19,col=rainbow(length(sampleSizes),start=0,end=.8),xlab="p-value",ylab="Bayes Factor",yaxt="n",xaxt="n",ylim=c(-4,5),xlim=c(-10,0),cex.lab=1.5)
  abline(h=log(1))
  abline(v=log(.05))
  ya <- c(.3,3,30)
  axis(side=2,at=log(ya),labels = ya,cex.axis=1.3)
  xa <- c(.0000001,.001,.01,.05,.5)
  axis(side=1,at=log(xa),labels = xa,cex.axis=1.3)
  la <- c(sampleSizes[1],sampleSizes[length(sampleSizes)/3],sampleSizes[length(sampleSizes)/3*2],sampleSizes[length(sampleSizes)])
  la2 <- c(1,length(sampleSizes)/3,length(sampleSizes)/3*2,length(sampleSizes))
  legend("topright",pch=rep(19,4),col=rainbow(length(sampleSizes),start=0,end=.8)[la2],legend = la,bty="n",cex=1.5)
  
  plot(log(allPs), log(allBFs),pch=19,bty="l",col=rainbow(length(sampleSizes),start=0,end=.8),xlab="log(p-value)",ylab="log(Bayes Factor)",yaxt="n",xaxt="n",cex.lab=1.5)
  rect(-10,-4,0,5,lty=2)  #matches xlim and ylim in plot above
  legend("topright",pch=rep(1,4),col=rainbow(4),legend = la,bty="n",cex=1.2)

  #add shading to inconclusive BF
  plot(log(allPs), log(allBFs),bty="l",pch=19,col=rainbow(length(sampleSizes),start=0,end=.8),xlab="p-value",ylab="Bayes Factor",yaxt="n",xaxt="n",ylim=c(-4,5),xlim=c(-10,0),cex.lab=1.5)
  abline(h=log(1))
  abline(v=log(.05))
  ya <- c(.3,3,30)
  axis(side=2,at=log(ya),labels = ya,cex.axis=1.3)
  xa <- c(.0000001,.001,.01,.05,.5)
  axis(side=1,at=log(xa),labels = xa,cex.axis=1.3)
  la <- c(sampleSizes[1],sampleSizes[length(sampleSizes)/3],sampleSizes[length(sampleSizes)/3*2],sampleSizes[length(sampleSizes)])
  la2 <- c(1,length(sampleSizes)/3,length(sampleSizes)/3*2,length(sampleSizes))
  legend("topright",pch=rep(19,4),col=rainbow(length(sampleSizes),start=0,end=.8)[la2],legend = la,bty="n",cex=1.5)
  rect(log(.00000001),log(.3), log(1.1),log(3), col= rgb(.5,.5,.5,alpha=0.3),border=NA)
  
  
  plot(log(allPs), log(allBFs),col="white",xlab="log(p-value)",ylab="log(Bayes Factor")
  for (compIndex in 1:length(sampleSizes)) {
    points(log(allPs[compIndex,,]), log(allBFs[compIndex,,]),col=rainbow(compIndex))
    
    ty <- lm(log(allBFs[compIndex,,]) ~ log(allPs[compIndex,,]))
    saveLM[compIndex,1] <- sampleSizes[compIndex]
    saveLM[compIndex,2] <- ty$coefficients[1]
    saveLM[compIndex,3] <- ty$coefficients[2]
    saveLM[compIndex,4] <- summary.lm(ty)$r.squared
    saveLM[compIndex,5] <- exp(ty$coefficients[1] + ty$coefficients[2] * log(.1))
    saveLM[compIndex,6] <- exp(ty$coefficients[1] + ty$coefficients[2] * log(.05))
    saveLM[compIndex,7] <- exp(ty$coefficients[1] + ty$coefficients[2] * log(.005))
    saveLM[compIndex,8] <- exp(ty$coefficients[1] + ty$coefficients[2] * log(.44))
    
  } #end for compIndex
  
  #### Plot Fig 4c and d ####
  
  if(1>0) {  #plot results of LM between p and BFs
    plot(saveLM[,1],saveLM[,2],bty="l",xlab="N (per group)", ylab="Intercept",pch=19)
    plot(saveLM[,1],exp(saveLM[,2]),bty="l",xlab="N (per group)", ylab="Intercept (as a p-value)",pch=19)
    plot(saveLM[,1],saveLM[,3],bty="l",xlab="N (per group)", ylab="Slope",pch=19,ylim=c(-1,-.7))
    plot(saveLM[,1],saveLM[,4],bty="l",xlab="N (per group)", ylab="R^2",pch=19)
    
    plot(saveLM[,1],exp(saveLM[,2]),bty="l",xlab="N (per group)", ylab="Intercept (as a p-value)",cex=2,pch=19,col=rainbow(length(sampleSizes),start=0,end=.8))
    plot(saveLM[,1],saveLM[,3],bty="l",xlab="N (per group)", ylab="Slope",cex=2,pch=19,ylim=c(-1,-.7),col=rainbow(length(sampleSizes),start=0,end=.8))

    somePs <- seq(.00001,1,length.out = 1000)
#    plot(log(somePs),saveLM$intercept[1] + saveLM$slope[1]*log(somePs),type="l",bty="l",xlab="p-value",ylab="Bayes factor")
#    for (ipi in 2:length(saveLM[,1])) {
#      lines(log(somePs),saveLM$intercept[ipi] + saveLM$slope[ipi]*log(somePs),col=ipi)
#    }
#    abline(v=log(.05),lwd=2)
#    abline(v=log(.005))
#    abline(h=log(1))
#    legend("topright",lwd=1,col=seq(1:length(saveLM$N)),legend = sort(unique(saveLM$N)))
    
  } #end if 
} #end if plot LM p vs BF (fig 4)



#### Min/Max value plots (Fig 6) ####
n2 <- 20

#for studies modeled as null effect
plot(sampleSizes[1:n2],seq(0,max(1/allBFs[1:n2,1,]),length.out = n2),col="white",bty="l",ylab="1 / Bayes factor",xlab= "N")
for (i in 1:n2) {
  points(rep(sampleSizes[i],1000),1/allBFs[i,1,],col=rainbow(n2,start=0,end=.8)[i])
}

plot(sampleSizes[1:n2],seq(0,1,length.out = n2),col="white",bty="l",ylab="p-value",xlab= "N")
for (i in 1:n2) {
  points(rep(sampleSizes[i],1000),allPs[i,1,],col=rainbow(n2,start=0,end=.8)[i])
}

#for studies modeled as null effect w jitter
plot(sampleSizes[1:n2],seq(0,max(1/allBFs[1:n2,1,]),length.out = n2),col="white",bty="l",ylab="1 / Bayes factor",xlab= "N")
for (i in 1:n2) {
  points(jitter(rep(sampleSizes[i],1000)),1/allBFs[i,1,],col=rainbow(n2,start=0,end=.8)[i],cex=.2)
}

plot(sampleSizes[1:n2],seq(0,1,length.out = n2),col="white",bty="l",ylab="p-value",xlab= "N")
for (i in 1:n2) {
  points(jitter(rep(sampleSizes[i],1000)),allPs[i,1,],col=rainbow(n2,start=0,end=.8)[i], cex=.2)
}


#modeled as real effects
plot(sampleSizes[1:n2],seq(0,max(allBFs[1:n2,2,]),length.out = n2),col="white",bty="l",ylab="Bayes factor",xlab= "N")
for (i in 1:n2) {
  points(rep(sampleSizes[i],1000),allBFs[i,2,],col=rainbow(n2,start=0,end=.8)[i])
}

plot(sampleSizes[1:n2],seq(0,1,length.out = n2),col="white",bty="l",ylab="p-value",xlab= "N")
for (i in 1:n2) {
  points(rep(sampleSizes[i],1000),allPs[i,2,],col=rainbow(n2,start=0,end=.8)[i])
}

#modeled as real effects w jitter
plot(sampleSizes[1:n2],seq(0,max(log(allBFs[1:n2,2,])),length.out = n2),col="white",bty="l",ylab="log(Bayes factor)",xlab= "N")
for (i in 1:n2) {
  points(jitter(rep(sampleSizes[i],1000)),log(allBFs[i,2,]),col=rainbow(n2,start=0,end=.8)[i],cex=.2)
}
abline(h=log(3))

plot(sampleSizes[1:n2],seq(0,1,length.out = n2),col="white",bty="l",ylab="p-value",xlab= "N")
for (i in 1:n2) {
  points(jitter(rep(sampleSizes[i],1000)),allPs[i,2,],col=rainbow(n2,start=0,end=.8)[i],cex=.2)
}


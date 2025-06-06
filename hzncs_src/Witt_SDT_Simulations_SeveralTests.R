##
## Jessica K. Witt - Colorado State University - Psychology
##
## November 15, 2017 (last revised 1/11/18)
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
## Strategy: Simulate means from a study and compute relevant statistics
## Studies include two-sample t-test, one-sample t-test, uneven sample sizes two-sample t-test, and correlation
##
## Two outcome measures:
##  (1) Area under the curve (AUC) - measure of discriminability
##  (2) Distance to perfection - Euclidean distance between a given point on the ROC curve and perfect performance (100% hits and 0% false alarms)



########## Packages #################
require(BayesFactor)
require(pROC)
require(ggplot2)
require(reshape2)
require(pwr)
require(MASS)



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

runSimsUneven <- function(sampleSizes, effectSizes, mean1, mySDs, numStudies) {
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
        group2 <- rnorm(N*1.2,mean2,mySDs[2])  #20% more participants in one group
        
        if  (mySDs[1]==mySDs[2]) {
          a <- t.test(group2, group1,var.equal=T)
        } else {
          a <- t.test(group2, group1)
        }
        
        allPs[i,j,k] <- a$p.value
        allTs[i,j,k] <- a$statistic
        allBFs[i,j,k] <- exp(ttest.tstat(t=a$statistic, n1=length(group1), n2=length(group2), rscale = 0.707)[['bf']])
        allESs[i,j,k] <- (mean(group1) - mean(group2))/sqrt((sd(group1)^2 + sd(group2)^2)/2)
        
        
      } #end for k numstudies  
    } #end for j effectSizes
  } #end for i sampleSizes
  
  return(list(allPs,allBFs,allESs))
} #end function



runSimsOne <- function(sampleSizes, effectSizes, mean1, mySDs, numStudies) {
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
        
        group2 <- rnorm(N,mean2,mySDs[2])
        
        a <- t.test(group2, mu = mean1)
        
        allPs[i,j,k] <- a$p.value
        allTs[i,j,k] <- a$statistic
        allBFs[i,j,k] <- exp(ttest.tstat(t=a$statistic, n1=N, n2=N, rscale = 0.707)[['bf']])
        allESs[i,j,k] <- a$statistic / sqrt(N)
        
        
      } #end for k numstudies  
    } #end for j effectSizes
  } #end for i sampleSizes
  
  return(list(allPs,allBFs,allESs))
} #end function


runSimsCorr <- function(sampleSizes, effectSizes, mean1, mySDs, numStudies) {
  sdPooled <- sqrt((mySDs[1]^2 + mySDs[2]^2) / 2)
  allPs <- array(0,c(length(sampleSizes),length(effectSizes),numStudies))
  allBFs <- array(0,c(length(sampleSizes),length(effectSizes),numStudies))
  allESs <- array(0,c(length(sampleSizes),length(effectSizes),numStudies))
  allTs <- array(0,c(length(sampleSizes),length(effectSizes),numStudies))
  
  for (i in 1:length(sampleSizes)) {
    N <- sampleSizes[i]
    
    for (j in 1:length(effectSizes)) {
      for (k in 1:numStudies) {
        
        gg <- mvrnorm(n=N,mu=c(mean1,mean1),Sigma = matrix(c(1,effectSizes[j],effectSizes[j],1),nrow=2),empirical=F)
        group1 <- gg[,1]
        group2 <- gg[,2]
        
        a <- cor.test(group2, group1,var.equal=T)
        
        allPs[i,j,k] <- a$p.value
        allTs[i,j,k] <- a$statistic
        allBFs[i,j,k] <- exp(ttest.tstat(t=a$statistic, n1=N, rscale = 0.707)[['bf']])
        allESs[i,j,k] <- a$estimate
        
        
      } #end for k numstudies  
    } #end for j effectSizes
  } #end for i sampleSizes
  
  return(list(allPs,allBFs,allESs))
} #end function



#AUC Function
#tips taken from:
#https://www.r-bloggers.com/roc-curves-in-two-lines-of-r-code/

## example call:
## myAUCs <- plotAUC(allPs,allBFs,dt,numStudies,TRUE)

plotAUC <- function(allPs,allBFs,dt,numStudies,plotIt) {

  auc <- c(0,0,0)  #auc for p-value and Bayes factor
  
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

  sg3 <- c(allESs[1,2,],allESs[1,1,])
  sf3 <- c(rep(1,numStudies),rep(0,numStudies))
  
  fit_glm3 <- glm(sf3 ~ sg3, family=binomial(link="logit"))
  glm_response_scores3 <- predict(fit_glm3, data.frame(sg3,sf3), type="response")
  
  auc[3] <- auc(sf3, glm_response_scores3)

  
  if(plotIt) {  
#    plot(roc(sf, glm_response_scores, direction="<"), col="green", lwd=4)
    plot.roc(sf, glm_response_scores, direction="<",legacy.axes = T, col="green", lwd=4,xlab="False Alarm Rate",ylab="Hit Rate",auc.polygon=T)
    
    legend("bottomright",col=seq(1:8),legend = crits,pch=19)
    abline(v=1)
    abline(v=0)
  
    lines(roc(sf2, glm_response_scores2, direction="<"), col="blue", lwd=1)
    lines(roc(sf3, glm_response_scores3, direction="<"), col="red", lwd=.5)
    
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
        } else if (critCat[k] == 2) {  #Bayes factor
          dt$hits[da] <- length(which(allBFs[i,jj,] >= critValue[k]))
          dt$fa[da] <- length(which(allBFs[i,1,] >= critValue[k]))
          dt$miss[da] <- length(which(allBFs[i,jj,] <= (1/critValue[k])))
          dt$corrRej[da] <- length(which(allBFs[i,1,] <= (1/critValue[k])))
          dt$numAmbigEff[da] <- length(which(allBFs[i,jj,] < critValue[k] & allBFs[i,jj,] > (1/critValue[k]) ))
          dt$numAmbigNull[da] <- length(which(allBFs[i,1,] < critValue[k] & allBFs[i,1,] > (1/critValue[k]) ))
          dt$criterionF[da] <- paste("BF",critValue[k])
        } else {  #effectSize
          dt$hits[da] <- length(which(allESs[i,jj,] >= critValue[k]))
          dt$fa[da] <-length(which(allESs[i,1,] >= critValue[k]))
          dt$miss[da] <- numStudies - dt$hits[da]
          dt$corrRej[da] <- numStudies - dt$fa[da]
          dt$criterionF[da] <- paste("ES",critValue[k])
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
numStudies <- 10
numTimes <- 100

#Sample sizes and Effect sizes (note: sample size = sample size per group; assumes equal samples per group)
sSizesAll <- c(50,176,290,64,105,26,45)
effSzAll <- c(.3,.3,.3,.5,.5,.8,.8)
#sSizesAll <- c(176)
#effSzAll <- c(.3)
# for d=.5, n=64 = 80%; n=105; 95% power  [n = num per group for 2 independent samples t test]
# for d=.3, n=176 = 80%; n=290; 95% power
# for d=.8, n=26 = 80%; n=42; 95% power
numES <- 3
effSzAll <- rep(seq(.1,1,length.out = numES),3)
effSzAll <- rep(c(.3,.5,.8),3)
ppower <- c(rep(.8,numES),rep(.9,numES),rep(.95,numES))
sSizesAll <- effSzAll
sSizesAllOne <- effSzAll
sSizesAllCorr <- effSzAll
for(i in 1:length(effSzAll)) {
  a <- pwr.t.test(n=NULL,d=effSzAll[i],sig.level = .05,power=ppower[i],type="two.sample",alternative="two.sided")
  sSizesAll[i] <- ceiling(a$n)

  a <- pwr.t.test(n=NULL,d=effSzAll[i],sig.level = .05,power=ppower[i],type="one.sample",alternative="two.sided")
  sSizesAllOne[i] <- ceiling(a$n)
  
  a <- pwr.r.test(n=NULL,r=effSzAll[i]/2,sig.level = .05,power=ppower[i],alternative="two.sided")
  sSizesAllCorr[i] <- ceiling(a$n)
}

testTypes <- c("two-sample t test","one-sample t test","uneven two-sample t","correlation")

#data frame to save SDT analyses
allDt <- as.data.frame(matrix(0,ncol=15,nrow=1))
colnames(allDt) <- c("sampleSize","effectSizes","criterion","hits","fa","actEffectSize","criterionF","miss","corrRej","numAmbigEff","numAmbigNull", "ROCdist","runNum","testType","power")

#criteria for significance to be used
critCat <- c(1,1,1,1,2,2,2,2) # 1 = p-value; 2 = BF
critValue <- c(.1,.05,.005,.001,1,2,3,10)  #critical value for significance,, aligns with critCat
crits <- c("p<.10","p<.05","p<.005","p<.001","BF>1","BF>2","BF>3","BF>10")  #label not needed for function

critCat <- c(1,1,1,2,2) # 1 = p-value; 2 = BF
critValue <- c(.1,.05,.005,1,3)  #critical value for significance,, aligns with critCat
crits <- c("p<.10","p<.05","p<.005","BF>1","BF>3")  #label not needed for function

#critCat <- c(1,1,1,2,2,3,3,3,3) # 1 = p-value; 2 = BF; 3 = ES
#critValue <- c(.1,.05,.005,1,3,.1,.3,.5,.8)  #critical value for significance,, aligns with critCat
#crits <- c("p<.10","p<.05","p<.005","BF>1","BF>3","ES>.1","ES>.3","ES>.5","ES>.8")  #label not needed for function

#Group mean and SDs
mean1 <- .5
mySDs <- c(.1, .1)

#AUC estimates
aucs <- as.data.frame(matrix(0,ncol=8,nrow=length(sSizesAll)*numTimes))
colnames (aucs) <- c("sampleSize","effectSize","runTime","aucP","aucBF","aucES","testType","power")
ai <- 0
plotIt <- FALSE

###### Loops for Data simulation and SDT calculations ####
for (ssa in 1:length(sSizesAll)) {
  effectSizes <- c(0,effSzAll[ssa])  #sets first effect size to 0; necessary for SDT analyses
  
  for(i in 1:numTimes) {
    
    ## two-sample t-tests
    ## Run simulations
    sampleSizes <- sSizesAll[ssa]
    myOut <- runSims(sampleSizes,effectSizes, mean1, mySDs, numStudies)
    allPs <- myOut[[1]]
    allBFs <- myOut[[2]]
    allESs <- myOut[[3]]
    ## Run SDT analyses
    ##critical: Assumes the first effect size is 0
    dt <- evalSig(effectSizes,critCat,critValue,allPs,allBFs,allESs)
    dt$runNum <- i
    dt$testType <- 1
    dt$power <- ppower[ssa]
    allDt <- rbind(allDt,dt)
    # Plot ROC curves and caculate AUC (if desired)
    currAUC <- plotAUC(allPs,allBFs,dt,numStudies,plotIt)  #currently does not save the AUCs
    ai <- ai+1
    aucs[ai,1] <- sampleSizes[1]
    aucs[ai,2] <- effectSizes[2]
    aucs[ai,3] <- i #runNumber
    aucs[ai,4] <- currAUC[1]
    aucs[ai,5] <- currAUC[2]
    aucs[ai,6] <- currAUC[3]
    aucs[ai,7] <- 1
    aucs[ai,8] <- ppower[ssa]
    
    
    ## one-sample t-tests
    ## Run simulations
    sampleSizes <- sSizesAllOne[ssa]
    myOut <- runSimsOne(sampleSizes,effectSizes, mean1, mySDs, numStudies)
    allPs <- myOut[[1]]
    allBFs <- myOut[[2]]
    allESs <- myOut[[3]]
    ## Run SDT analyses
    ##critical: Assumes the first effect size is 0
    dt <- evalSig(effectSizes,critCat,critValue,allPs,allBFs,allESs)
    dt$runNum <- i
    dt$testType <- 2
    dt$power <- ppower[ssa]
    allDt <- rbind(allDt,dt)
    # Plot ROC curves and caculate AUC (if desired)
    currAUC <- plotAUC(allPs,allBFs,dt,numStudies,plotIt)  #currently does not save the AUCs
    ai <- ai+1
    aucs[ai,1] <- sampleSizes[1]
    aucs[ai,2] <- effectSizes[2]
    aucs[ai,3] <- i #runNumber
    aucs[ai,4] <- currAUC[1]
    aucs[ai,5] <- currAUC[2]
    aucs[ai,6] <- currAUC[3]
    aucs[ai,7] <- 2
    aucs[ai,8] <- ppower[ssa]
    

    
    ## uneven sample sizes
    ## Run simulations
    sampleSizes <- sSizesAll[ssa]  #uses power test that assumes equal sample sizes even though I didn't simulate equal
    myOut <- runSimsUneven(sampleSizes,effectSizes, mean1, mySDs, numStudies)
    allPs <- myOut[[1]]
    allBFs <- myOut[[2]]
    allESs <- myOut[[3]]
    ## Run SDT analyses
    ##critical: Assumes the first effect size is 0
    dt <- evalSig(effectSizes,critCat,critValue,allPs,allBFs,allESs)
    dt$runNum <- i
    dt$testType <- 3
    dt$power <- ppower[ssa]
    allDt <- rbind(allDt,dt)
    # Plot ROC curves and caculate AUC (if desired)
    currAUC <- plotAUC(allPs,allBFs,dt,numStudies,plotIt)  #currently does not save the AUCs
    ai <- ai+1
    aucs[ai,1] <- sampleSizes[1]
    aucs[ai,2] <- effectSizes[2]
    aucs[ai,3] <- i #runNumber
    aucs[ai,4] <- currAUC[1]
    aucs[ai,5] <- currAUC[2]
    aucs[ai,6] <- currAUC[3]
    aucs[ai,7] <- 3
    aucs[ai,8] <- ppower[ssa]
    

    ## correlations
    ## Run simulations
    sampleSizes <- sSizesAllCorr[ssa]
    myOut <- runSimsCorr(sampleSizes,effectSizes/2, mean1, mySDs, numStudies)
    allPs <- myOut[[1]]
    allBFs <- myOut[[2]]
    allESs <- myOut[[3]]
    ## Run SDT analyses
    ##critical: Assumes the first effect size is 0
    dt <- evalSig(effectSizes,critCat,critValue,allPs,allBFs,allESs)
    dt$runNum <- i
    dt$testType <- 4
    dt$power <- ppower[ssa]
    allDt <- rbind(allDt,dt)
    # Plot ROC curves and caculate AUC (if desired)
    currAUC <- plotAUC(allPs,allBFs,dt,numStudies,plotIt)  #currently does not save the AUCs
    ai <- ai+1
    aucs[ai,1] <- sampleSizes[1]
    aucs[ai,2] <- effectSizes[2]
    aucs[ai,3] <- i #runNumber
    aucs[ai,4] <- currAUC[1]
    aucs[ai,5] <- currAUC[2]
    aucs[ai,6] <- currAUC[3]
    aucs[ai,7] <- 4
    aucs[ai,8] <- ppower[ssa]
    

  } #end for i runTimes
} #end for ssa

allDt <- allDt[-1,]  #get rid of initial dummy row that was just used to init allDt

######################
###### Analyses ######
######################

#### Compare AUCs ####
plot(aucs$aucP,aucs$aucES,col=aucs$effectSize*10,bty="l",ylim=c(.5,1),xlab="Area Under Curve for P-Value",ylab="AUC for Bayes factor (or Effect Size)")
points(aucs$aucP,aucs$aucBF,pch=19)
legend("bottomright",pch=c(1,19),legend = c("Effect size","Bayes factor"))
print(paste("AUC for effect size is better",length(which(aucs$aucES > aucs$aucP)),"times out of",length(aucs$aucES)))
ap <- aggregate(aucP ~ sampleSize + effectSize,aucs,mean)
ap2 <- aggregate(aucBF ~ sampleSize + effectSize,aucs,mean)
auc <- merge(ap,ap2,by=c("sampleSize","effectSize"))
ap2 <- aggregate(aucES ~ sampleSize + effectSize,aucs,mean)
auc <- merge(auc,ap2,by=c("sampleSize","effectSize"))
auc$x <- rank(auc$effectSize*10000000 + auc$sampleSize)
auc$xLab <- NA
for (i in 1:length(auc$x)) {
  auc$xLab[i] <- paste(auc$effectSize[i],"(",auc$sampleSize[i],")",sep="")
}
plot(auc$x,auc$aucP,col=auc$effectSize*10,bty="l",ylim=c(.7,1),pch=1,cex=2,xaxt="n",xlab="Effect Size (N)",ylab="Area Under the Curve")
axis(side=1,at=auc$x,auc$xLab)
points(auc$x,auc$aucBF,col=auc$effectSize*10,pch=19)
points(auc$x,auc$aucES,col=auc$effectSize*10,pch=18,cex=2)
legend("bottomright",pch=c(18,1,19),legend = c("for effect sizes","for p-values","for Bayes factors"))

#### Figure 9 ####
pHH <- ifelse(aucs$power==.80,0,ifelse(aucs$power==.90,3,4))
pCX <- ifelse(aucs$power==.80,1,ifelse(aucs$power==.90,1.3,1.5))
ccol <- c("red","orange","green","blue")
plot(aucs$aucP,aucs$aucBF,col=ccol[aucs$testType],cex = 1+aucs$effectSize, pch = pHH,bty="l",ylim=c(.5,1),xlim=c(.5,1),xlab="Area Under Curve for P-Value",ylab="AUC for Bayes factor")
legend("topleft",pch=19,col=ccol,legend = testTypes)
#legend("bottom",horiz = TRUE, pch = c(0,3,4),legend = c("80%","90%","95%"))
legend(.6,.511,yjust=.5,horiz = TRUE, pch = c(0,3,4),legend = c("80%","90%","95%"))
legend("bottomright",pch = 19,pt.cex=unique(aucs$effectSize),legend = unique(aucs$effectSize))


#### Plot Hits vs FA ####
ap <- aggregate(hits ~ criterion+sampleSize+effectSizes, data=allDt,mean)
ap2 <- aggregate(fa ~ criterion+sampleSize+effectSizes, data=allDt,mean)
hFA <- merge(ap,ap2,by=c("criterion","sampleSize","effectSizes"))
hFA$hits <- hFA$hits / numStudies
hFA$fa <- hFA$fa / numStudies
plot(hFA$fa,hFA$hits,col=hFA$criterion,pch=hFA$effectSizes*10,xlim=c(0,1),ylim=c(0,1),xlab="False Alarm Rate",ylab="Hit Rate")
legend("topright",pch=sort(unique(hFA$effectSizes))*10,legend = unique(hFA$effectSizes))
legend("bottom",horiz=T,pch=rep(19,length(hFA$criterion)),col=sort(unique(hFA$criterion)),legend = crits)
abline(a=1,b=-1,lty=2)


#plot false alarm for each effect size/sample size combo
if(1>2) { #plot initial stacked plots
for (i in 1:length(sSizesAll)) {
  #  readline(i)
  dr <- hFA[which(hFA$sampleSize == sSizesAll[i] & hFA$effectSizes == effSzAll[i]),c(1,5)]
  dr <- melt(dr,id = "criterion")
  
  print(ggplot(data = dr, aes(x = criterion, y = value, fill = variable)) + 
          geom_bar(stat = "identity") +
          ylab("Proportion of Outcomes\n") +
          scale_x_continuous(breaks=seq(1,8), labels=crits) +
          labs(title = paste("N=",sSizesAll[i],"d=",effSzAll[i])))
}

#plot hit vs false alarm for each effect size/sample size combo
for (i in 1:length(sSizesAll)) {
#  readline(i)
  dr <- hFA[which(hFA$sampleSize == sSizesAll[i] & hFA$effectSizes == effSzAll[i]),c(1,4,5)]
  dr <- melt(dr,id = "criterion")
  dr$value <- dr$value / 2

  print(ggplot(data = dr, aes(x = criterion, y = value, fill = variable)) + 
    geom_bar(stat = "identity") +
    ylab("Proportion of Outcomes\n") +
    scale_x_continuous(breaks=seq(1,8), labels=crits) +
    labs(title = paste("N=",sSizesAll[i],"d=",effSzAll[i])))
}

} #end initial stacked plots


#### Stacked plots of outcomes ####
if(1>2) {
#plot hit, FA, miss, corrRejection for each effect size/sample size combo
ap <- aggregate(hits ~ criterion+sampleSize+effectSizes, data=allDt,mean)
ap2 <- aggregate(fa ~ criterion+sampleSize+effectSizes, data=allDt,mean)
hFA <- merge(ap,ap2,by=c("criterion","sampleSize","effectSizes"))
ap2 <- aggregate(miss ~ criterion+sampleSize+effectSizes, data=allDt,mean)
hFA <- merge(hFA,ap2,by=c("criterion","sampleSize","effectSizes"))
ap2 <- aggregate(corrRej ~ criterion+sampleSize+effectSizes, data=allDt,mean)
hFA <- merge(hFA,ap2,by=c("criterion","sampleSize","effectSizes"))
hFA$hits <- hFA$hits / numStudies
hFA$fa <- hFA$fa / numStudies
hFA$corrRej <- hFA$corrRej / numStudies
hFA$miss <- hFA$miss / numStudies

#plot hit vs false alarm for each effect size/sample size combo
for (i in 1:length(sSizesAll)) {
  #  readline(i)
  dr <- hFA[which(hFA$sampleSize == sSizesAll[i] & hFA$effectSizes == effSzAll[i]),c(1,4,5,6,7)]
  dr <- melt(dr,id = "criterion")
  dr$value <- dr$value / 2
  
  print(ggplot(data = dr, aes(x = criterion, y = value, fill = variable)) + 
          geom_bar(stat = "identity") +
          ylab("Proportion of Outcomes\n") +
          xlab("Criterion for Statistical Significance\n") +
          scale_x_continuous(breaks=seq(1,8), labels=crits) +
          labs(title = paste("N=",sSizesAll[i],"d=",effSzAll[i])))

}
} #end if stacked outcomes


#### plot errors only (FA and misses, separately) ####




#### Plot ROC dist ####
aa <- aggregate(ROCdist ~ criterion, data=allDt,mean)
ab <- aggregate(ROCdist ~ criterion, data=allDt,sd)
ab$ci <- qnorm(.975) * ab$ROCdist / sqrt(numStudies)
titleText <- ifelse(length(sSizesAll) < 2, paste("N=",sSizesAll[1],"d=",effSzAll[1]), paste("firstES:",effSzAll[1]))
plot(seq(1:length(crits)),aa$ROCdist,bty="l",xaxt="n",xlab="Criterion for Statistical Significance",ylab="Distance to Perfection",pch=19,col=rainbow(length(crits)),cex=2,ylim=c(0,max(aa$ROCdist)+max(ab$ci)),main = titleText)
axis(side=1,at=seq(1:length(crits)),labels = crits)
for(i in 1:length(crits)) {
  segments(i,aa$ROCdist[i] - ab$ci[i],i,aa$ROCdist[i]+ab$ci[i])
}

#### plot ROC by effect size/sample size pair ####
if (1>2) {  
  par(mfrow=c(4,2))
  for (i in 1:length(sSizesAll)) {
    aa <- aggregate(ROCdist ~ criterion, data=allDt[which(allDt$sampleSize == sSizesAll[i] & allDt$effectSizes == effSzAll[i]),],mean)
    ab <- aggregate(ROCdist ~ criterion, data=allDt[which(allDt$sampleSize == sSizesAll[i] & allDt$effectSizes == effSzAll[i]),],sd)
    ab$ci <- qnorm(.975) * ab$ROCdist / sqrt(numStudies)
    titleText <- paste("N=",sSizesAll[i],"d=",effSzAll[i])
    plot(seq(1:length(crits)),aa$ROCdist,bty="l",xaxt="n",xlab="Criterion for Statistical Significance",ylab="Distance to Perfection",pch=19,col=rainbow(length(crits)),cex=2,ylim=c(0,max(aa$ROCdist)+max(ab$ci)),main = titleText)
    axis(side=1,at=seq(1:length(crits)),labels = crits)
    for(j in 1:length(crits)) {
      segments(j,aa$ROCdist[j] - ab$ci[j],j,aa$ROCdist[j]+ab$ci[j])
    }
  }
  par(mfrow=c(1,1))
}

#### plot dist2Perfection (Figure 10) ####
if(1>0) {
aa <- aggregate(ROCdist ~ criterion+effectSizes+power+testType,data=allDt,mean)
es <- sort(unique(aa$effectSizes))
pw <- sort(unique(ppower))
aa$eRank <- 0
aa$pRank <- 0
aa$critCat <- ifelse(aa$criterion<min(which(critCat==2)),1,2)
for(i in 1:length(effSzAll)) {
  aa$eRank[which(aa$effectSizes==es[i])] <- i
  aa$pRank[which(aa$power==pw[i])] <- i*.2+.4
}
aa$xRank <- aa$eRank + aa$pRank - 1 + (aa$critCat-1)*.1
aa$critCat <- ifelse(aa$criterion<min(which(critCat==2)),1,19)
ct <- ifelse(critCat<2,1,19)
xLabs <- c("Effect Size (d)","Effect Size (dz)","Effect Size (d)","Correlation (r)")
xx <- c(1,1,1,2)
par(mfrow=c(2,2))
for(i in 1:length(testTypes)) {
  ab <- aa[which(aa$testType==i),]
#  plot(ab$xRank,ab$ROCdist,col=rainbow(length(crits))[ab$criterion],pch=aa$critCat,cex=-4+(6*ab$power),ylim=c(0,.5),xaxt="n",bty="l",ylab="Distance to Perfection",xlab = "Effect size",main=testTypes[i])
  plot(ab$xRank,ab$ROCdist,col=rainbow(length(crits))[ab$criterion],pch=19,cex=-4+(6*ab$power),ylim=c(0,.6),xaxt="n",bty="l",ylab="Distance to Perfection",xlab = xLabs[i],main=testTypes[i])
  axis(side=1,at=seq(1,length(es))-.2,labels = es/xx[i])
#  if(i ==1) {
#    legend("top",horiz=T,pch=ct[1:4],col=rainbow(length(crits))[1:4],legend = crits[1:4])
#  } else if(i == 2) {
#    legend("top",horiz=T,pch=ct[5:8],col=rainbow(length(crits))[5:8],legend = crits[5:8])
#  }
}
par(mfrow=c(1,1))
plot(ab$xRank,ab$ROCdist,col=rainbow(length(crits))[ab$criterion],pch=19,cex=-4+(6*ab$power),ylim=c(0,.5),xaxt="n",bty="l",ylab="Distance to Perfection",xlab = "Effect size",main=testTypes[i],xlim=c(0,1.5))
legend("topright",pch=19,col=rainbow(length(crits)),legend=crits)
} #plot dist 2 perf


#### plot BF vs p-values ####
#plot(log(allPs),log(allBFs),bty="l")

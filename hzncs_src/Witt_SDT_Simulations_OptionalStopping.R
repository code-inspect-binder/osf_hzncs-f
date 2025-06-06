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
## Strategy: Simulate means from two groups and compute relevant statistics
## Mean for group 1 is .5 (as in 50% on a memory test)
## Mean for group 2 is a function of the specified effect size
##
## Two outcome measures:
##  (1) Area under the curve (AUC) - measure of discriminability
##  (2) Distance to perfection - Euclidean distance between a given point on the ROC curve and perfect performance (100% hits and 0% false alarms)
##
## NEW for optionalStopping (12/18/17):
## Only p<.05 and BF > 3
## If p-value is close but not there, run 20 more Ss
## Could run X more Ss; could define what close is
##



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
  allPs2 <- array(NA,c(length(sampleSizes),length(effectSizes),numStudies))
  allBFs <- array(0,c(length(sampleSizes),length(effectSizes),numStudies))
  allBFs2 <- array(NA,c(length(sampleSizes),length(effectSizes),numStudies))
  allESs <- array(0,c(length(sampleSizes),length(effectSizes),numStudies))
  allESs2 <- array(NA,c(length(sampleSizes),length(effectSizes),numStudies))
  allTs <- array(0,c(length(sampleSizes),length(effectSizes),numStudies))
  
  #Optional stopping rule params
  numExtraSs <- 10
  numTimesRunExtra <- 10
  
  for (i in 1:length(sampleSizes)) {
    N <- sampleSizes[i]

    for (j in 1:length(effectSizes)) {
      mean2 <- mean1 - (effectSizes[j] * sdPooled)

      for (k in 1:numStudies) {
        
        group1 <- rnorm(N,mean1,mySDs[1])
        group2 <- rnorm(N,mean2,mySDs[2])

        a <- t.test(group2, group1,var.equal=T)

        allPs[i,j,k] <- a$p.value
        allTs[i,j,k] <- a$statistic
        allBFs[i,j,k] <- exp(ttest.tstat(t=a$statistic, n1=N, n2=N, rscale = 0.707)[['bf']])
        allESs[i,j,k] <- (mean(group1) - mean(group2))/sqrt((sd(group1)^2 + sd(group2)^2)/2)
        BFcurr <- allBFs[i,j,k]

        #p-hack Optional Stopping Rule
        if (a$p.value < .2 & a$p.value > .05) {
          for(keepGoing in 1:numTimesRunExtra) {
            if (a$p.value > .05) {
              group1b <- c(group1,rnorm(numExtraSs,mean1,mySDs[1]))
              group2b <- c(group2,rnorm(numExtraSs,mean2,mySDs[2]))
              a <- t.test(group2b, group1b,var.equal=T)
            }
          }
          allPs2[i,j,k] <- a$p.value
          allESs2[i,j,k] <- (mean(group1b) - mean(group2b))/sqrt((sd(group1b)^2 + sd(group2b)^2)/2)
        }
        
        #BF-hack (run up to 4 * 20 more participants until significant but only if p < .20)
        if (BFcurr < 3 & BFcurr > 1) {
          for(keepGoing in 1:numTimesRunExtra) {
            if (BFcurr < 3) {
              group1b <- c(group1,rnorm(numExtraSs,mean1,mySDs[1]))
              group2b <- c(group2,rnorm(numExtraSs,mean2,mySDs[2]))
              a <- t.test(group2b, group1b,var.equal=T)
              BFcurr <- exp(ttest.tstat(t=a$statistic, n1=length(group1b), n2=length(group2b), rscale = 0.707)[['bf']])
            }
          }
          allBFs2[i,j,k] <- BFcurr
        }
        
        
        
      } #end for k numstudies  
    } #end for j effectSizes
  } #end for i sampleSizes
  
  return(list(allPs,allBFs,allESs, allPs2, allBFs2, allESs2))
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
  
  sg <- c(allPs[1,2,],allPs[1,1,])
  sf <- c(rep(1,numStudies),rep(0,numStudies))

  fit_glm <- glm(sf ~ sg, family=binomial(link="logit"))
  glm_response_scores <- predict(fit_glm, data.frame(sg,sf), type="response")

  
  sg2 <- c(allBFs[1,2,],allBFs[1,1,])
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
numStudies <- 10
numTimes <- 100

#Sample sizes and Effect sizes (note: sample size = sample size per group; assumes equal samples per group)
sSizesAll <- c(50,176,290,64,105,26,45)
effSzAll <- c(.3,.3,.3,.5,.5,.8,.8)
sSizesAll <- c(30)
effSzAll <- c(.5)
# for d=.5, n=64 = 80%; n=105; 95% power  [n = num per group for 2 independent samples t test]
# for d=.3, n=176 = 80%; n=290; 95% power
# for d=.8, n=26 = 80%; n=42; 95% power

#data frame to save SDT analyses
allDt <- as.data.frame(matrix(0,ncol=14,nrow=1))
colnames(allDt) <- c("sampleSize","effectSizes","criterion","hits","fa","actEffectSize","criterionF","miss","corrRej","numAmbigEff","numAmbigNull", "ROCdist","runNum","withHack")

#criteria for significance to be used
critCat <- c(1,1,1,2,2,2) # 1 = p-value; 2 = BF
critValue <- c(.1,.05,.01,2,3,5)  #critical value for significance,, aligns with critCat
crits <- c("p<.10","p<.05","p<.01","BF>2","BF>3","BF>5")  #label not needed for function

#Group mean and SDs
mean1 <- .5
mySDs <- c(.1, .1)

#AUC estimates
aucs <- as.data.frame(matrix(0,ncol=7,nrow=length(sSizesAll)*numTimes))
colnames (aucs) <- c("sampleSize","effectSize","runTime","aucP","aucBF","aucPhacked","aucBFhacked")
ai <- 0

# p-hacking estimates
numHacked <- as.data.frame(matrix(0,ncol=2,nrow=numTimes))
colnames(numHacked) <- c("pHack","BFhack")

###### Loops for Data simulation and SDT calculations ####
for (ssa in 1:length(sSizesAll)) {
  sampleSizes <- sSizesAll[ssa]
  effectSizes <- c(0,effSzAll[ssa])  #sets first effect size to 0; necessary for SDT analyses
  
  for(i in 1:numTimes) {
    
    ## Run simulations
    myOut <- runSims(sampleSizes,effectSizes, mean1, mySDs, numStudies)
    allPs <- myOut[[1]]
    allBFs <- myOut[[2]]
    allESs <- myOut[[3]]
    allPs2 <- myOut[[4]]
    allBFs2 <- myOut[[5]]
    allESs2 <- myOut[[6]]
    

    ## Run SDT analyses
    ##critical: Assumes the first effect size is 0
    dt <- evalSig(effectSizes,critCat,critValue,allPs,allBFs,allESs)
    dt$runNum <- i
    dt$withHack <- 0
    allDt <- rbind(allDt,dt)

    #SDT analyses with p- and BF-hacked outcomes
    notHacked <- which(is.na(allPs2))
    numHacked[i,1] <- length(which(!is.na(allPs2)))
    allPs2[notHacked] <- allPs[notHacked]
    
    notHacked <- which(is.na(allBFs2))
    numHacked[i,2] <- length(which(!is.na(allBFs2)))
    allBFs2[notHacked] <- allBFs[notHacked]
    
    notHacked <- which(is.na(allESs2))
    allESs2[notHacked] <- allESs[notHacked]
    
    dt <- evalSig(effectSizes,critCat,critValue,allPs2,allBFs2,allESs2)
    dt$runNum <- i
    dt$withHack <- 1
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
    plotIt <- FALSE
    currAUC <- plotAUC(allPs2,allBFs2,dt,numStudies,plotIt)  #currently does not save the AUCs
    aucs$aucPhacked[ai] <- currAUC[1]
    aucs$aucBFhacked[ai] <- currAUC[2]
    
    if (currAUC[1] != currAUC[2]) {
      print(paste(effectSizes,i))
      print(currAUC)
    }
  
  } #end for i runTimes
} #end for ssa

allDt <- allDt[-1,]  #get rid of initial dummy row that was just used to init allDt

######################
###### Analyses ######
######################


#### Compare AUCs ####
aucs$hackDiff <- aucs$aucPhacked - aucs$aucBFhacked
numHacked$hackDiff <- numHacked$pHack - numHacked$BFhack
plot(numHacked$hackDiff,aucs$hackDiff)
abline(lm(aucs$hackDiff ~ numHacked$hackDiff))

plot(aucs$aucP,aucs$aucBF,xlim=c(.6,1),ylim=c(.6,1),col=1,pch=19,bty="l",xlab="Area Under Curve for P-Value",ylab="AUC for Bayes factor")
points(aucs$aucPhacked,aucs$aucBFhacked,pch=1,col=2,cex=2)
points(aucs$aucPhacked[which(numHacked$hackDiff>0)],aucs$aucBFhacked[which(numHacked$hackDiff>0)],pch=19,col=4)
abline(a=0,b=1)
legend("topleft",pch=c(19,1,19),col=c(1,2,4),legend = c("No Hacking","More P hacks","More BF hacks"))


plot(aucs$aucP,aucs$aucPhacked,xlim=c(.5,1),ylim=c(.5,1),col=1,pch=19,bty="l",xlab="Area Under Curve for P-Value",ylab="AUC for Hacked P-Value")
abline(a=0,b=1)

plot(aucs$aucBF,aucs$aucBFhacked,xlim=c(.5,1),ylim=c(.5,1),col=1,pch=19,bty="l",xlab="Area Under Curve for Bayes Factor",ylab="AUC for Hacked Bayes Factor")
abline(a=0,b=1)


print(t.test(aucs$aucPhacked,aucs$aucBFhacked,paired=T,var.equal = T))
ap <- aggregate(aucP ~ sampleSize + effectSize,aucs,mean)
ap2 <- aggregate(aucBF ~ sampleSize + effectSize,aucs,mean)
auc <- merge(ap,ap2,by=c("sampleSize","effectSize"))
ap <- aggregate(aucPhacked ~ sampleSize + effectSize,aucs,mean)
auc <- merge(auc,ap,by=c("sampleSize","effectSize"))
ap2 <- aggregate(aucBFhacked ~ sampleSize + effectSize,aucs,mean)
auc <- merge(auc,ap2,by=c("sampleSize","effectSize"))
plot(seq(1:length(sSizesAll)),auc$aucP,col=1,bty="l",ylim=c(.7,1),pch=19,cex=2,xaxt="n",xlab="",ylab="Area Under the Curve")
points(seq(1:length(sSizesAll)),auc$aucBF,col=2,pch=19)
points(seq(1:length(sSizesAll)),auc$aucPhacked,col=3,pch=19)
points(seq(1:length(sSizesAll)),auc$aucBFhacked,col=4,pch=19)
legend("bottomleft",pch=c(1,19),legend = c("for p-values","for Bayes factors"))
axis(side=1,at=seq(1:length(sSizesAll)),auc$sampleSize)


##### ROC plots ####
ap <- aggregate(hits ~ criterion+sampleSize+effectSizes+withHack, data=allDt,mean)
ap2 <- aggregate(fa ~ criterion+sampleSize+effectSizes+withHack, data=allDt,mean)
hFA <- merge(ap,ap2,by=c("criterion","sampleSize","effectSizes","withHack"))
hFA$hits <- hFA$hits / numStudies
hFA$fa <- hFA$fa / numStudies
hFA$isP <- ifelse(critCat[hFA$criterion] == 1,2,4)
plot(hFA$fa,hFA$hits,col=rainbow(8,start=0,end=.8)[hFA$criterion],pch=hFA$withHack*18+1,xlim=c(0,1),ylim=c(0,1),cex=2,xlab="False Alarm Rate",ylab="Hit Rate")
legend("topright",pch=c(1,19,19,19),col=c(1,1,2,4),legend = c("No Hack","With Hack","p-value","Bayes factor"))
legend("bottom",horiz=T,pch=rep(19,length(hFA$criterion)),col=sort(unique(hFA$criterion)),legend = crits)
abline(a=1,b=-1,lty=2)

##### Stacked plots by SDT outcome ####
# Plot by whether hacked or not
#plot hit, FA, miss, corrRejection for each effect size/sample size combo
ap <- aggregate(hits ~ criterion+sampleSize+effectSizes+withHack, data=allDt,mean)
ap2 <- aggregate(fa ~ criterion+sampleSize+effectSizes+withHack, data=allDt,mean)
hFA <- merge(ap,ap2,by=c("criterion","sampleSize","effectSizes","withHack"))
ap2 <- aggregate(miss ~ criterion+sampleSize+effectSizes+withHack, data=allDt,mean)
hFA <- merge(hFA,ap2,by=c("criterion","sampleSize","effectSizes","withHack"))
ap2 <- aggregate(corrRej ~ criterion+sampleSize+effectSizes+withHack, data=allDt,mean)
hFA <- merge(hFA,ap2,by=c("criterion","sampleSize","effectSizes","withHack"))
hFA$hits <- hFA$hits / numStudies
hFA$fa <- hFA$fa / numStudies
hFA$corrRej <- hFA$corrRej / numStudies
hFA$miss <- hFA$miss / numStudies

#plot hit vs false alarm for each effect size/sample size combo
for (i in 1:length(sSizesAll)) {
  #  readline(i)
  dr <- hFA[which(hFA$sampleSize == sSizesAll[i] & hFA$effectSizes == effSzAll[i]),c(1,4,5,6,7,8)]
  dr <- melt(dr,id = c("criterion","withHack"))
  dr$value <- dr$value / 2
  dr$cond <- dr$criterion*2+dr$withHack - 1
  xLabs <- c("pNo","pHack","BFno","BFhack")
  xl <- 1
  for (j in 1:max(dr$criterion)) {
    xLabs[xl] <- paste(crits[j],"_no",sep="")
    xLabs[xl+1] <- paste(crits[j],"_hack",sep="")
    xl <- xl+2
  }
  print(ggplot(data = dr, aes(x = cond, y = value, fill = variable)) + 
          geom_bar(stat = "identity") +
          ylab("Proportion of Outcomes\n") +
          xlab("Criterion for Statistical Significance\n") +
          scale_x_continuous(breaks=seq(1,2*length(unique(dr$criterion))), labels=xLabs) +
          labs(title = paste("N=",sSizesAll[i],"d=",effSzAll[i])))
  
}


#### Plot ROC dist ####
aa <- aggregate(ROCdist ~ criterion, data=allDt,mean)
ab <- aggregate(ROCdist ~ criterion, data=allDt,sd)
ab$ci <- qnorm(.975) * ab$ROCdist / sqrt(numStudies)
plot(seq(1:length(crits)),aa$ROCdist,bty="l",xaxt="n",xlab="Criterion for Statistical Significance",ylab="Distance to Perfection",pch=19,col=seq(1:8),cex=2,ylim=c(0,max(aa$ROCdist)+max(ab$ci)))
axis(side=1,at=seq(1:length(crits)),labels = crits)
for(i in 1:length(crits)) {
  segments(i,aa$ROCdist[i] - ab$ci[i],i,aa$ROCdist[i]+ab$ci[i])
}


#### plot BF vs p-values ####
plot(log(allPs),log(allBFs),bty="l")


#### plot number of hacks ####
mycol <- rgb(0, 0, 255, max = 255, alpha = 125)
hist(numHacked$BFhack,xlim=c(0,10),main = "",col="black",xlab="Number of Studies Hacked")
hist(numHacked$pHack,add=T,col=mycol)
legend("topright",pch=c(15,15),pt.cex=2,col=c(1,mycol),legend = c("Bayes factor","P-value"))


#### Descriptives ####
print(apply(numHacked,2,mean)/20*100)
print(apply(numHacked,2,sd)/20*100)
print(apply(aucs,2,mean))


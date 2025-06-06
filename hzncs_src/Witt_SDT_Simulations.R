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



########## Packages #################
require(BayesFactor)
require(pROC)
require(ggplot2)
require(reshape2)



########## Functions #################

## Assumes 2 groups, so two in mySDs
runSims <- function(sampleSizes, effectSizes, mean1, mySDs, numStudiesNull,numStudiesEff) {
  
  sdPooled <- sqrt((mySDs[1]^2 + mySDs[2]^2) / 2)
  saveSims <- as.data.frame(matrix(0,ncol=7, nrow=length(sampleSizes) * (numStudiesNull+numStudiesEff)))
  colnames(saveSims) <- c("groupNum","modeledD","sampleSize","actualD","t","p","BF")
  sa <- 0
  
  for (i in 1:length(sampleSizes)) {
    N <- sampleSizes[i]

    for (j in 1:length(effectSizes)) {
      mean2 <- mean1 - (effectSizes[j] * sdPooled)
      
      numStudies3 <- ifelse(effectSizes[j] == 0, numStudiesNull,numStudiesEff)  #only run a quarter of real effects

      for (k in 1:numStudies3) {
        
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
        saveSims$modeledD[sa] <- effectSizes[j]
        saveSims$sampleSize[sa] <- N
        saveSims$actualD[sa] <- (mean(group1) - mean(group2))/sqrt((sd(group1)^2 + sd(group2)^2)/2)

        saveSims$p[sa] <- a$p.value
        saveSims$t[sa] <- a$statistic
        saveSims$BF[sa] <- exp(ttest.tstat(t=a$statistic, n1=N, n2=N, rscale = 0.707)[['bf']])

      } #end for k numstudies  
    } #end for j effectSizes
  } #end for i sampleSizes
  
  return(saveSims)
} #end function



#AUC Function
# Calculates the area under the curve (AUC) for a given set of simulations
plotAUC <- function(saveSims,dt,plotIt,crits,priorOdds) {

  auc <- rep(0,2+length(priorOdds))  #auc for p-value and Bayes factor
  
  sg <- c(log(saveSims$p[which(saveSims$modeledD > 0)]),log(saveSims$p[which(saveSims$modeledD == 0)]))
  sf <- c(rep(1,length(which(saveSims$modeledD > 0))),rep(0,length(which(saveSims$modeledD == 0))))

  fit_glm <- glm(sf ~ sg, family=binomial(link="logit"))
  glm_response_scores <- predict(fit_glm, data.frame(sg,sf), type="response")

  auc[1] <- auc(sf, glm_response_scores)
  
  sg3 <- c(saveSims$actualD[which(saveSims$modeledD > 0)],saveSims$actualD[which(saveSims$modeledD == 0)])
#  sg3 <- abs(c(saveSims$actualD[which(saveSims$modeledD > 0)],saveSims$actualD[which(saveSims$modeledD == 0)]))
  
  fit_glm3 <- glm(sf ~ sg3, family=binomial(link="logit"))
  glm_response_scores3 <- predict(fit_glm3, data.frame(sg3,sf), type="response")
  
  auc[2] <- auc(sf, glm_response_scores3)

  for(ii in 1:length(priorOdds)) {
    sg2 <- c(log(priorOdds[ii] * saveSims$BF[which(saveSims$modeledD > 0)]),log(priorOdds[ii] * saveSims$BF[which(saveSims$modeledD == 0)]))
    fit_glm2 <- glm(sf ~ sg2, family=binomial(link="logit"))
    glm_response_scores2 <- predict(fit_glm2, data.frame(sg2,sf), type="response")
    auc[2 + ii] <- auc(sf, glm_response_scores2)

}
  
  if(plotIt) {  
#    plot(roc(sf, glm_response_scores, direction="<"), col="green", lwd=4)
    plot.roc(sf, glm_response_scores, direction="<",legacy.axes = T, col="green", lwd=4,xlab="False Alarm Rate",ylab="Hit Rate",auc.polygon=T)
    lines(roc(sf, glm_response_scores2, direction="<"), col="blue", lwd=1)
    #lines(roc(sf, glm_response_scores3, direction="<"), col="red", lwd=.5)  #if want to include ROC by effect size
    
    rCol <- rainbow(length(crits))
    dtCol <- rCol[dt$criterion]
    pSz <- ifelse(dt$critCat==1,1.5,.8)+1
    pCH <- ifelse(dt$critCat==1,1,19)
    points(1-(dt$fa/numStudiesNull),dt$hits/numStudiesEff,col=dtCol,pch=19,cex=pSz)
    legend("bottomright",col=rCol,legend = crits,pch=19)
  } #end if plotIt
  
  return(auc)
}  # end function plotAUC



plotAUCBFonly <- function(saveSims,dt,plotIt,crits,priorOdds) {
  
  auc <- rep(0,length(priorOdds))  #auc for BF 50/50; BF 1/10; BF 10/1
  dt <- dt[which(dt$critCat == 2),]
  dt$hits <- dt$hits / dt$numStudiesEff
  dt$fa <- dt$fa / dt$numStudiesNull
  
  sf <- c(rep(1,length(which(saveSims$modeledD > 0))),rep(0,length(which(saveSims$modeledD == 0))))

  for(i in 1:length(priorOdds)) {
    sg2 <- c(log(priorOdds[i] * saveSims$BF[which(saveSims$modeledD > 0)]),log(priorOdds[i] * saveSims$BF[which(saveSims$modeledD == 0)]))
    fit_glm2 <- glm(sf ~ sg2, family=binomial(link="logit"))
    glm_response_scores2 <- predict(fit_glm2, data.frame(sg2,sf), type="response")
    auc[i] <- auc(sf, glm_response_scores2)
      
    if(i == 1) {
      plot.roc(sf, glm_response_scores2, direction="<",legacy.axes = T, col=rainbow(length(priorOdds))[i], lwd=length(priorOdds)-i,xlab="False Alarm Rate",ylab="Hit Rate",auc.polygon=T)
    } else {
      lines(roc(sf, glm_response_scores2, direction="<"), col=rainbow(length(priorOdds))[i], lwd=length(priorOdds)-i)  
    }
  }
    mm <- aggregate(hits ~ criterion + priorOdds, dt, mean)
    m2 <- aggregate(fa ~ criterion + priorOdds, dt, mean)
    mm$col <- "black"
    mm$cex <- 1

    rCol <- rainbow(length(priorOdds))
    for (i in 1:length(priorOdds)) {
      mm$col[which(mm$priorOdds == priorOdds[i])] <- rCol[i]
      mm$cex[which(mm$priorOdds == priorOdds[i])] <- 1 + ((i*4)/10)
    }
    points(1-(m2$fa),mm$hits,col=mm$col,pch=19,cex=mm$cex)
    
    #legend("bottomright",lty=1,col=rainbow(length(priorOdds)),legend = priorOdds)  #lines only
    legend("bottomright",pch=19,col=rainbow(length(priorOdds)),legend = priorOdds)
    
  return(auc)
}  # end function plotAUCBFonly



## sample call:
# dt <- evalSig(effectSizes,critCat,critValue,allPs,allBFs,allESs)

evalSig <- function(critCat,critValue,saveSims,priorOdds) {

  effectSizes <- sort(unique(saveSims$modeledD))
  sampleSizes <- sort(unique(saveSims$sampleSize))
  numEffectSizes <- length(effectSizes) - 1
  numCriteria <- length(critCat)
  nc2 <- length(critCat[which(critCat == 1)]) + (length(critCat[which(critCat == 2)]) * length(priorOdds))
  allC <- length(sampleSizes) * numEffectSizes * nc2 
  dt <- as.data.frame(matrix(NA,ncol=14,nrow=allC))
  colnames(dt) <- c("sampleSize","effectSizes","criterion","hits","fa","actEffectSize","criterionF","miss","corrRej","numAmbigEff","numAmbigNull", "ROCdist","critCat","priorOdds")
  da <- 0  #index for dt

  for(i in 1:length(sampleSizes)) {
    for(j in 1:numEffectSizes) {
      for (k in 1:numCriteria) {
      
        jj <- j+1  #assumes first effect size is 0, working with next one
        
        if (critCat[k] == 1) {  #p-value
          da <- da+1
          dt$sampleSize[da] <- sampleSizes[i]
          dt$effectSizes[da] <- effectSizes[j+1]
          dt$criterion[da] <- k

          dt$hits[da] <- length(which(saveSims$p[which(saveSims$modeledD>0)] <= critValue[k]))
          dt$fa[da] <-length(which(saveSims$p[which(saveSims$modeledD==0)] <= critValue[k]))
          dt$miss[da] <- length(which(saveSims$modeledD > 0)) - dt$hits[da]
          dt$corrRej[da] <- length(which(saveSims$modeledD == 0)) - dt$fa[da]
          dt$criterionF[da] <- paste("p",critValue[k])
          dt$critCat[da] <- 1
        } else {  #Bayes factor
          for (jj in 1:length(priorOdds)) {
            da <- da+1
            dt$sampleSize[da] <- sampleSizes[i]
            dt$effectSizes[da] <- effectSizes[j+1]
            dt$criterion[da] <- k
            
            dt$hits[da] <- length(which(saveSims$BF[which(saveSims$modeledD > 0)]*priorOdds[jj] >= critValue[k]))
            dt$fa[da] <- length(which(saveSims$BF[which(saveSims$modeledD == 0)]*priorOdds[jj] >= critValue[k]))
            dt$miss[da] <- length(which(saveSims$BF[which(saveSims$modeledD > 0)]*priorOdds[jj] <= (1/critValue[k])))
            dt$corrRej[da] <- length(which(saveSims$BF[which(saveSims$modeledD == 0)*priorOdds[jj]] <= (1/critValue[k])))
            dt$numAmbigEff[da] <- length(which(saveSims$BF[which(saveSims$modeledD > 0)]*priorOdds[jj] < critValue[k] & saveSims$BF[which(saveSims$modeledD > 0)]*priorOdds[jj] > (1/critValue[k]) ))
            dt$numAmbigNull[da] <- length(which(saveSims$BF[which(saveSims$modeledD == 0)]*priorOdds[jj] < critValue[k] & saveSims$BF[which(saveSims$modeledD == 0)]*priorOdds[jj] > (1/critValue[k]) ))
            dt$criterionF[da] <- paste("BF",critValue[k])
            dt$critCat[da] <- 2
            dt$priorOdds[da] <- priorOdds[jj]
          }
        } 
      
      h1 <- 1 - (dt$hits[da]/length(which(saveSims$modeledD > 0)))  #vertical distance
      h2 <- dt$fa[da]/length(which(saveSims$modeledD == 0))  #horizDistance
      dt$ROCdist[da] <- sqrt(h1^2 + h2^2)
      
      dt$actEffectSize[da] <- mean(saveSims$actualD[which(saveSims$modeledD > 0)])
      
      }
   }
  }
  
  return(dt)
}



############################
####### Start Code #########
############################


#Num of studies to be run and number of times that number of studies should be run
numStudiesNull <- 10 #number of studies to run with d = 0
numStudiesEff <- 10  #number of studies to run with d > 0
numTimes <- 100  #number of times to run each set of studies


#If running replications
runReplication <- FALSE

#prior Odds: for use when calculating posterior odds; keep as 1 for using Bayes factor only
priorOdds <- c(1)  
priorOdds <- sort(priorOdds, decreasing = T) #useful for plotting but not necessary

#Sample sizes and Effect sizes (note: sample size = sample size per group; assumes equal samples per group)
sSizesAll <- c(64)
effSzAll <- c(.5)


#data frame to save sims outcomes
saveSims <- as.data.frame(matrix(0,ncol=7, nrow=1))
colnames(saveSims) <- c("groupNum","modeledD","sampleSize","actualD","t","p","BF")

#data frame to save SDT analyses
allDt <- as.data.frame(matrix(0,ncol=15,nrow=1))
colnames(allDt) <- c("sampleSize","effectSizes","criterion","hits","fa","actEffectSize","criterionF","miss","corrRej","numAmbigEff","numAmbigNull", "ROCdist","critCat","priorOdds","runNum")

#criteria for significance to be used
critCat <- c(1,1,1,1,2,2,2,2) # 1 = p-value; 2 = BF
critValue <- c(.1,.05,.005,.001,1,2,3,10)  #critical value for significance,, aligns with critCat
crits <- c("p<.10","p<.05","p<.005","p<.001","BF>1","BF>2","BF>3","BF>10")  #label not needed for function

#Group mean and SDs
mean1 <- .5  #mean for group 1 (e.g. 50% correct on memory test)
mySDs <- c(.1, .1) #standard deviation for both groups

#AUC estimates
if (length(priorOdds) == 1) {
  aa <- 6
  ab <- "aucBF1"
} else {
  aa <- 5+length(priorOdds)
  ab <- rep("BF",length(priorOdds))
  for(i in 1:length(priorOdds)) {
    ab[i] <- paste("aucBF",priorOdds[i],sep="")
  }
}
aucs <- as.data.frame(matrix(0,ncol=aa,nrow=length(sSizesAll)*numTimes))
colnames(aucs) <- c("sampleSize","effectSize","runTime","aucP","aucES",ab)
ai <- 0


###### Loops for Data simulation and SDT calculations ####
for (ssa in 1:length(sSizesAll)) {
  sampleSizes <- sSizesAll[ssa]  #can also run all at once, although code needs to be changed during graphing/analysis
  effectSizes <- c(0,effSzAll[ssa])  #sets first effect size to 0; necessary for SDT analyses
  
  for(i in 1:numTimes) {
    
    ## Run simulations
    ss <- runSims(sampleSizes,effectSizes, mean1, mySDs, numStudiesNull, numStudiesEff)

    if (runReplication) {
      #run again (e.g. run replications)
      ss2 <- runSims(sampleSizes,effectSizes, mean1, mySDs, numStudiesNull, numStudiesEff)
      ss$p <- pmax(ss$p,ss2$p)
      ss$BF <- pmin(ss$BF,ss2$BF)
      ss$actualD <- pmin(ss$actualD,ss2$actualD)
    }
    
  
    ss$groupNum <- i
    saveSims <- rbind(saveSims,ss)
    
    ## Run SDT analyses
    ##critical: Assumes the first effect size is 0
    dt <- evalSig(critCat,critValue,ss,priorOdds)
    dt$runNum <- i
    allDt <- rbind(allDt,dt)
    
    # Plot ROC curves and caculate AUC (if desired)
    plotIt <- FALSE
    currAUC <- plotAUC(ss,dt,plotIt,crits,priorOdds)  
    ai <- ai+1
    aucs[ai,1] <- sampleSizes[1]
    aucs[ai,2] <- effectSizes[2]
    aucs[ai,3] <- i #runNumber
    aucs[ai,4] <- currAUC[1]
    aucs[ai,5] <- currAUC[2]
    for(ii in 1:length(priorOdds)) {
      aucs[ai,5 + ii] <- currAUC[2 + ii]
    }    
  
  } #end for i runTimes
} #end for ssa

allDt <- allDt[-1,]  #get rid of initial dummy row that was just used to init allDt
saveSims <- saveSims[-1,]  #get rid of initial dummy row that was just used to init allDt

######################
###### Analyses ######
######################

####Compare AUCs ####
plot(aucs$aucP,aucs$aucBF1,col=aucs$effectSize*10,cex=1.5,bty="l",ylim=c(.5,1),xlab="Area Under Curve for P-Value",ylab="AUC for Bayes factor or Effect Size")
if(length(priorOdds) > 1) {
  for(i in 1:length(priorOdds)) {
    points(aucs$aucP,aucs[,i+5],pch=19)
  }
}

#### AUCs for Cohen's d ####
plot(aucs$aucP,aucs$aucES,pch=1,cex=1.5,bty="l",main=NA,ylim=c(.5,1),xlim=c(.5,1),xlab="AUC for P-Value and Bayes Factor",ylab="AUC for Cohen's d")
abline(a=0,b=1,lty=2)

#### Plot hits by FA ####
ap <- aggregate(hits ~ criterion+sampleSize+effectSizes, data=allDt,mean)
ap2 <- aggregate(fa ~ criterion+sampleSize+effectSizes, data=allDt,mean)
hFA <- merge(ap,ap2,by=c("criterion","sampleSize","effectSizes"))
hFA$hits <- hFA$hits / numStudiesEff
hFA$fa <- hFA$fa / numStudiesNull
pCh <- ifelse(critCat[hFA$criterion]==1,1,19)
pSz <- ifelse(critCat[hFA$criterion]==1,1.5,1)+1
plot(hFA$fa,hFA$hits,col=rainbow(length(crits)),cex=pSz,cex.lab=1.2,pch=pCh,xlim=c(0,1),ylim=c(0,1),xlab="False Alarm Rate",ylab="Hit Rate")
legend("bottomright",horiz=F,cex=1.2,pch=pCh,col=rainbow(length(crits)),legend = crits)
points(0,1,pch=3)


#### Comp AUCs across prior odds ####
if(length(priorOdds) > 1) {
  allDt$numStudiesEff <- numStudiesEff
  allDt$numStudiesNull <- numStudiesNull
  compAucs <- plotAUCBFonly(saveSims,allDt,TRUE,crits,priorOdds)
}

#### Plot hit, FA, miss, corrRejection for each effect size/sample size combo ####
ap <- aggregate(hits ~ criterion+sampleSize+effectSizes, data=allDt,mean)
ap2 <- aggregate(fa ~ criterion+sampleSize+effectSizes, data=allDt,mean)
hFA <- merge(ap,ap2,by=c("criterion","sampleSize","effectSizes"))
ap2 <- aggregate(miss ~ criterion+sampleSize+effectSizes, data=allDt,mean)
hFA <- merge(hFA,ap2,by=c("criterion","sampleSize","effectSizes"))
ap2 <- aggregate(corrRej ~ criterion+sampleSize+effectSizes, data=allDt,mean)
hFA <- merge(hFA,ap2,by=c("criterion","sampleSize","effectSizes"))
hFA$hits <- hFA$hits / numStudiesEff
hFA$fa <- hFA$fa / numStudiesNull
hFA$corrRej <- hFA$corrRej / numStudiesNull
hFA$miss <- hFA$miss / numStudiesEff

for (i in 1:length(sSizesAll)) {
  #  readline(i)
  dr <- hFA[which(hFA$sampleSize == sSizesAll[i] & hFA$effectSizes == effSzAll[i]),c(1,4,5,6,7)]
  dr <- melt(dr,id = "criterion")
  dr$value <- dr$value / 2
  
  print(ggplot(data = dr, aes(x = criterion, y = value, fill = variable)) + 
          geom_bar(stat = "identity") +
          ylab("Proportion of Outcomes\n") +
          xlab("Criterion for Statistical Significance\n") +
          scale_x_continuous(breaks=seq(1,length(crits)), labels=crits) +
          labs(title = paste("N=",sSizesAll[i],"d=",effSzAll[i])))

}


#### Plot distance to perfection ####
aa <- aggregate(ROCdist ~ criterion, data=allDt,mean)
ab <- aggregate(ROCdist ~ criterion, data=allDt,sd)
ab$ci <- qnorm(.975) * ab$ROCdist / sqrt(length(ab$ROCdist))
titleText <- ifelse(length(sSizesAll) < 2, paste("N=",sSizesAll[1],"d=",effSzAll[1]), paste("firstES:",effSzAll[1]))
plot(seq(1:length(crits)),aa$ROCdist,bty="l",xaxt="n",xlab="Criterion for Statistical Significance",ylab="Distance to Perfection",pch=19,col=rainbow(length(crits)),cex=2,ylim=c(0,max(aa$ROCdist)+max(ab$ci)),main=titleText)
axis(side=1,at=seq(1:length(crits)),labels = crits)
for(i in 1:length(crits)) {
  segments(i,aa$ROCdist[i] - ab$ci[i],i,aa$ROCdist[i]+ab$ci[i])
}


#### plot BF vs p-values ####
plot(log(saveSims$p),log(saveSims$BF),bty="l",xlab="log(p-value)", ylab="log(Bayes factor)")


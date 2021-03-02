library(R.matlab)
library(car)
source('C:/Users/daumail/OneDrive - Vanderbilt/Documents/R/useful_functions/rowMaxs.R')
allData = readMat("C:/Users/daumail/OneDrive - Vanderbilt/Documents/LGN_data/single_units/binocular_adaptation/all_units/all_orig_bs_zscore_trials.mat")

pvals <- array(0, c(length(allData$peak.aligned.trials),6))
r2 <- array(0, c(length(allData$peak.aligned.trials),3))



for(i in 1:length(allData$peak.aligned.trials)){
  #dim 1 = cell number
  #dim 2 = data type: 1 = origin, 2 =bs origin, 3 = z-score
  #dim 3 = mono or binocular stim: 1 =monocular, 2 = binocular
  monoOrigin = allData$peak.aligned.trials[[i]][[1]][[1]]
  binoOrigin = allData$peak.aligned.trials[[i]][[1]][[2]]
  
  monopk1<- rowMaxs(t(monoOrigin[[1]]))
  monopk2<- rowMaxs(t(monoOrigin[[2]]))
  monopk3<- rowMaxs(t(monoOrigin[[3]]))
  monopk4<- rowMaxs(t(monoOrigin[[4]]))
  
  binopk1<- rowMaxs(t(binoOrigin[[1]]))
  binopk2<- rowMaxs(t(binoOrigin[[2]]))
  binopk3<- rowMaxs(t(binoOrigin[[3]]))
  binopk4<- rowMaxs(t(binoOrigin[[4]]))
  

col_bumps_mono <- cbind(monopk1,monopk2,monopk3,monopk4)
col_bumps_bino <- cbind(binopk1,binopk2,binopk3,binopk4)

col_bumps_all <- rbind(col_bumps_mono,col_bumps_bino)
bumps_df.split = data.frame(col_bumps_all)
colnames(bumps_df.split) = c("Pk1","Pk2","Pk3","Pk4")

#create trial index
#monotrIdx = t(t(1:length(monoOrigin[[1]][1,])))
#binotrIdx = t(t(1:length(binoOrigin[[1]][1,])))+length(monoOrigin[[1]][1,])
#alltrIdx = rbind(monotrIdx,binotrIdx)
#bumps_df.split$Trial <- alltrIdx

#condition label
monoCond <- c("Monocular")
monoLabel<- t(t(rep(monoCond, each =length(monoOrigin[[1]][1,]))))
binoCond <- c("Binocular")
binoLabel<- t(t(rep(binoCond, each =length(binoOrigin[[1]][1,]))))

condLabel <- rbind(monoLabel, binoLabel)
bumps_df.split$Condition <- condLabel
bumps_df.split <- bumps_df.split[,c("Condition","Pk1","Pk2","Pk3","Pk4")]

##Omnibus test
bumps_df.split$Condition <- factor(bumps_df.split$Condition)
#attach(bumps_df.split)

mult.fit<-lm(as.matrix(bumps_df.split[,2:5])~bumps_df.split$Condition)
library(car)

measures<-as.factor(1:4)
mult.ANOVA<- Anova(mult.fit, idata = data.frame(measures),idesign = ~measures, type = "III")
summary.mult.ANOVA <-summary(mult.ANOVA)
## CONTRASTS TESTS ## MAIN EFFECTS ##SIMPLE EFFECTS ##INTERACTIONS


Average<-rowMeans(bumps_df.split[,2:5])
#tapply(Average,Condition,FUN = mean)

##contrast in the condition MAIN EFFECT
Average.aov<-aov(Average~bumps_df.split$Condition)
conditionMain <- summary(Average.aov)
pvals[i,1] <- conditionMain[[1]][1,5]
#TukeyHSD(Average.aov) we don't need it here as we only have two groups to compare (monocular/binocular)

##test user defined contrasts (doesn't work here as only two levels in Condition)
library(gmodels)
#levels(Condition)

#cmat<-rbind(" Mono vs Bino" = c(1,-1),
#            "Bino vs mono" = c(-1,1))
#fit.contrast(Average.aov,Condition,cmat,df=TRUE)

#SIMPLE EFFECT of condition and contrast in Pk
#Pk1.aov<-aov(Pk1~Condition)
#summary(Pk1.aov)
#TukeyHSD(Pk1.aov) we don't need it here as we only have two groups to compare (monocular/binocular)

##Within subject effect

#1.First part of the analysis: test linear trend and interaction with condition 
#MAIN EFFECT of Linear Trend (just to make sure that most linear trends that show interaction are significant linear trends)
linear<- -3*bumps_df.split[,2]-bumps_df.split[,3]+bumps_df.split[,4]+3*bumps_df.split[,5]
#Condition <- factor(Condition, levels = c("Monocular","Binocular"))
options(contrasts = c("contr.sum","contr.poly"))
peakMain <- coef(summary(lm(linear~bumps_df.split$Condition)))
#tapply(linear,Condition,FUN=mean)
#PK main effect
pvals[i,2] <- peakMain[1,4]

##Interaction contrast (most important test of the script)
aov.Interaction <- aov(linear~bumps_df.split$Condition)
pvals[i,6] <- summary(aov.Interaction)[[1]][1,5]
#TukeyHSD(aov(linear~Condition))

#2. Second part of the analysis: How much does the linear trend explain the variance together with quadratic and cubic trends?
library(reshape2)
long_dat <- melt(bumps_df.split, id.var=c('Condition'), variable.name = "Peak", value.name = "Response")
long_dat.mono <- long_dat[long_dat$Condition == 'Monocular',]
long_dat.bino <- long_dat[long_dat$Condition == 'Binocular',]

contrasts(long_dat$Peak) <- contr.poly(4)
long_datX <- model.matrix(~1+long_dat$Peak, data = long_dat)
long_dat[,c("cLin","cQuad","cCub")] <- long_datX[,2:4]
regressX <-lm(Response~1+cLin+cQuad+cCub, data = long_dat)
aovModel <-anova(regressX)
pvals[i,3] <-aovModel[1,"Pr(>F)"]
pvals[i,4] <-aovModel[2,"Pr(>F)"]
pvals[i,5] <-aovModel[3,"Pr(>F)"]

SumSq <- aovModel[1:3,"Sum Sq"]
names(SumSq) <- c("cLinr","cQuad","cCub")
r2[i,] <- round(SumSq / sum(SumSq),2)



 

##simple effect of Pk
#monocular.fit<-lm(as.matrix(bumps_df.split[Condition=="Monocular",2:5])~1)

#monocular.ANOVA<- Anova(monocular.fit, idata = data.frame(measures),idesign = ~measures)
#summary.monocular.ANOVA<-summary(monocular.ANOVA)
#summary.monocular.ANOVA$univariate.tests
#summary.monocular.ANOVA$pval.adjustments


#binocular.fit<-lm(as.matrix(bumps_df.split[Condition=="Binocular",2:5])~1)
#binocular.ANOVA<- Anova(binocular.fit, idata = data.frame(measures),idesign = ~measures)
#summary.binocular.ANOVA<-summary(binocular.ANOVA)
#summary.binocular.ANOVA$univariate.tests
#summary.binocular.ANOVA$pval.adjustments



#summary(mult.ANOVA)$univariate.tests
#F<- (summary(mult.ANOVA)$univariate.tests[3,1]/summary(mult.ANOVA)$univariate.tests[3,2])/(summary(mult.ANOVA)$univariate.tests[3,3]/summary(mult.ANOVA)$univariate.tests[3,4]) 
#F
#summary(mult.ANOVA)$pval.adjustments
#qf(df1 = summary(mult.ANOVA)$univariate.tests[3,2]*summary(mult.ANOVA)$pval.adjustments[1,3],df2 = summary(mult.ANOVA)$univariate.tests[3,4]*summary(mult.ANOVA)$pval.adjustments[1,3],p=1-0.05/2)



}
pvals_df <- data.frame(pvals)
colnames(pvals_df) = c("MainCondition","MainTrend","Lin","Quad","Cub","Interaction")

filename <- paste("C:/Users/daumail/Documents/LGN_data/single_units/binocular_adaptation/all_units/mixedmodel_pvals_anova_linearTrend", ".csv", sep = "")
write.csv(pvals_df, file = filename)

filename <- paste("C:/Users/daumail/Documents/LGN_data/single_units/binocular_adaptation/all_units/mixedmodel_r2_anova_linearTrend", ".csv", sep = "")
write.csv(r2, file = filename)

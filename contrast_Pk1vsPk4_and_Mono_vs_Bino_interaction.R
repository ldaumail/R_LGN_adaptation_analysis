library(R.matlab)
library(car)
source('C:/Users/daumail/Documents/R/useful_functions/rowMaxs.R')
allData = readMat("C:/Users/daumail/Documents/LGN_data/single_units/binocular_adaptation/all_units/all_orig_bs_zscore_trials.mat")

stats <- array(0, c(length(allData$peak.aligned.trials),2))

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
  monotrIdx = t(t(1:length(monoOrigin[[1]][1,])))
  binotrIdx = t(t(1:length(binoOrigin[[1]][1,])))+length(monoOrigin[[1]][1,])
  alltrIdx = rbind(monotrIdx,binotrIdx)
  bumps_df.split$Trial <- alltrIdx
  
  #condition label
  monoCond <- c("Monocular")
  monoLabel<- t(t(rep(monoCond, each =length(monoOrigin[[1]][1,]))))
  binoCond <- c("Binocular")
  binoLabel<- t(t(rep(binoCond, each =length(binoOrigin[[1]][1,]))))
  
  condLabel <- rbind(monoLabel, binoLabel)
  bumps_df.split$Condition <- condLabel
  bumps_df.split <- bumps_df.split[,c("Trial","Condition","Pk1","Pk2","Pk3","Pk4")]
  
  attach(bumps_df.split)
  
  ##Interaction contrast
  #analysis <- aov( I(Pk1 - Pk4) ~ Condition + Error( Trial), data =bumps_df.split)
  #pvals[i] <- summary(analysis)[[2]][[1]][1,"Pr(>F)"]
 # analysis <- aov( I(Pk1 - Pk4) ~ Condition, data =bumps_df.split)
#  pvals[i] <- summary(analysis)[[1]][1,"Pr(>F)"]
#  mult.fit<-lm(as.matrix(bumps_df.split[,3:6])~Condition)
#  measures<-as.factor(1:4)
#  mult.ANOVA<- Anova(mult.fit, idata = data.frame(measures),idesign = ~measures, type = "III")
#  summary.mult.ANOVA <-summary(mult.ANOVA)
#  pvals[i] <-summary.mult.ANOVA$pval.adjustments[2,"Pr(>F[GG])"]
  #Compute difference and perform independent samples t-test
  
  x = bumps_df.split$Pk1[bumps_df.split$Condition == 'Monocular'] -bumps_df.split$Pk4[bumps_df.split$Condition == 'Monocular']
  y = bumps_df.split$Pk1[bumps_df.split$Condition == 'Binocular'] -bumps_df.split$Pk4[bumps_df.split$Condition == 'Binocular']
  
 indepSampTest <- t.test(x, y, alternative = "two.sided", var.equal = FALSE)
  stats[i,1] <-indepSampTest$statistic
  stats[i,2] <-indepSampTest$p.value
}
stats_df <- data.frame(stats)
colnames(stats_df) = c("T stat Pk1vsPk4 and MonovsBino", "pvalue")
stats_df

filename <- paste("C:/Users/daumail/Documents/LGN_data/single_units/binocular_adaptation/all_units/interaction_contrast_indepSampleTtest", ".csv", sep = "")
write.csv(stats_df, file = filename)

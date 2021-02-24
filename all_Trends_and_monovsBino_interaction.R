library(R.matlab)
library(car)
library(magrittr)
library(MASS)
library(dplyr)
source('C:/Users/daumail/OneDrive - Vanderbilt/Documents/R/useful_functions/rowMaxs.R')
allData = readMat("C:/Users/daumail/OneDrive - Vanderbilt/Documents/LGN_data/single_units/binocular_adaptation/all_units/all_orig_bs_zscore_trials.mat")

pvals <- array(0, c(length(allData$peak.aligned.trials),7))
r2 <- array(0, c(length(allData$peak.aligned.trials),7))

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


  ## CONTRASTS TESTS ## MAIN EFFECTS ##SIMPLE EFFECTS ##INTERACTIONS
  library(gmodels)

  #Is the linear trend enough to explain the group difference?
  library(reshape2)
  long_dat <- melt(bumps_df.split, id.var=c('Condition'), variable.name = "Peak", value.name = "Response")
  
  contrasts(long_dat$Peak) <- contr.poly(4)
  long_dat$Condition <- factor(long_dat$Condition, levels = c("Monocular","Binocular"))
  contrasts(long_dat$Condition) <- contr.treatment(2)
  long_datX <- model.matrix(~1+Peak*Condition, data = long_dat)
  long_dat[,c("cLin","cQuad","cCub","Condition", "cLin:Condition","cQuad:Condition","cCub:Condition")] <- long_datX[,2:8]
  regressX <-lm(Response~1+(cLin+cQuad+cCub)*Condition, data = long_dat)
  aovModel <-anova(regressX)
  pvals[i,1] <-aovModel[1,"Pr(>F)"]
  pvals[i,2] <-aovModel[2,"Pr(>F)"]
  pvals[i,3] <-aovModel[3,"Pr(>F)"]
  pvals[i,4] <-aovModel[4,"Pr(>F)"]
  pvals[i,5] <-aovModel[5,"Pr(>F)"]
  pvals[i,6] <-aovModel[6,"Pr(>F)"]
  pvals[i,7] <-aovModel[7,"Pr(>F)"]
  SumSq <- aovModel[1:7,"Sum Sq"]
  names(SumSq) <- c("cLinr","cQuad","cCub","Condition", "cLin:Condition","cQuad:Condition","cCub:Condition")
  r2[i,] <- round(SumSq / sum(SumSq),2)
  
  #get the contrast matrix
  Xc <- long_dat %>%
    group_by(Peak, Condition) %>%
    summarise() %>%
    model.matrix(~1+Peak*Condition, .) %>%
    as.data.frame() %>% as.matrix()
   rownames(Xc) <- c("Pk1_Mono","Pk1_Bino", "Pk2_Mono","Pk2_Bino", "Pk3_Mono", "Pk3_Bino", "Pk4_Mono", "Pk4_Bino" )

}
pvals_df <- data.frame(pvals)
colnames(pvals_df) = c("Lin","Quad","Cub","Condition","cLin:Condition","cQuad:Condition","cCub:Condition")

filename <- paste("C:/Users/daumail/Documents/LGN_data/single_units/binocular_adaptation/all_units/linearmodel_pvals_allTrends_02242021", ".csv", sep = "")
write.csv(pvals_df, file = filename)

filename <- paste("C:/Users/daumail/Documents/LGN_data/single_units/binocular_adaptation/all_units/linearmodel_r2_anova_allTrends_02242021", ".csv", sep = "")
write.csv(r2, file = filename)
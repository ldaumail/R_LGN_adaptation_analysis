library(R.matlab)
library(lmerTest)
library(car)
library(lme4)
#library(pbkrtest)
library(emmeans)
library(pander)
source('C:/Users/daumail/OneDrive - Vanderbilt/Documents/R/useful_functions/get_ddf_Lb.R')
source('C:/Users/daumail/OneDrive - Vanderbilt/Documents/R/useful_functions/rowMaxs.R')
#library(BayesFactor)

##
## create a for loop to treat all the data at once

#file names of all 71 units in the folder
#ignore for the data analysis
folderfilenames <-list.files("C:/Users/daumail/OneDrive - Vanderbilt/Documents/LGN_data_042021/single_units/inverted_power_channels/good_single_units_data_4bumps_more/new_peak_alignment_anal/su_peaks_03032020_corrected/orig_peak_values/", pattern=NULL, all.files=FALSE,  full.names=FALSE)#ignore
#ignore for the data analysis
folderfilenames <- folderfilenames[1:71] #ignore for the data analysis

#file names of all 71 units with cell class label
filenames <- read.csv("C:/Users/daumail/OneDrive - Vanderbilt/Documents/LGN_data_042021/single_units/inverted_power_channels/good_single_units_data_4bumps_more/new_peak_alignment_anal/su_peaks_03032020_corrected/orig_peak_values/all_units/filenames_layers.csv")
library(tidyr)
#keep only labelled units file names
row.has.na <- apply(filenames, 1, function(x){any(is.na(x))})
filenames.filt <- filenames[!row.has.na,]


#data_filename_list <- Sys.glob("C:/Users/daumail/OneDrive - Vanderbilt/Documents/LGN_data_042021/single_units/inverted_power_channels/good_single_units_data_4bumps_more/new_peak_alignment_anal/su_peaks_03032020_corrected/orig_peak_values/*.mat")
path <- 'C:/Users/daumail/OneDrive - Vanderbilt/Documents/LGN_data_042021/single_units/inverted_power_channels/good_single_units_data_4bumps_more/new_peak_alignment_anal/su_peaks_03032020_corrected/orig_peak_values/all_units'
bumpsfile <-file.path(path,"all_data_peaks.mat")
bumps_data <- readMat(bumpsfile)
peaks =bumps_data[["peak.vals"]]
bfs <- array(0, c(length(folderfilenames),1))
pvals <- array(0, c(length(folderfilenames),3))
tvals <- 0
df.KR <- 0
p.KR <- 0

logical_filenames = logical(length(folderfilenames))

for(i in 1:length(folderfilenames)){
  
  #filename <- filenames[i]
  #path <- system.file("C:/Users/daumail/OneDrive - Vanderbilt/Documents/LGN_data/single_units/inverted_power_channels/good_single_units_data_4bumps_more/new_peak_alignment_anal/su_peaks_03032020_corrected/orig_peak_values", package = "R.matlab")
  #pathname <- file.path("C:/Users/daumail/OneDrive - Vanderbilt/Documents/LGN_data/single_units/inverted_power_channels/good_single_units_data_4bumps_more/new_peak_alignment_anal/su_peaks_03032020_corrected/orig_peak_values",folderfilenames[i])
  #bumps_data <- readMat(data_filename_list[i])
  #bumpsfile <-file.path(path,folderfilenames[i])
  #bumps_data <- readMat(bumpsfile)  
  
  if(length(bumps_data$"peak.vals"[[i]])>0) { 
    
    if(length(filenames$layer[i])>0){ #ignore for the data analysis
      logical_filenames[i] <- TRUE #ignore for the data analysis
    } #ignore for the data analysis
    
    peak1<- t(bumps_data$"peak.vals"[[i]][1,])
    peak2<- t(bumps_data$"peak.vals"[[i]][2,])
    peak3<- t(bumps_data$"peak.vals"[[i]][3,])
    peak4<- t(bumps_data$"peak.vals"[[i]][4,])
    col_bumps_data <- t(cbind(peak1,peak2,peak3,peak4))
    
    #create labels
    labels <- c("P1","P2","P3", "P4")
    bumpnb <- matrix(rep(labels, each =length(bumps_data$"peak.vals"[[i]][1,])), ncol=4)
    t_label <- t(cbind(t(bumpnb[,1]),t(bumpnb[,2]),t(bumpnb[,3]),t(bumpnb[,4])))
    
    #create trial index
    trial_idx = t(t(rep(1:length(bumps_data$"peak.vals"[[i]][1,]), times =4)))
    letterT =  t(t(rep("T", times =4*length(bumps_data$"peak.vals"[[i]][1,]))))
    trial_idx = paste(letterT, trial_idx)
    #create table
    org_data = data.frame(trial_idx, t_label, col_bumps_data)
    
    ##fit linear model
    linearPeaks = lmer(col_bumps_data ~ t_label +(1|trial_idx), org_data)
    linsummary <- summary(linearPeaks)
    
    ##Approximation with Kenward-Rogers = helps dfs for small samples
    #coefs <- data.frame(coef(summary(linearPeaks)))
    #df.KR <- get_ddf_Lb(linearPeaks, fixef(linearPeaks))
    #p.KR <- 2 * (1 - pt(abs(coefs$t.value), df.KR))
    #pvals[i,] <- p.KR
    #m <- lmer(Reaction~factor(Days)+(1|Subject), data=lme4::sleepstudy)
    
    ##Same, with pvalue adjustment for multiple comparisons using Dunnett's correction (treatment versus control)
    #linearPeaks.emm <- emmeans(linearPeaks, "t_label")
    #contrast(linearPeaks.emm, 'trt.vs.ctrl') %>%
    #broom::tidy() %>%
    #head %>%
    #pander
    
    #Kenward-Rogers approximation + Dunnetts adjustment of pvalues
    linearPeaks.emm <- emmeans(linearPeaks, specs = trt.vs.ctrl ~ t_label)
    res <- linearPeaks.emm$contrasts %>%
      as.data.frame()
    pvals[i,] <- res$p.value
    #Compute bayes factor
    ## Bayes factor of model against null
    #org_data$trial_idx= factor(org_data$trial_idx)
    #factor(org_data$t_label)
    #bf = lmBF(col_bumps_data ~ t_label + trial_idx, data = org_data, whichRandom = "trial_idx")#, posterior=TRUE, iterations=10000)
    #bfs[i] <- bf
    
    #for(j in 1:4){ 
    #if you wanna consider the distribution normal (not recommended if small sample)
    #Vcov <- vcov(linearPeaks, useScale = FALSE) 
    #betas <- fixef(linearPeaks) 
    #se <- sqrt(diag(Vcov)) 
    #zval <- betas / se 
    #pvals[i,] <- 2 * pnorm(abs(zval), lower.tail = FALSE)
    #pvals[i,] <- 2 * (1 - pnorm(abs(coefs$t.value)))
    #}
    
    ##save output in a mat file
    #filename2 <- paste("C:/Users/maier/Documents/LGN_data/single_units/inverted_power_channels/good_single_units_data_4bumps_more/new_peak_alignment_anal/lmer_results_peaks",filename, ".mat", sep = "")
    #writeMat(filename2, table = linsummary)
    
    print(res$p.value) 
    #print(filename)
    #print(linsummary)
  } else {
    pvals[i,] <- rep(NaN, times = 3)
    print('File is empty')
  }
}
#print(tvals)
#print(df.KR)
#print(p.KR)
print(pvals)
#save Bayes factors as a .mat file
#filenamebf <- paste(path,"/", "bayesFactors.mat", sep = "")
#writeMat(filenamebf, bayesfactors = bfs)

#save pvalues as a .csv file
filename4 <- paste("C:/Users/daumail/OneDrive - Vanderbilt/Documents/LGN_data_042021/single_units/inverted_power_channels/good_single_units_data_4bumps_more/new_peak_alignment_anal/lmer_results_peaks/lmer_results_orig_03032020_corrected_dunnett", ".csv", sep = "")
write.csv(pvals, file = filename4)

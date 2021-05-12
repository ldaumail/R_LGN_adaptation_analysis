library(R.matlab)
library(lmerTest)
library(car)
library(lme4)
library(pbkrtest)
source('C:/Users/daumail/OneDrive - Vanderbilt/Documents/R/useful_functions/get_ddf_Lb.R')
source('C:/Users/daumail/OneDrive - Vanderbilt/Documents/R/useful_functions/rowMaxs.R')

#define function colMax
colMax <- function(X) apply(X, 2, max, na.rm = T)
rowMin <- function(X) apply(X, 1, min, na.rm = T)

##
## create a for loop to treat all the data at once

#file names of all 71 units in the folder
#ignore for the data analysis
path = 'C:/Users/daumail/OneDrive - Vanderbilt/Documents/LGN_data_042021/multi_units/adaptation_analysis/all_channels/'
file = file.path(path,'all_orig_bs_zscore_trials_05122021_mono_bino.mat' )

TRresps = readMat(file) #peak values obtained with BinocularAdaptationTrialSelection.m

origPks = vector(mode = "list", length = length(TRresps$peak.aligned.trials))
tempList = list()
#cellclass = array(0, length(TRresps$peak.aligned.trials)) #let's ignore the cell class for now

for(i in 1:length(TRresps$peak.aligned.trials)){
  for(b in 1:2){
    for(p in 1:4){
      if (length(TRresps$peak.aligned.trials[[i]][[1]]) == 2){
        tempList[[p]] = colMax(TRresps$peak.aligned.trials[[i]][[1]][[b]][[p]])
        origPks[[i]][[b]] = tempList
        #if(length(TRresps$peak.aligned.trials[[i]][[4]])>0){  #let's ignore the cell class for now
          #cellclass[i] = TRresps$peak.aligned.trials[[i]][[4]]
        #}
      }
    }
  }
}

#x <- lapply(1:10, function(i) i)

#folderfilenames <-list.files("C:/Users/daumail/OneDrive - Vanderbilt/Documents/LGN_data/single_units/inverted_power_channels/good_single_units_data_4bumps_more/new_peak_alignment_anal/su_peaks_03032020_corrected/orig_peak_values/", pattern=NULL, all.files=FALSE,  full.names=FALSE)#ignore
#ignore for the data analysis
#folderfilenames <- folderfilenames[1:71] #ignore for the data analysis

#file names of all 71 units with cell class label
#filenames <- read.csv("C:/Users/daumail/OneDrive - Vanderbilt/Documents/LGN_data/single_units/inverted_power_channels/good_single_units_data_4bumps_more/new_peak_alignment_anal/su_peaks_03032020_corrected/orig_peak_values/all_units/filenames_layers.csv")
library(tidyr)
#keep only labelled units file names
#row.has.na <- apply(filenames, 1, function(x){any(is.na(x))})
#filenames.filt <- filenames[!row.has.na,]


#data_filename_list <- Sys.glob("C:/Users/daumail/OneDrive - Vanderbilt/Documents/LGN_data/single_units/inverted_power_channels/good_single_units_data_4bumps_more/new_peak_alignment_anal/su_peaks_03032020_corrected/orig_peak_values/*.mat")
pvals <- array(0, c(length(origPks),4,2))
tvals <- 0
df.KR <- 0
p.KR <- 0
cnt <- 0
#logical_filenames = logical(length(folderfilenames))

for(i in 1:length(origPks)){
  if(length(origPks[[i]]) == 2){
    cnt = cnt +1
    for(b in 1:2){
      #filename <- filenames[i]
      #path <- system.file("C:/Users/daumail/OneDrive - Vanderbilt/Documents/LGN_data/single_units/inverted_power_channels/good_single_units_data_4bumps_more/new_peak_alignment_anal/su_peaks_03032020_corrected/orig_peak_values", package = "R.matlab")
      #pathname <- file.path("C:/Users/daumail/OneDrive - Vanderbilt/Documents/LGN_data/single_units/inverted_power_channels/good_single_units_data_4bumps_more/new_peak_alignment_anal/su_peaks_03032020_corrected/orig_peak_values",folderfilenames[i])
      bumps_data <- origPks[[i]][[b]]
      #if(length(bumps_data)>0) { 
      #  if(length(filenames$layer[i])>0){ #ignore for the data analysis
      #    logical_filenames[i] <- TRUE #ignore for the data analysis
      #  } #ignore for the data analysis
      
      peak1<- t(bumps_data[[1]])
      peak2<- t(bumps_data[[2]])
      peak3<- t(bumps_data[[3]])
      peak4<- t(bumps_data[[4]])
      col_bumps_data <- t(cbind(peak1,peak2,peak3,peak4))
      
      #create labels
      labels <- c("P1","P2","P3", "P4")
      bumpnb <- matrix(rep(labels, each =length(bumps_data[[1]])), ncol=4)
      t_label <- t(cbind(t(bumpnb[,1]),t(bumpnb[,2]),t(bumpnb[,3]),t(bumpnb[,4])))
      
      #create trial index
      trial_idx = t(t(rep(1:length(bumps_data[[1]]), times =4)))
      
      #create table
      org_data = data.frame(trial_idx, t_label, col_bumps_data)
      
      ##fit linear model
      linearPeaks = lmer(col_bumps_data ~ t_label +(1|trial_idx), org_data)
      linsummary <- summary(linearPeaks)
      
      coefs <- data.frame(coef(summary(linearPeaks)))
      df.KR <- get_ddf_Lb(linearPeaks, fixef(linearPeaks))
      p.KR <- 2 * (1 - pt(abs(coefs$t.value), df.KR))
      pvals[i,1:4,b] <- p.KR
      
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
      
      print(p.KR) 
      #print(filename)
      #print(linsummary)
      #} else {
      #  pvals[i,1:4,1] <- rep(NaN, times = 4)
      #  print('File is empty')
      #}
    }
  }
}
#print(tvals)
#print(df.KR)
#print(p.KR)
print(pvals)
#save pvalues as a .mat file
#filename3 <- paste("C:/Users/maier/Documents/LGN_data/single_units/inverted_power_channels/good_single_units_data_4bumps_more/new_peak_alignment_anal/lmer_results_peaks/lmer_results_02212020", ".mat", sep = "")
#writeMat(filename3, pvalues = pvals)

#save pvalues as a .csv file

filename4 <- paste(path,"adaptation_pvalues_mono_bino_05122021", ".csv", sep = "")
write.csv(pvals, file = filename4)



## Count number of adapting units for each cell class


cntsupr = array(0,c(1,2))
cntfac = array(0,c(1,2))
cntn = array(0,c(1,2))
for(b in 1:2){ # condition
  for(i in 1:length(origPks)){ #nunits
   # cell = cellclass[i]
  #  if( cell == 'M'){
   #   c <- 1
  #  }
   # if( cell == 'P'){
  #    c <- 2
   # }
  #  if( cell == 'K'){
   #   c <- 3
  #  }
  #  if( cell == '0'){
  #   c <- 4
   # }
    if(is.null(origPks[[i]]) == F){
      #compute mean pk1
      meanpk1 = mean(origPks[[i]][[b]][[1]])
      #compute mean pk4
      meanpk4 = mean(origPks[[i]][[b]][[4]])
      #get pvalue

      if( meanpk1>meanpk4 && pvals[i,4,b] < 0.05){
        cntsupr[1,b] = cntsupr[1,b]+1
      }
      if( meanpk1<meanpk4 && pvals[i,4,b] < 0.05){
        cntfac[1,b] = cntfac[1,b]+1
      }
      if( pvals[i,4,b]> 0.05){
        cntn[1,b] = cntn[1,b]+1
      }
      
    }
    
  }
}


#create data frame with 8 rows and 5 columns
results_df <- data.frame(matrix(ncol = 5, nrow = 8))

#provide column names
colnames(results_df) <- c('Cell Class', 'Condition', 'Adapt Count', 'Facilit Count', 'NS Count')

#create labels
cell_labels <- c("M","P","K", "Unknown")
class_col <- matrix(rep(cell_labels, times =2), ncol=1)
condition <- c("Monocular", "Binocular")
cond_col <- matrix(rep(condition, each =4), ncol=1)

results_df$`Cell Class`=class_col
results_df$Condition =cond_col
results_df$`Adapt Count`=c(cntsupr[,1], cntsupr[,2])
results_df$`Facilit Count`=c(cntfac[,1], cntfac[,2])
results_df$`NS Count`=c(cntn[,1], cntn[,2])

print(results_df)

#Look at units that are significantly modulated between mono and binocular conditions
#Test for significant modulation with wilkoxon ranksum test (= independant samples)

Wpvals <- array(0, c(length(origPks)))
for(i in 1:length(origPks)){
  if(length(origPks[[i]]) == 2){
    peaks1_mono <- origPks[[i]][[b]][[1]]
    peaks1_bino <- origPks[[i]][[b]][[2]]
    W= wilcox.test(peaks1_mono,peaks1_bino)
    Wpvals[i]<- W$p.value
  } else {
    Wpvals[i] <- NA
    
  }
}

##Create summary table with all binocular modulation pvalues (Wpvals) together with adaptation pvalues and mean of peak values

sumTable <- data.frame(matrix(ncol = 10,nrow = 2*length(origPks)))
colnames(sumTable) <- c('Cell Class', 'Condition', 'Pk1Pk4pvalue', 'Wpvalue', 'Pk1', 'Pk2', 'Pk3', 'Pk4', 'Pk1Pk4supress','binosup')
sumTable$`Cell Class`<- c(cellclass,cellclass)
condition <- c("Monocular", "Binocular")
sumTable$`Condition`<- matrix(rep(condition, each =length(origPks)), ncol=1)
sumTable$Pk1Pk4pvalue<- c(pvals[,4,1],pvals[,4,2])
sumTable$Wpvalue<- c(Wpvals,Wpvals) #useless to fill twice but just to fill up the table

#### THIS code segement is to replace oriPks values by normalized (z-score or simple normalizing process)  peak values for plotting comparisons
origPks <- vector(mode = "list", length = length(TRresps$peak.aligned.trials))
tempList <- list()
for(i in 1:length(TRresps$peak.aligned.trials)){
  if (length(TRresps$peak.aligned.trials[[i]][[1]]) == 2){
    for(b in 1:2){
      for(p in 1:4){
        
        tempList[[p]] <- colMax(TRresps$peak.aligned.trials[[i]][[3]][[b]][[p]])
        
      }
      
      matTempList <- matrix(unlist(tempList), ncol = 4)
      normTempMat <- (matTempList - rowMin(matTempList))/(rowMaxs(matTempList)-rowMin(matTempList))
      normTempList <- split(normTempMat, rep(1:ncol(normTempMat), each = nrow(normTempMat)))
      origPks[[i]][[b]] <- normTempList
    }
  }
}

#### ENd of  normalized segment

#fill up mean peak values column
for(i in 1:length(origPks)){
  for(b in 1:2){
    if(is.null(origPks[[i]][[b]]) == FALSE){
      meanPks = colMeans(matrix(unlist(origPks[[i]][[b]]), ncol = 4, byrow = FALSE))
      if(b == 1){
        sumTable$Pk1[i] = meanPks[1]
        sumTable$Pk2[i] = meanPks[2]
        sumTable$Pk3[i] = meanPks[3]
        sumTable$Pk4[i] = meanPks[4]
      } else if (b == 2){
        sumTable$Pk1[length(origPks)+i] = meanPks[1]
        sumTable$Pk2[length(origPks)+i] = meanPks[2]
        sumTable$Pk3[length(origPks)+i] = meanPks[3]
        sumTable$Pk4[length(origPks)+i] = meanPks[4]
      }
    }
  }
}
sumTable$Pk1Pk4supress = (sumTable$Pk1 - sumTable$Pk4) > 0
sumTable$binosup = matrix(rep(sumTable$Pk1[sumTable$Condition == 'Monocular'] - sumTable$Pk1[sumTable$Condition == 'Binocular'] > 0, times = 2), ncol=1)

filename5 <- paste(path,"summary_table_pvalues__meanpks_mono_bino_05022021", ".csv", sep = "")
write.csv(sumTable, file = filename5)


#Plot meanPks of units with significant binocular modulation vs peaks that don't (all units showing adaptation).
#All cells
row.has.na <- apply(sumTable, 1, function(x){any(is.na(x))})
sumTable.filt <- sumTable[!row.has.na,] #remove rows with NA
row_sub = apply(sumTable.filt, 1, function(row) all(row !=0 ))
sumTable.filt <- sumTable.filt[row_sub,] #remove rows with 0 (no cell class label)


library("ggpubr")
library(ggplot2)
library(reshape2)
library(gridGraphics)
library(cowplot)


plot_dat <- melt(sumTable.filt, id.var=c('Condition','Cell Class','Pk1Pk4pvalue','Wpvalue','Pk1Pk4supress', 'binosup'), variable.name = 'Peak')


#############################################                                  1                                ################################################
#Plot units with binocular modulation (any type of modulation)
#make appropriate dataframe
mod_df = plot_dat[plot_dat$Wpvalue<0.05,]

#count number of units that are included in plot
count = paste('Units Count:', length(plot_dat$Peak[plot_dat$Wpvalue<=0.05 & plot_dat$Condition == 'Monocular' & plot_dat$Peak == 'Pk1']))
#means = aggregate(plot_dat[plot_dat$Wpvalue<=0.05,], list(plot_dat$Peak[plot_dat$Wpvalue<=0.05],plot_dat$Condition[plot_dat$Wpvalue<=0.05]), mean)


plot
# plot
ggplot(mod_df,aes(y = value, x = Peak,fill=Condition)) + 
  
  # box plots and jitter points, with modified x value
  geom_violin(width=0.6, trim=F, color = NA)+
  stat_summary(fun.y = mean, geom = "point", shape = 23, size = 2,width=0.6)+
  #geom_boxplot(aes(middle = median(value)),width=0.6) + #
  #geom_boxplot(aes(x=Peak, y=value,fill=Condition,middle = median(value)),width=0.6) +
  #geom_bar(aes(x=Peak,fill=Condition),width=0.6)
  
  theme_bw() +
  theme(
    plot.title   = element_text(color = "black", size = 15, face = "bold", hjust=0.5),
    axis.title.x = element_text(color = "black", size = 13),
    axis.title.y = element_text(color = "black", size = 13),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text    = element_text(size  = 9),
    axis.line = element_line(colour = "black", size = 1) 
  )+
  
  labs(
    title = paste("Comparison of modulated units responses \n in monocular and binocular conditions \n",count),
    x = "Peak",
    #y = "Spike rate (spikes/sec)" )
    y = "Spike rate (normalized)" )
ggsave(filename= "C:/Users/daumail/OneDrive - Vanderbilt/Documents/LGN_data_042021/single_units/binocular_adaptation/plots/comparison_mono_bino_signif_modul_mean_normalized.svg")
ggsave(filename= "C:/Users/daumail/OneDrive - Vanderbilt/Documents/LGN_data_042021/single_units/binocular_adaptation/plots/comparison_mono_bino_signif_modul_mean_normalized.png")

#grid.newpage()
#grid.draw(as_grob(plot))+
#ggdraw(count) + draw_label("bottom left at (0, 0)", x = 0, y = 0, hjust = 0, vjust = 0)

##########################################                                        2                                        ######################################
#Same plot of units with binocular modulation but with significantly adaptating units only
mod_df = plot_dat[plot_dat$Wpvalue<0.05,]
mono_adapt <- mod_df$Condition == 'Monocular' & mod_df$Pk1Pk4pvalue <= 0.05 & mod_df$Peak == 'Pk1'
bino_adapt <- mod_df$Condition == 'Binocular' & mod_df$Pk1Pk4pvalue <= 0.05 & mod_df$Peak == 'Pk1'
mono_bino_adapt <- rep(mono_adapt[mod_df$Condition == 'Monocular' & mod_df$Peak == 'Pk1'], times =8) | rep(bino_adapt[mod_df$Condition == 'Binocular' & mod_df$Peak == 'Pk1'], times =8)
mod_adapt_df = mod_df[mono_bino_adapt,]
#count number of units that are included in plot
count = paste('Units Count:', length(mod_adapt_df$Peak[mod_adapt_df$Peak == 'Pk1' & mod_adapt_df$Condition == 'Monocular'] ))

plot
# plot
ggplot(mod_adapt_df,aes(y = value, x = Peak,fill=Condition)) + 
  
  # box plots and jitter points, with modified x value
  geom_violin(aes(middle = mean(value)),width=0.6, trim=F, color = NA)+
  stat_summary(fun.y = mean, geom = "point", shape = 23, size = 2,width=0.6)+
  #geom_boxplot(aes(middle = mean(value)),width=0.6) +
  #geom_boxplot(aes(x=Peak, y=value,fill=Condition),width=0.6) +
  #geom_text(data = means, aes(label = round(x, 1), y = x + 1), size = 3) + #adds average labels
  #geom_text(data = medians, aes(label = round(x, 1), y = x - 0.5), size = 3) + #adds median labels
  
  theme_bw() +
  theme(
    plot.title   = element_text(color = "black", size = 15, face = "bold", hjust=0.5),
    axis.title.x = element_text(color = "black", size = 13),
    axis.title.y = element_text(color = "black", size = 13),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text    = element_text(size  = 9),
    axis.line = element_line(colour = "black", size = 1) 
  )+
  
  labs(
    title = paste("Comparison of modulated adapting units \n responses in monocular and binocular \n conditions \n",count),
    x = "Peak",
    #y = "Spike rate (spikes/sec)" )
    y = "Spike rate (normalized)" )
ggsave(filename= "C:/Users/daumail/OneDrive - Vanderbilt/Documents/LGN_data_042021/single_units/binocular_adaptation/plots/comparison_mono_bino_signif_modul_adapt_mean_normalized.svg")
ggsave(filename= "C:/Users/daumail/OneDrive - Vanderbilt/Documents/LGN_data_042021/single_units/binocular_adaptation/plots/comparison_mono_bino_signif_modul_adapt_mean_normalized.png")

#################################                               3                                 ####################################################
#Make same plot units with significant binocular modulation and SUPPRESSED units in monocular condition only 

#make approriate plot dataframe
mod_df <- plot_dat[plot_dat$Wpvalue<0.05,]
mono_supr <- rep(mod_df$Pk1Pk4supress[mod_df$Condition == 'Monocular' & mod_df$Peak == 'Pk1'], times = 8) #supression in monocular condition
sig_mono <- mod_df$Condition == 'Monocular' & mod_df$Peak == 'Pk1' & mod_df$Pk1Pk4pvalue <= 0.05 #significant adaptation in monocular condition
sig_mono_supr <- mono_supr & rep(sig_mono[mod_df$Condition == 'Monocular' & mod_df$Peak == 'Pk1'], times =8) #combine both
mod_adaptSupr_df = mod_df[sig_mono_supr,]  #only keep corresponding units
#count = paste('Units Count:', length(plot_dat$Peak[plot_dat$Wpvalue<=0.05 & plot_dat$Condition == 'Monocular' & plot_dat$Peak == 'Pk1' & plot_dat$Pk1Pk4pvalue<=0.05 & (plot_dat$value[plot_dat$Peak == 'Pk1']>plot_dat$value[plot_dat$Peak == 'Pk4'])]))
count = paste('Units Count:', length(mod_adaptSupr_df$Peak[mod_adaptSupr_df$Condition == 'Monocular' & mod_adaptSupr_df$Peak == 'Pk1']))

#mod_adapt_df$id <- rep(1:10,times = 4)
#dcast(data = mod_adapt_df,formula = id~Peak,fun.aggregate = sum,value.var = "value")
#mod_adaptSupr_df = mod_adapt_df[mod_adapt_df$Pk1Pk4supress == TRUE,]

plot
# plot
ggplot(mod_adaptSupr_df,aes(y = value, x = Peak,fill=Condition)) + 
  
  # box plots and jitter points, with modified x value
  geom_violin(aes(middle = mean(value)),width=0.6, trim=F, color = NA)+
  stat_summary(fun.y = mean, geom = "point", shape = 23, size = 2,width=0.6)+
  #geom_boxplot(aes(middle = mean(value)),width=0.6) +
  #geom_boxplot(aes(x=Peak, y=value,fill=Condition),width=0.6) +
  #geom_text(data = means, aes(label = round(x, 1), y = x + 1), size = 3) + #adds average labels
  #geom_text(data = medians, aes(label = round(x, 1), y = x - 0.5), size = 3) + #adds median labels
  
  theme_bw() +
  theme(
    plot.title   = element_text(color = "black", size = 15, face = "bold", hjust=0.5),
    axis.title.x = element_text(color = "black", size = 13),
    axis.title.y = element_text(color = "black", size = 13),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text    = element_text(size  = 9),
    axis.line = element_line(colour = "black", size = 1) 
  )+
  
  labs(
    title = paste("Comparison of bino modulated sig suppressed \n units responses in monocular and binocular \n conditions \n",count),
    x = "Peak",
    #y = "Spike rate (spikes/sec)" )
    y = "Spike rate (normalized)" )
ggsave(filename= "C:/Users/daumail/OneDrive - Vanderbilt/Documents/LGN_data_042021/single_units/binocular_adaptation/plots/comparison_mono_bino_signif_modul_suppressed_mean_normalized.svg")
ggsave(filename= "C:/Users/daumail/OneDrive - Vanderbilt/Documents/LGN_data_042021/single_units/binocular_adaptation/plots/comparison_mono_bino_signif_modul_suppressed_mean_normalized.png")



################################                                              4                                         ###################################
#Make same plot units with significant binocular suppression and SUPPRESSED adaptation in monocular condition  
#count = paste('Units Count:', length(plot_dat$Peak[plot_dat$Wpvalue<=0.05 & plot_dat$Condition == 'Monocular' & plot_dat$Peak == 'Pk1' & plot_dat$Pk1Pk4pvalue<=0.05 & (plot_dat$value[plot_dat$Peak == 'Pk1'  & plot_dat$Condition == 'Monocular']>plot_dat$value[plot_dat$Peak == 'Pk4' & plot_dat$Condition == 'Monocular']) & (plot_dat$value[plot_dat$Peak == 'Pk1' & plot_dat$Condition == 'Monocular']>plot_dat$value[plot_dat$Peak == 'Pk1' & plot_dat$Condition == 'Binocular'])]))
bino_supr_df <- plot_dat[plot_dat$Wpvalue<0.05 & plot_dat$binosup == TRUE,]
mono_supr <- rep(bino_supr_df$Pk1Pk4supress[bino_supr_df$Condition == 'Monocular' & bino_supr_df$Peak == 'Pk1'], times = 8) #supression in monocular condition
sig_mono <- bino_supr_df$Condition == 'Monocular' & bino_supr_df$Peak == 'Pk1' & bino_supr_df$Pk1Pk4pvalue <= 0.05 #significant adaptation in monocular condition
sig_mono_supr <- mono_supr & rep(sig_mono[bino_supr_df$Condition == 'Monocular' & bino_supr_df$Peak == 'Pk1'], times =8) #combine both
binosupr_adaptSupr_df = bino_supr_df[sig_mono_supr,]  #only keep corresponding units
#count = paste('Units Count:', length(plot_dat$Peak[plot_dat$Wpvalue<=0.05 & plot_dat$Condition == 'Monocular' & plot_dat$Peak == 'Pk1' & plot_dat$Pk1Pk4pvalue<=0.05 & (plot_dat$value[plot_dat$Peak == 'Pk1']>plot_dat$value[plot_dat$Peak == 'Pk4'])]))
count = paste('Units Count:', length(binosupr_adaptSupr_df$Peak[binosupr_adaptSupr_df$Condition == 'Monocular' & binosupr_adaptSupr_df$Peak == 'Pk1']))

#bino_supr_df = plot_dat[plot_dat$Wpvalue<=0.05 & (plot_dat$value[plot_dat$Condition == 'Monocular' & plot_dat$Peak == 'Pk1']>plot_dat$value[plot_dat$Condition == 'Binocular' & plot_dat$Peak == 'Pk1']),]
#bss_supr_df = bino_supr_df[bino_supr_df$Pk1Pk4pvalue<=0.05 & (bino_supr_df$value[bino_supr_df$Peak[bino_supr_df$Condition == 'Monocular'] == 'Pk1' ]>bino_supr_df$value[bino_supr_df$Peak[bino_supr_df$Condition == 'Monocular'] == 'Pk4']),]

plot
# plot
ggplot(binosupr_adaptSupr_df,aes(y = value, x = Peak,fill=Condition)) + 
  
  # box plots and jitter points, with modified x value
  geom_violin(width=0.6, trim=F, color = NA)+
  stat_summary(fun.y = mean, geom = "point", shape = 23, size = 2)+
  #geom_boxplot(aes(middle = mean(value)),width=0.6) +
  #geom_boxplot(aes(x=Peak, y=value,fill=Condition),width=0.6) +
  #geom_text(data = means, aes(label = round(x, 1), y = x + 1), size = 3) + #adds average labels
  #geom_text(data = medians, aes(label = round(x, 1), y = x - 0.5), size = 3) + #adds median labels
  
  theme_bw() +
  theme(
    plot.title   = element_text(color = "black", size = 15, face = "bold", hjust=0.5),
    axis.title.x = element_text(color = "black", size = 13),
    axis.title.y = element_text(color = "black", size = 13),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text    = element_text(size  = 9),
    axis.line = element_line(colour = "black", size = 1) 
  )+
  
  labs(
    title = paste("Comparison of binocularly supressed, suppressed \n units responses in monocular and binocular \n conditions \n",count),
    x = "Peak",
    #y = "Spike rate (spikes/sec)" )
    y = "Spike rate (normalized)" )
ggsave(filename= "C:/Users/daumail/OneDrive - Vanderbilt/Documents/LGN_data_042021/single_units/binocular_adaptation/plots/comparison_mono_bino_signif_binosupr_suppressed_mean_normalized.svg")
ggsave(filename= "C:/Users/daumail/OneDrive - Vanderbilt/Documents/LGN_data_042021/single_units/binocular_adaptation/plots/comparison_mono_bino_signif_binosupr_suppressed_mean_normalized.png")




#           Less Important Plots
#Make same plot with significantly modulated and FACILITATED units in monocular condition only 
count = paste('Units Count:', length(plot_dat$Peak[plot_dat$Wpvalue<=0.05 & plot_dat$Condition == 'Monocular' & plot_dat$Peak == 'Pk1' & plot_dat$Pk1Pk4pvalue<=0.05 & (plot_dat$value[plot_dat$Peak == 'Pk1']<plot_dat$value[plot_dat$Peak == 'Pk4'])]))

plot
# plot
ggplot(plot_dat[plot_dat$Wpvalue<=0.05 & plot_dat$Pk1Pk4pvalue<=0.05 & (plot_dat$value[plot_dat$Condition == 'Monocular' & plot_dat$Peak == 'Pk1']<plot_dat$value[plot_dat$Condition == 'Monocular' & plot_dat$Peak == 'Pk4']),],aes(y = value, x = Peak,fill=Condition)) + 
  
  # box plots and jitter points, with modified x value
  geom_violin(aes(middle = mean(value)),width=0.6, trim=F, color = NA)+
  stat_summary(fun.y = mean, geom = "point", shape = 23, size = 2,width=0.6)+
  #geom_boxplot(aes(middle = mean(value)),width=0.6) +
  #geom_boxplot(aes(x=Peak, y=value,fill=Condition),width=0.6) +
  #geom_text(data = means, aes(label = round(x, 1), y = x + 1), size = 3) + #adds average labels
  #geom_text(data = medians, aes(label = round(x, 1), y = x - 0.5), size = 3) + #adds median labels
  
  theme_bw() +
  theme(
    plot.title   = element_text(color = "black", size = 15, face = "bold", hjust=0.5),
    axis.title.x = element_text(color = "black", size = 13),
    axis.title.y = element_text(color = "black", size = 13),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text    = element_text(size  = 9),
    axis.line = element_line(colour = "black", size = 1) 
  )+
  
  labs(
    title = paste("Comparison of modulated facilitated (mono) \n units responses in monocular and binocular \n conditions \n",count),
    x = "Peak",
    y = "Spike rate (spikes/sec)" )
ggsave(filename= "C:/Users/daumail/OneDrive - Vanderbilt/Documents/LGN_data_042021/single_units/binocular_adaptation/plots/comparison_mono_bino_signif_modul_facilit_mono_mean.svg")
ggsave(filename= "C:/Users/daumail/OneDrive - Vanderbilt/Documents/LGN_data_042021/single_units/binocular_adaptation/plots/comparison_mono_bino_signif_modul_facilit_mono_mean.png")

##
#Make same plot with significantly modulated and SUPPRESSED units in BInocular condition only 
count = paste('Units Count:', length(plot_dat$Peak[plot_dat$Wpvalue<=0.05 & plot_dat$Condition == 'Binocular' & plot_dat$Peak == 'Pk1' & plot_dat$Pk1Pk4pvalue<=0.05 & (plot_dat$value[plot_dat$Peak == 'Pk1']>plot_dat$value[plot_dat$Peak == 'Pk4'])]))

plot
# plot
ggplot(plot_dat[plot_dat$Wpvalue<=0.05 & plot_dat$Pk1Pk4pvalue<=0.05 & (plot_dat$value[plot_dat$Condition == 'Binocular' & plot_dat$Peak == 'Pk1']>plot_dat$value[plot_dat$Condition == 'Binocular' & plot_dat$Peak == 'Pk4']),],aes(y = value, x = Peak,fill=Condition)) + 
  
  # box plots and jitter points, with modified x value
  geom_violin(aes(middle = mean(value)),width=0.6, trim=F, color = NA)+
  stat_summary(fun.y = mean, geom = "point", shape = 23, size = 2,width=0.6)+
  #geom_boxplot(aes(middle = mean(value)),width=0.6) +
  #geom_boxplot(aes(x=Peak, y=value,fill=Condition),width=0.6) +
  #geom_text(data = means, aes(label = round(x, 1), y = x + 1), size = 3) + #adds average labels
  #geom_text(data = medians, aes(label = round(x, 1), y = x - 0.5), size = 3) + #adds median labels
  
  theme_bw() +
  theme(
    plot.title   = element_text(color = "black", size = 15, face = "bold", hjust=0.5),
    axis.title.x = element_text(color = "black", size = 13),
    axis.title.y = element_text(color = "black", size = 13),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text    = element_text(size  = 9),
    axis.line = element_line(colour = "black", size = 1) 
  )+
  
  labs(
    title = paste("Comparison of modulated suppressed (bino) \n units responses in monocular and binocular \n conditions \n",count),
    x = "Peak",
    y = "Spike rate (spikes/sec)" )
ggsave(filename= "C:/Users/daumail/OneDrive - Vanderbilt/Documents/LGN_data_042021/single_units/binocular_adaptation/plots/comparison_mono_bino_signif_modul_suppressed_bino_mean.svg")
ggsave(filename= "C:/Users/daumail/OneDrive - Vanderbilt/Documents/LGN_data_042021/single_units/binocular_adaptation/plots/comparison_mono_bino_signif_modul_suppressed_bino_mean.png")


cnt = 0
invcnt = 0
intcnt = 0

for(i in 1:length(origPks)){
  if(is.na(Wpvals[i])==F){
    if(Wpvals[i] <0.05){ #&& (pvals[i,4,1]<0.05 && pvals[i,4,2]<0.05)){
      cnt = cnt+1
    }
    if(Wpvals[i] <0.05 && ((pvals[i,4,1]<0.05 && pvals[i,4,2]>0.05) || (pvals[i,4,1]>0.05 && pvals[i,4,2]<0.05))){
      invcnt = invcnt+1
    }
    if(Wpvals[i] <0.05 && (pvals[i,4,1]<0.05 || pvals[i,4,2]>0.05)){
      intcnt = intcnt +1
    }
  }
}

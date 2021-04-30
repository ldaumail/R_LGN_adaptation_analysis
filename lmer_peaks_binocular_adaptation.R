library(R.matlab)
library(lmerTest)
library(car)
library(lme4)
library(pbkrtest)
source('C:/Users/daumail/OneDrive - Vanderbilt/Documents/R/useful_functions/get_ddf_Lb.R')
source('C:/Users/daumail/OneDrive - Vanderbilt/Documents/R/useful_functions/rowMaxs.R')

#define function colMax
colMax <- function(X) apply(X, 2, max, na.rm = T)

##
## create a for loop to treat all the data at once

#file names of all 71 units in the folder
#ignore for the data analysis
path = 'C:/Users/daumail/OneDrive - Vanderbilt/Documents/LGN_data_042021/single_units/binocular_adaptation/all_units/'
file = file.path(path,'all_orig_bs_zscore_trials.mat' )
TRresps = readMat(file) #peak values obtained with BinocularAdaptationTrialSelection.m

origPks = vector(mode = "list", length = length(TRresps$peak.aligned.trials))
tempList = list()
cellclass = array(0, length(TRresps$peak.aligned.trials))

for(i in 1:length(TRresps$peak.aligned.trials)){
  for(b in 1:2){
    for(p in 1:4){
     tempList[[p]] = colMax(TRresps$peak.aligned.trials[[i]][[1]][[b]][[p]])
     origPks[[i]][[b]] = tempList
     if(length(TRresps$peak.aligned.trials[[i]][[4]])>0){
     cellclass[i] = TRresps$peak.aligned.trials[[i]][[4]]
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

#logical_filenames = logical(length(folderfilenames))

for(i in 1:length(origPks)){
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
#print(tvals)
#print(df.KR)
#print(p.KR)
print(pvals)
#save pvalues as a .mat file
#filename3 <- paste("C:/Users/maier/Documents/LGN_data/single_units/inverted_power_channels/good_single_units_data_4bumps_more/new_peak_alignment_anal/lmer_results_peaks/lmer_results_02212020", ".mat", sep = "")
#writeMat(filename3, pvalues = pvals)

#save pvalues as a .csv file

filename4 <- paste(path,"all_pvalues_mono_bino", ".csv", sep = "")
write.csv(pvals, file = filename4)



## Count number of adapting units for each cell class

#create data frame with 0 rows and 3 columns
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
cntsupr = array(0,c(4,2))
cntfac = array(0,c(4,2))
cntn = array(0,c(4,2))
for(b in 1:2){ # condition
  for(i in 1:length(origPks)){ #nunits
    cell = cellclass[i]
    if( cell == 'M'){
      c <- 1
    }
    if( cell == 'P'){
      c <- 2
    }
    if( cell == 'K'){
      c <- 3
    }
    if( cell == '0'){
      c <- 4
    }
    #compute mean pk1
    meanpk1 = mean(origPks[[i]][[b]][[1]])
    #compute mean pk4
    meanpk4 = mean(origPks[[i]][[b]][[4]])
    #get pvalue
    
    
    if( meanpk1>meanpk4 && pvals[i,4,b] < 0.09){
      cntsupr[c,b] = cntsupr[c,b]+1
    }
    if( meanpk1<meanpk4 && pvals[i,4,b] < 0.09){
      cntfac[c,b] = cntfac[c,b]+1
    }
    if( pvals[i,4,b]> 0.09){
      cntn[c,b] = cntn[c,b]+1
    }
    
      
    
  }
}













#Get filenames of non empty units (e.g. the remaining processed data== clean data)
selectedfilenames <- folderfilenames[logical_filenames] #this line doesn't work as the order of the files in the folder 
#(190205..uclust6 = is 61st in the folder, and 190205...uclust61= is 62nd in the folder) are different from the order
#in the original list from the file containing the whole dataset (filenames_layer.csv lists the names in this order with 190205...uclust61 = 61st
#and 190205...uclust6= 62nd). 
#Therefore the logical_filenames vector picks the filenames in the order of the original dataset, not in the folder filenames order
#ending up picking the wrong names. A fix should be done to adjust this error, lets see when I get some time. The order is right for the analysis 
#of the data.
#As a conclusion: a difference should be noted in the way Sys.glob lists the files (order of the files in alphabetical order : right here, 
#and is the one determining the order of the analysis,
#and allows to give the right pvalue order), in comparison with list.files ( order of files in the folder : wrong here, but not impacting 
#the order of the pvalues results, only used
#to try to select the filenames, but didn't work out in this situation where the order of the filename maters)
#As a lesson to take from this: the use of logicals should be done with certainty that the data elements are picked in the right order (either alphabetical or known order such as in the folder)
selectedfilenames<- gsub(".mat","",selectedfilenames) #ignore

cellclassfn <- gsub("gmat","", filenames.filt$filename)#ignore

#get filenames of cell class labelled units

finalfilename <- list() #ignore
for( i in 1:length(selectedfilenames)){ #ignore
  for(j in 1:length(filenames.filt$filename)) #ignore
    if(selectedfilenames[i]==cellclassfn[j]){ #ignore
      finalfilename[i] <- selectedfilenames[i] #ignore
    } #ignore
} #ignore
finalfilename <-finalfilename[lengths(finalfilename) > 0] #ignore
filename5 <- paste("C:/Users/daumail/Documents/LGN_data/single_units/inverted_power_channels/good_single_units_data_4bumps_more/new_peak_alignment_anal/su_peaks_03032020_corrected/orig_peak_values/all_units/selected_units_filenames", ".csv", sep = "") #ignore
write.csv(finalfilename, file = filename5) #ignore


library(R.matlab)
library("ggpubr")
library("colorspace")

path = 'C:/Users/daumail/OneDrive - Vanderbilt/Documents/LGN_data_042021/single_units/inverted_power_channels/good_single_units_data_4bumps_more/new_peak_alignment_anal/su_peaks_03032020_corrected/orig_peak_values/all_units'
file = file.path(path,'all_data_peaks.mat')
bumps_data <- readMat(file)
peaks =bumps_data[["peak.vals"]]
#file_names contains both the filenames and the layer
file2 = file.path(path,'filenames_layers.csv')
file_names_data <- read.csv(file2, header = TRUE,sep = ',', na.strings="")

#pvalues of adaptation of peak1 vs peak4
pvalues_data <- read.csv("lmer_results_orig_03032020_corrected.csv",header = TRUE,sep = ',', na.strings="")

effectsize<-  array(NaN, c(length(file_names_data[[1]])))
peak1<- array(NaN, c(length(file_names_data[[1]])))
peak4<- array(NaN, c(length(file_names_data[[1]])))

for(i in 1:length(file_names_data[[1]])){
  if(length(peaks[[i]])>0){
    unit_peaks1 = peaks[[i]][1,]
    unit_peaks4 =peaks[[i]][4,]
    
    mean_pk1 = mean(unit_peaks1)
    mean_pk4 = mean(unit_peaks4)
    peak1[i]=mean_pk1
    peak4[i]=mean_pk4
    
    #calculate cohen's d (cd)
    sd_diff = sd(unit_peaks4-unit_peaks1)
    cd = (mean_pk4 - mean_pk1)/sd_diff
    effectsize[i]=cd

  }
}
channel_idx = t(t(rep(1:71, times =1)))

#data frame

org_data = data.frame("Unit Index"=channel_idx,"P-value"=pvalues_data[,5], "DeltaPeak"=peak4-peak1, peak1, peak4,effectsize, file_names_data)


library(tidyr)
#no_na_data <- drop_na(org_data)

row.has.na <- apply(org_data, 1, function(x){any(is.na(x))})
org_data.filt <- org_data[!row.has.na,]


mpeak1 = org_data.filt$peak1[org_data.filt$layer=="M"]
mpeak4 = org_data.filt$peak4[org_data.filt$layer=="M"]

ppeak1 = org_data.filt$peak1[org_data.filt$layer=="P"]
ppeak4 = org_data.filt$peak4[org_data.filt$layer=="P"]

kpeak1 = org_data.filt$peak1[org_data.filt$layer=="K"]
kpeak4 = org_data.filt$peak4[org_data.filt$layer=="K"]



# dependencies
library(ggplot2)
library(tibble)
source('./bootstrap_functions/theme_gar.txt')
source('./bootstrap_functions/fun.txt')
library(beepr)
library(WRS)
#library(WRS2)

##First option: using the WRS package
#yuend1<-WRS::yuend(mpeak1,mpeak4,tr=0,alpha=.05)

#here no trimming on M cell class as skewness is <.5
#mbootstraptrim <- WRS::ydbt(mpeak4,mpeak1,tr=.15,alpha=0.05,nboot=1000,side=FALSE,plotit=TRUE,op=1)
mbootstrap <- WRS::ydbt(mpeak4,mpeak1,tr=0,alpha=0.05,nboot=5000,side=FALSE,plotit=TRUE,op=1)
#trimci(mpeak1,tr=.2,alpha=.05,null.value=0)

#here trimming on P cell class as skewness is >0.5
pbootstraptrim <- WRS::ydbt(ppeak4,ppeak1,tr=.15,alpha=0.05,nboot=5000,side=FALSE,plotit=TRUE,op=1)

#pbootstrap <- WRS::ydbt(ppeak4,ppeak1,tr=0,alpha=0.05,nboot=5000,side=FALSE,plotit=TRUE,op=1)
#not enough data in K sample
#kbootstraptrim <- WRS::ydbt(kpeak1,kpeak4,tr=0,alpha=0.05,nboot=1000,side=FALSE,plotit=TRUE,op=1)


###other way to compute the bootstrap 95%CI (manually)
##The following method works for a one-sample test, but need to use ydbt function to perform the correct bootstrap on dependent sample data.
#sample
#samp = ppeak4 - ppeak1
#samp.m <- mean(samp, trim = 0.15) # sample mean
#samp.sem <- trimse(samp, tr = 0.15) # sample estimation of the standard error of the mean
# centred data
#mcsamp <-  samp - mean(samp)
#n = length(samp)
#nboot <- 5000
#alpha <- 0.05
#tval <- function(x,nv,tr=0.15){
#  se <- sqrt(winvar(x,tr)) / ((1-2*tr) * sqrt(length(x)))
#  tval <- (mean(x,tr)-nv) / se
#  tval
#}

#bootT <- apply(matrix(sample(mcsamp, nboot*n, replace = TRUE), nrow = nboot), 1, tval, nv = 0)
#bootT.q <- quantile(bootT, probs = c(alpha/2, 1-alpha/2), type = 6)

# Asymmetric CI
#CI.lo <- round(samp.m - bootT.q[2] * samp.sem, digits = 2)
#CI.up <- round(samp.m - bootT.q[1] * samp.sem, digits = 2) 
#print(paste0("Asymmetric CI = [",CI.lo,", ",CI.up,"]"))


#Plot Peak1-P4 difference histograms along with bootstrap-t 95% CIs for each cell class

#pvalues for each cell class
mpval = org_data.filt$P.value[org_data.filt$layer=="M"]
ppval = org_data.filt$P.value[org_data.filt$layer=="P"]
kpval = org_data.filt$P.value[org_data.filt$layer=="K"]

#delta for each cell class
mdelta = org_data.filt$DeltaPeak[org_data.filt$layer=="M"]
pdelta = org_data.filt$DeltaPeak[org_data.filt$layer=="P"]
kdelta = org_data.filt$DeltaPeak[org_data.filt$layer=="K"]


#means and medians for each cell class
mdelta.m <- mean(mdelta, trim = 0) # sample mean
mdelta.med <- median(mdelta,trim =0)
pdelta.m <- mean(pdelta, trim = 0.15)
pdelta.med <- median(pdelta,trim =0.15)
kdelta.m <- mean(kdelta, trim = 0)
kdelta.med <- mean(kdelta, trim = 0)
##M cells

#t.test(mpeak4, mpeak1, paired = TRUE, alternative = "two.sided")
#qt(0.975,17)


p_idx =org_data.filt$Unit.Index[org_data.filt$layer=="M"]

m_data = data.frame(p_idx,mdelta,mpval)

#only get significant decreasing
m_data$sigdata[m_data$mpval<0.05 &mdelta<0] <- "Adapting Decrease"
m_data$sigdata[m_data$mpval>0.05] <- "Non Adapting"
m_data$sigdata[m_data$mpval<0.05 &mdelta>0] <- "Adapting Increase"
#m_data[ which(m_data$mpval<0.05 & mcohen<0), ]
library(reshape2)
plot_dat <- melt(m_data, id.var=c('p_idx','sigdata','mpval'))

ggplot(plot_dat, aes(y=value))+
  geom_density(aes(y=value))+
  geom_histogram(aes(fill =sigdata),bins=20)+
  scale_fill_manual(values=c("#93C5DE","#93C5DE","#B5B5B5"))+
  scale_x_discrete(limits=c(1,2,3,4))+
  geom_hline(aes(yintercept= mbootstrap$ci[1], linetype = "Bootstrap-T 95%CI"), color = "red")+
  geom_hline(aes(yintercept=mbootstrap$ci[2],linetype="Bootstrap-T 95%CI"), color = "red")+
  geom_hline(aes(yintercept=mdelta.m, linetype="Sample Mean"), color = "black")+
  geom_hline(aes(yintercept=mbootstrap$dif, linetype="Bootstrap Mean"), color = "blue")+
 
  
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
    title = "M cell class delta P1 - P4",
    x = "Unit Count",
    y = "Spike rate difference (spikes/sec)" )
ggsave(filename= "/Users/loicdaumail/Documents/Research/adaptation_LGN_project/R/plots/m_bootstrapt_ydbt_trim0_asym95CI_deltaspikerate_histogram.png")


## M density plot

p_idx =org_data.filt$Unit.Index[org_data.filt$layer=="M"]
m_data = data.frame(p_idx,mdelta)
library(reshape2)
plot_dat <- melt(m_data, id.var=c('p_idx'))

ggplot(plot_dat, aes(y=value))+
  #geom_histogram(aes(x=..density..), alpha =.5)+
  geom_density(aes(y=value),color ="#93C5DE", fill ="#93C5DE", alpha =.3)+
  #geom_curve()
  geom_hline(aes(yintercept= mbootstrap$ci[1], linetype = "Bootstrap-T 95%CI"), color = "red")+
  geom_hline(aes(yintercept=mbootstrap$ci[2],linetype="Bootstrap-T 95%CI"), color = "red")+
  geom_hline(aes(yintercept=mdelta.m, linetype="Sample Mean"), color = "black")+
  geom_hline(aes(yintercept=mbootstrap$dif, linetype="Bootstrap Mean"), color = "#93C5DE")+
  geom_hline(aes(yintercept=mdelta.med, linetype="Sample Median"), color = "black")+
  scale_y_continuous(limits = c(-80, 35),breaks = seq(-80,35, by=10)) +
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
    title = "M cell class delta P1 - P4",
    x = "Probability",
    y = "Spike rate difference (spikes/sec)" )
ggsave(filename= "/Users/loicdaumail/Documents/Research/Maier_lab/adaptation_LGN_project/R/plots/m_bootstrapt_ydbt_trim0_asym95CI_deltaspikerate_density.svg")

##P cells

#t.test(mpeak4, mpeak1, paired = TRUE, alternative = "two.sided")
#qt(0.975,17)


p_idx =org_data.filt$Unit.Index[org_data.filt$layer=="P"]

p_data = data.frame(p_idx,pdelta,ppval)

#only get significant decreasing
p_data$sigdata[p_data$ppval<0.05 &pdelta<0] <- "Adapting Decrease"
p_data$sigdata[p_data$ppval>0.05] <- "Non Adapting"
p_data$sigdata[p_data$ppval<0.05 &pdelta>0] <- "Adapting Increase"
#m_data[ which(m_data$mpval<0.05 & mcohen<0), ]
library(reshape2)
plot_dat <- melt(p_data, id.var=c('p_idx','sigdata','ppval'))

ggplot(plot_dat, aes(y=value))+
  geom_histogram(aes(fill=sigdata),bins=20)+
  scale_fill_manual(values=c("#F6999A","#F6999A","#B5B5B5"))+
  scale_x_discrete(limits=c(1,2,3,4))+
  geom_hline(aes(yintercept=pbootstraptrim$ci[1], linetype = "Bootstrap-T 95%CI"), color = "red")+
  geom_hline(aes(yintercept=pbootstraptrim$ci[2],linetype="Bootstrap-T 95%CI"), color = "red")+
  geom_hline(aes(yintercept=pdelta.m, linetype="Sample Mean"), color = "black")+
  
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
    title = "P cell class delta P1 - P4",
    x = "Unit Count",
    y = "Spike rate difference (spikes/sec)" )
ggsave(filename= "/Users/loicdaumail/Documents/Research/adaptation_LGN_project/R/plots/p_bootstrapt_ydbt_asym95CI_deltaspikerate_histogram.png")

##P density Plot
p_idx =org_data.filt$Unit.Index[org_data.filt$layer=="P"]
p_data = data.frame(p_idx,pdelta)
library(reshape2)
plot_dat <- melt(p_data, id.var=c('p_idx'))

ggplot(plot_dat, aes(y=value))+
  #geom_histogram(aes(x=..density..), alpha =.5)+
  geom_density(aes(y=value),color ="#F6999A", fill ="#F6999A", alpha =.3)+
  #geom_curve()
  geom_hline(aes(yintercept= pbootstraptrim$ci[1], linetype = "Bootstrap-T 95%CI"), color = "#B5B5B5")+
  geom_hline(aes(yintercept=pbootstraptrim$ci[2],linetype="Bootstrap-T 95%CI"), color = "#B5B5B5")+
  geom_hline(aes(yintercept=pdelta.m, linetype="Sample Mean"), color = "black")+
  geom_hline(aes(yintercept=pbootstraptrim$dif, linetype="Bootstrap Mean"), color = "#F6999A")+
  geom_hline(aes(yintercept=pdelta.med, linetype="Sample Median"), color = "black")+
  scale_y_continuous(limits = c(-80, 35),breaks = seq(-80,35, by=10)) +
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
    title = "P cell class delta P1 - P4",
    x = "Probability",
    y = "Spike rate difference (spikes/sec)" )
ggsave(filename= "/Users/loicdaumail/Documents/Research/Maier_lab/adaptation_LGN_project/R/plots/p_bootstrapt_ydbt_trim15_asym95CI_deltaspikerate_density.svg")

##K density Plot
p_idx =org_data.filt$Unit.Index[org_data.filt$layer=="K"]
k_data = data.frame(p_idx,kdelta)
library(reshape2)
plot_dat <- melt(k_data, id.var=c('p_idx'))

ggplot(plot_dat, aes(y=value))+
  #geom_histogram(aes(x=..density..), alpha =.5)+
  geom_density(aes(y=value),color ="#93C5DE", fill ="#93C5DE", alpha =.3)+
  geom_hline(aes(yintercept=kdelta.m, linetype="Sample Mean"), color = "black")+
  geom_hline(aes(yintercept=kdelta.med, linetype="Sample Median"), color = "black")+
  scale_y_continuous(limits = c(-80, 35),breaks = seq(-80,35, by=10)) +
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
    title = "K cell class delta P1 - P4",
    x = "Probability",
    y = "Spike rate difference (spikes/sec)" )
ggsave(filename= "/Users/loicdaumail/Documents/Research/Maier_lab/adaptation_LGN_project/R/plots/k_bootstrapt_ydbt_trim0_asym95CI_deltaspikerate_density.svg")

######IGNORE THE REST AS MIGHT NOT BE RELEVANT############################

###### Repeat analysis with Cohen's d  ## this part doesn't seem valid (bootstrap-t on Cohen's d might not be feasible as it isn't the direct data)

mcohen = org_data.filt$effectsize[org_data.filt$layer=="M"]
pcohen = org_data.filt$effectsize[org_data.filt$layer=="P"]
kcohen = org_data.filt$effectsize[org_data.filt$layer=="K"]

##The following method works for a one-sample test, but need to use ydbt function to perform the correct bootstrap on dependent sample data.
#sample
#samp = pcohen
#samp.m <- mean(samp, trim = 0.15) # sample mean
#samp.sem <- trimse(samp, tr = 0.15) # sample estimation of the standard error of the mean
# centred data
#mcsamp <-  samp - mean(samp)
#n = length(samp)
#nboot <- 5000
#alpha <- 0.05
#tval <- function(x,nv,tr=0.15){
#  se <- sqrt(winvar(x,tr)) / ((1-2*tr) * sqrt(length(x)))
#  tval <- (mean(x,tr)-nv) / se
#  tval
#}

#bootT <- apply(matrix(sample(mcsamp, nboot*n, replace = TRUE), nrow = nboot), 1, tval, nv = 0)
#bootT.q <- quantile(bootT, probs = c(alpha/2, 1-alpha/2), type = 6)

# Asymmetric CI
#CI.lo <- round(samp.m - bootT.q[2] * samp.sem, digits = 2)
#CI.up <- round(samp.m - bootT.q[1] * samp.sem, digits = 2) 
#print(paste0("Asymmetric CI = [",CI.lo,", ",CI.up,"]"))

##M cells

p_idx =org_data.filt$Unit.Index[org_data.filt$layer=="M"]

m_data = data.frame(p_idx,mcohen,mpval)

#only get significant decreasing
m_data$sigdata[m_data$mpval<0.05 &mcohen<0] <- "Adapting Decrease"
m_data$sigdata[m_data$mpval>0.05] <- "Non Adapting"
m_data$sigdata[m_data$mpval<0.05 &mcohen>0] <- "Adapting Increase"

library(reshape2)
plot_dat <- melt(m_data, id.var=c('p_idx','sigdata','mpval'))

ggplot(plot_dat, aes(y=value))+
  geom_histogram(aes(fill=sigdata),bins=20)+
  scale_fill_manual(values=c("#93C5DE","#93C5DE","#B5B5B5"))+
  scale_x_discrete(limits=c(1,2,3,4))+
  geom_hline(aes(yintercept=CI.lo, linetype = "Bootstrap-T 95%CI"), color = "red")+
  geom_hline(aes(yintercept=CI.up,linetype="Bootstrap-T 95%CI"), color = "red")+
  geom_hline(aes(yintercept=samp.m, linetype="Sample Mean"), color = "black")+
  
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
    title = "M cell class delta P1 - P4 Cohen's d",
    x = "Unit Count",
    y = "Cohen's d" )
ggsave(filename= "/Users/loicdaumail/Documents/Research/adaptation_LGN_project/R/plots/m_bootstrapt_trim15_asym95CI_cohensd_histogram.png")

##P cells

p_idx =org_data.filt$Unit.Index[org_data.filt$layer=="P"]

p_data = data.frame(p_idx,pcohen,ppval)

#only get significant decreasing
p_data$sigdata[p_data$ppval<0.05 &pcohen<0] <- "Adapting Decrease"
p_data$sigdata[p_data$ppval>0.05] <- "Non Adapting"
p_data$sigdata[p_data$ppval<0.05 &pcohen>0] <- "Adapting Increase"

library(reshape2)
plot_dat <- melt(p_data, id.var=c('p_idx','sigdata','ppval'))

ggplot(plot_dat, aes(y=value))+
  geom_histogram(aes(fill=sigdata),bins=20)+
  scale_fill_manual(values=c("#F6999A","#F6999A","#B5B5B5"))+
  scale_x_discrete(limits=c(1,2,3,4))+
  geom_hline(aes(yintercept=CI.lo, linetype = "Bootstrap-T 95%CI"), color = "red")+
  geom_hline(aes(yintercept=CI.up,linetype="Bootstrap-T 95%CI"), color = "red")+
  geom_hline(aes(yintercept=samp.m, linetype="Sample Mean"), color = "black")+
  
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
    title = "P cell class delta P1 - P4 Cohen's d",
    x = "Unit Count",
    y = "Cohen's d" )
ggsave(filename= "/Users/loicdaumail/Documents/Research/adaptation_LGN_project/R/plots/p_bootstrapt_trim15_asym95CI_cohensd_histogram.png")

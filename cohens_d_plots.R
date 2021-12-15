library(R.matlab)
library("ggpubr")
library("colorspace")

setwd("~/LGN_data_042021/single_units/inverted_power_channels/good_single_units_data_4bumps_more/new_peak_alignment_anal/su_peaks_03032020_corrected/orig_peak_values/all_units")
bumps_data <- readMat("all_data_peaks.mat")
peaks =bumps_data[["peak.vals"]]
bfs_data <- readMat("bayesFactors.mat")
bfs =bfs_data[["bayesfactors"]]

#file_names contains both the filenames and the layer
file_names_data <- read.csv("filenames_layers.csv", header = TRUE,sep = ',', na.strings="")

#pvalues of adaptation of peak1 vs peak4
setwd("~/LGN_data_042021/single_units/inverted_power_channels/good_single_units_data_4bumps_more/new_peak_alignment_anal/lmer_results_peaks")
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
    
    sd_diff = sd(unit_peaks4-unit_peaks1)
    
    
    cd = (mean_pk4 - mean_pk1)/sd_diff
    
    effectsize[i]=cd
    peak1[i]=mean_pk1
    peak4[i]=mean_pk4
  }
}
channel_idx = t(t(rep(1:71, times =1)))

#data frame

org_data = data.frame("Unit Index"=channel_idx,"Cohen's D"=effectsize,"P-value"=pvalues_data[,5], "Peak1"=peak1,"Peak4"= peak4, file_names_data)


library(tidyr)
#no_na_data <- drop_na(org_data)

row.has.na <- apply(org_data, 1, function(x){any(is.na(x))})
org_data.filt <- org_data[!row.has.na,]


library(reshape2)
library(ggplot2)

## all cell classes effect sizes
mcohen = org_data.filt$Cohen.s.D[org_data.filt$layer=="M"]
pcohen = org_data.filt$Cohen.s.D[org_data.filt$layer=="P"]
kcohen = org_data.filt$Cohen.s.D[org_data.filt$layer=="K"]

cohen_dat <- c(mcohen,pcohen,kcohen)

##percent change
mpeak1 = org_data.filt$Peak1[org_data.filt$layer=="M"]
mpeak4 = org_data.filt$Peak4[org_data.filt$layer=="M"]
m_percent = -100*(mpeak1 - mpeak4)/mpeak1


ppeak1 = org_data.filt$Peak1[org_data.filt$layer=="P"]
ppeak4 = org_data.filt$Peak4[org_data.filt$layer=="P"]
p_percent = -100*(ppeak1 - ppeak4)/ppeak1


kpeak1 = org_data.filt$Peak1[org_data.filt$layer=="K"]
kpeak4 = org_data.filt$Peak4[org_data.filt$layer=="K"]
k_percent = -100*(kpeak1 - kpeak4)/kpeak1

# convert from wide to long

percent_dat <- c(m_percent,p_percent,k_percent)

m_cells = rep('M', times =18)
p_cells = rep('P', times = 15)
k_cells = rep('K', times= 3)

cell_class = c(m_cells, p_cells, k_cells)

all_cd = data.frame(cell_class, cohen_dat, percent_dat)

#Plot cds for all cell classes
plot
# plot
ggplot(all_cd) + 
  
  # box plots and jitter points, with modified x value
  geom_violin(aes(x=cell_class, y=cohen_dat,fill=cell_class), width=0.6, trim=F, color = NA)+
  scale_fill_manual(values=c("#A7D5A0","#93C5DE","#F69999"))+
  geom_boxplot(aes(x=cell_class, y=cohen_dat),width=0.1) +
  geom_jitter(aes(x=cell_class, y=cohen_dat),alpha =.5, position = position_jitter(width = 0)) +
  #scale_colour_gradient(high="#FAAE62", low="#333333" )+
  # specify x value order
  scale_x_discrete(limits=c('M', 'P', 'K'))+
  
  scale_y_continuous(limits = c(-15, 15),breaks = seq(-10,15, by=5)) +
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
    title = "Cell class comparison of P1 - P4 \n Cohen's D of spiking activity",
    x = "Cell Class",
    y = "Spike rate Cohen's D" )
#ggsave(filename= "C:/Users/daumail/Documents/LGN_data/single_units/inverted_power_channels/good_single_units_data_4bumps_more/new_peak_alignment_anal/su_peaks_03032020_corrected/orig_peak_values/plots/cell_class_plot_spikerate_percentchange.png")
#ggsave(filename= "C:/Users/daumail/Documents/LGN_data/single_units/inverted_power_channels/good_single_units_data_4bumps_more/new_peak_alignment_anal/su_peaks_03032020_corrected/orig_peak_values/plots/cell_class_plot_spikerate_percentchange.svg")



###plot just cell class per cell class effect size

#M data
#p_idx =org_data.filt$channel_idx[org_data.filt$layer=="M"]
p_idx =org_data.filt$Unit.Index[org_data.filt$layer=="M"]
colgrad <-abs(m_percent/100)
m_data = data.frame(p_idx,colgrad,mcohen)

plot_dat <- melt(m_data, id.var=c('p_idx','colgrad'))



ggplot(plot_dat) + 
  
  #scale_colour_gradientn(colours=heat.colors(50, alpha=.7, rev = TRUE))+
  # box plots and jitter points, with modified x value
  geom_violin(aes(x=variable, y=value,fill=variable), width=0.6, trim=F, color = NA)+
  scale_fill_manual(values=c("#93C5DE"))+
  geom_boxplot(aes(x=variable, y=value),width=0.1) +
  # simple lines
  #geom_jitter(alpha =.5, aes(x=paste0('dots ',variable), y=value,color=colgrad), position = position_jitter(width = 0)) +
  
  # specify x value order
  scale_x_discrete(limits=c('mcohen'))+
  scale_y_continuous(limits = c(.5, 1),breaks = seq(.5,1, by=.1)) +
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
    title = "M cell class comparison: percent change vs cohen's d",
    x = "Measure",
    y = "Spike rate (spikes/sec)" )

#### Histograms of effect sizes


#pvalues
mpval = org_data.filt$P.value[org_data.filt$layer=="M"]
ppval = org_data.filt$P.value[org_data.filt$layer=="P"]
kpval = org_data.filt$P.value[org_data.filt$layer=="K"]


##M cells

#t.test(mpeak4, mpeak1, paired = TRUE, alternative = "two.sided")
#qt(0.975,17)


p_idx =org_data.filt$Unit.Index[org_data.filt$layer=="M"]
colgrad <-abs(m_percent/100)
m_data = data.frame(p_idx,colgrad,mcohen,mpval)
mmcohen= mean(mcohen, trim = 0)
#only get significant decreasing
m_data$sigdata[m_data$mpval<0.05 &mcohen<0] <- "Adapting Decrease"
m_data$sigdata[m_data$mpval>0.05] <- "Non Adapting"
m_data$sigdata[m_data$mpval<0.05 &mcohen>0] <- "Adapting Increase"
#m_data[ which(m_data$mpval<0.05 & mcohen<0), ]

plot_dat <- melt(m_data, id.var=c('p_idx','colgrad','sigdata','mpval'))

ggplot(plot_dat, aes(y=value))+
  geom_histogram(aes(fill=sigdata),bins=20)+
  scale_fill_manual(values=c("#93C5DE","#93C5DE","#B5B5B5"))+
  scale_x_discrete(limits=c(1,2,3,4))+
  geom_hline(aes(yintercept=mmcohen, linetype="Sample Mean"), color = "black")+
  #geom_hline(yintercept=qt(0.975,18), linetype="dashed", color = "red")
  #geom_hline(yintercept=-(qt(0.975,9)/sqrt(10)), linetype="dashed", color = "red")
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
  )
ggsave(filename= "/Users/daumail/OneDrive - Vanderbilt/Documents/LGN_data_042021/single_units/inverted_power_channels/good_single_units_data_4bumps_more/new_peak_alignment_anal/plots/m_cells_cohenD_hist_mean.svg")

#P cells

p_idx =org_data.filt$Unit.Index[org_data.filt$layer=="P"]
colgrad <-abs(p_percent/100)
p_data = data.frame(p_idx,colgrad,pcohen,ppval)
mpcohen= mean(pcohen, trim = 0)
#only get significant decreasing
p_data$sigdata[p_data$ppval<0.05 &pcohen<0] <- "Adapting Decrease"
p_data$sigdata[p_data$ppval>0.05] <- "Non Adapting"
p_data$sigdata[p_data$ppval<0.05 &pcohen>0] <- "Adapting Increase"
#m_data[ which(m_data$mpval<0.05 & mcohen<0), ]

plot_dat <- melt(p_data, id.var=c('p_idx','colgrad','sigdata','ppval'))

ggplot(plot_dat, aes(y=value))+
  geom_histogram(aes(fill=sigdata),bins=20)+
  scale_fill_manual(values=c("#F6999A","#F6999A","#B5B5B5"))+
  scale_x_discrete(limits=c(1,2,3))+
  geom_hline(aes(yintercept=mpcohen, linetype="Sample Mean"), color = "black")+
  #geom_hline(yintercept=qt(0.975,18), linetype="dashed", color = "red")
  #geom_hline(yintercept=-(qt(0.975,9)/sqrt(10)), linetype="dashed", color = "red")
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
  )
ggsave(filename= "/Users/daumail/OneDrive - Vanderbilt/Documents/LGN_data_042021/single_units/inverted_power_channels/good_single_units_data_4bumps_more/new_peak_alignment_anal/plots/p_cells_cohenD_hist_mean.svg")

#K cells
p_idx =org_data.filt$Unit.Index[org_data.filt$layer=="K"]
colgrad <-abs(k_percent/100)
k_data = data.frame(p_idx,colgrad,kcohen,kpval)
mkcohen= mean(kcohen, trim = 0)
#only get significant decreasing
k_data$sigdata[k_data$kpval<0.05 &kcohen<0] <- "Adapting Decrease"
k_data$sigdata[k_data$kpval>0.05] <- "Non Adapting"
k_data$sigdata[k_data$kpval<0.05 &kcohen>0] <- "Adapting Increase"
#m_data[ which(m_data$mpval<0.05 & mcohen<0), ]

plot_dat <- melt(k_data, id.var=c('p_idx','colgrad','sigdata','kpval'))

ggplot(plot_dat, aes(y=value))+
  geom_histogram(aes(fill=sigdata),bins=20)+
  scale_fill_manual(values=c("#93C5DE","#B5B5B5"))+ #93C5DE #c("#A8D5A0","#B5B5B5"))+
  scale_x_discrete(limits=c(1))+
  geom_hline(aes(yintercept=mkcohen, linetype="Sample Mean"), color = "black")+
  #geom_hline(yintercept=qt(0.975,18), linetype="dashed", color = "red")
  #geom_hline(yintercept=-(qt(0.975,9)/sqrt(10)), linetype="dashed", color = "red")
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
  )
ggsave(filename= "/Users/daumail/OneDrive - Vanderbilt/Documents/LGN_data_042021/single_units/inverted_power_channels/good_single_units_data_4bumps_more/new_peak_alignment_anal/plots/k_cells_cohenD_hist_mean.svg")

## same plots for Bayes Factor ##############################

library(reshape2)
library(ggplot2)
#data frame

org_data = data.frame("Unit Index"=channel_idx,"Bayes Factor"=bfs,"Cohen's D"=effectsize,"P-value"=pvalues_data[,5], "Peak1"=peak1,"Peak4"= peak4, file_names_data)

library(tidyr)

row.has.na <- apply(org_data, 1, function(x){any(is.na(x))})
org_data.filt <- org_data[!row.has.na,]

## all cell classes bayes factors
mbayes = org_data.filt$Bayes.Factor[org_data.filt$layer=="M"]
pbayes = org_data.filt$Bayes.Factor[org_data.filt$layer=="P"]
kbayes = org_data.filt$Bayes.Factor[org_data.filt$layer=="K"]

#pvalues
mpval = org_data.filt$P.value[org_data.filt$layer=="M"]
ppval = org_data.filt$P.value[org_data.filt$layer=="P"]
kpval = org_data.filt$P.value[org_data.filt$layer=="K"]

##percent change
mpeak1 = org_data.filt$Peak1[org_data.filt$layer=="M"]
mpeak4 = org_data.filt$Peak4[org_data.filt$layer=="M"]
m_percent = -100*(mpeak1 - mpeak4)/mpeak1


ppeak1 = org_data.filt$Peak1[org_data.filt$layer=="P"]
ppeak4 = org_data.filt$Peak4[org_data.filt$layer=="P"]
p_percent = -100*(ppeak1 - ppeak4)/ppeak1


kpeak1 = org_data.filt$Peak1[org_data.filt$layer=="K"]
kpeak4 = org_data.filt$Peak4[org_data.filt$layer=="K"]
k_percent = -100*(kpeak1 - kpeak4)/kpeak1

## all cell classes effect sizes
mcohen = org_data.filt$Cohen.s.D[org_data.filt$layer=="M"]
pcohen = org_data.filt$Cohen.s.D[org_data.filt$layer=="P"]
kcohen = org_data.filt$Cohen.s.D[org_data.filt$layer=="K"]

##M cells

#t.test(mpeak4, mpeak1, paired = TRUE, alternative = "two.sided")
#qt(0.975,17)


p_idx =org_data.filt$Unit.Index[org_data.filt$layer=="M"]
#colgrad <-abs(m_percent/100)
m_data = data.frame(p_idx,log(mbayes),mpval)
mmbayes= mean(log(mbayes), trim = 0)
#only get significant decreasing
m_data$sigdata[m_data$mpval<0.05 &mcohen<0] <- "Adapting Decrease"
m_data$sigdata[m_data$mpval>0.05] <- "Non Adapting"
m_data$sigdata[m_data$mpval<0.05 &mcohen>0] <- "Adapting Increase"
#m_data[ which(m_data$mpval<0.05 & mcohen<0), ]

plot_dat <- melt(m_data, id.var=c('p_idx','sigdata','mpval'))

ggplot(plot_dat, aes(y=value))+
  geom_histogram(aes(fill=sigdata),bins=20)+
  scale_fill_manual(values=c("#93C5DE","#93C5DE","#B5B5B5"))+
  scale_x_discrete(limits=c(1,2,3,4))+
  geom_hline(aes(yintercept=mmbayes, linetype="Sample Mean"), color = "black")+
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
    title = paste("Bayes Factors \n of M cells"),
    x = "Unit Count",

    y = "Bayes Factor (log)" )
ggsave(filename= "/Users/daumail/OneDrive - Vanderbilt/Documents/LGN_data_042021/single_units/inverted_power_channels/good_single_units_data_4bumps_more/new_peak_alignment_anal/plots/m_cells_bayesfactor_hist_mean.svg")

#P cells

p_idx =org_data.filt$Unit.Index[org_data.filt$layer=="P"]
colgrad <-abs(p_percent/100)
p_data = data.frame(p_idx,log(pbayes),ppval)
mpbayes= mean(log(pbayes), trim = 0)
#only get significant decreasing
p_data$sigdata[p_data$ppval<0.05 &pcohen<0] <- "Adapting Decrease"
p_data$sigdata[p_data$ppval>0.05] <- "Non Adapting"
p_data$sigdata[p_data$ppval<0.05 &pcohen>0] <- "Adapting Increase"
#m_data[ which(m_data$mpval<0.05 & mcohen<0), ]

plot_dat <- melt(p_data, id.var=c('p_idx','sigdata','ppval'))

ggplot(plot_dat, aes(y=value))+
  geom_histogram(aes(fill=sigdata),bins=20)+
  scale_fill_manual(values=c("#F6999A","#F6999A","#B5B5B5"))+
  scale_x_discrete(limits=c(1,2,3))+
  geom_hline(aes(yintercept=mpbayes, linetype="Sample Mean"), color = "black")+
  #geom_hline(yintercept=qt(0.975,18), linetype="dashed", color = "red")
  #geom_hline(yintercept=-(qt(0.975,9)/sqrt(10)), linetype="dashed", color = "red")
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
    title = paste("Bayes Factors \n of P cells"),
    x = "Unit Count",
    y = "Bayes Factor (log)" )
ggsave(filename= "/Users/daumail/OneDrive - Vanderbilt/Documents/LGN_data_042021/single_units/inverted_power_channels/good_single_units_data_4bumps_more/new_peak_alignment_anal/plots/p_cells_bayesfactors_hist_mean.svg")


#K cells
p_idx =org_data.filt$Unit.Index[org_data.filt$layer=="K"]
k_data = data.frame(p_idx,log(kbayes),kpval)
mkbayes= mean(log(kbayes), trim = 0)
#only get significant decreasing
k_data$sigdata[k_data$kpval<0.05 &kcohen<0] <- "Adapting Decrease"
k_data$sigdata[k_data$kpval>0.05] <- "Non Adapting"
k_data$sigdata[k_data$kpval<0.05 &kcohen>0] <- "Adapting Increase"

plot_dat <- melt(k_data, id.var=c('p_idx','sigdata','kpval'))

ggplot(plot_dat, aes(y=value))+
  geom_histogram(aes(fill=sigdata),bins=20)+
  scale_fill_manual(values=c("#93C5DE","#B5B5B5"))+ #93C5DE #c("#A8D5A0","#B5B5B5"))+
  scale_x_discrete(limits=c(1))+
  geom_hline(aes(yintercept=mkbayes, linetype="Sample Mean"), color = "black")+
  #geom_hline(yintercept=qt(0.975,18), linetype="dashed", color = "red")
  #geom_hline(yintercept=-(qt(0.975,9)/sqrt(10)), linetype="dashed", color = "red")
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
    title = paste("Bayes Factors \n of K cells"),
    x = "Unit Count",
    y = "Bayes Factor (log)" )
ggsave(filename= "/Users/daumail/OneDrive - Vanderbilt/Documents/LGN_data_042021/single_units/inverted_power_channels/good_single_units_data_4bumps_more/new_peak_alignment_anal/plots/k_cells_bayesfactors_hist_mean.svg")


library(R.matlab)
library("ggpubr")
library("colorspace")
bumps_data <- readMat("all_norm_mean_data_peaks.mat")
peaks =bumps_data[["mean.peaks"]]
parts_data <- readMat("part1_part2_norm_power.mat")

#file_names contains both the filenames and the layer
file_names_data <- read.csv("filenames_layers.csv", header = TRUE,sep = ',', na.strings="")

peak1<- t(t(peaks[1,]))
peak4<- t(t(peaks[4,]))



#create channel index
channel_idx = t(t(rep(1:71, times =1)))

#data frame

org_data = data.frame(channel_idx,peak1, peak4, file_names_data)


library(tidyr)
#no_na_data <- drop_na(org_data)

row.has.na <- apply(org_data, 1, function(x){any(is.na(x))})
org_data.filt <- org_data[!row.has.na,]

##add normalized diff between peak 1 and peak 4 for colour gradient in plot multiplied by .045 for them to be between -1 and 1
org_data.filt$diff <- ((org_data.filt$peak1 - org_data.filt$peak4))/(min(org_data.filt$peak1 - org_data.filt$peak4))



## plot each cell class % change together


library(reshape2)
library(ggplot2)

## plot the all cell classes percent change
mpeak1 = org_data.filt$peak1[org_data.filt$layer=="M"]
mpeak4 = org_data.filt$peak4[org_data.filt$layer=="M"]
m_percent = -100*(mpeak1 - mpeak4)/mpeak1


ppeak1 = org_data.filt$peak1[org_data.filt$layer=="P"]
ppeak4 = org_data.filt$peak4[org_data.filt$layer=="P"]
p_percent = -100*(ppeak1 - ppeak4)/ppeak1


kpeak1 = org_data.filt$peak1[org_data.filt$layer=="K"]
kpeak4 = org_data.filt$peak4[org_data.filt$layer=="K"]
k_percent = -100*(kpeak1 - kpeak4)/kpeak1

# convert from wide to long

percent_dat <- c(m_percent,p_percent,k_percent)
m_cells = rep('M', times =18)
p_cells = rep('P', times = 15)
k_cells = rep('K', times= 3)

cell_class = c(m_cells, p_cells, k_cells)

all_percent = data.frame(cell_class, percent_dat)

plot
# plot
 ggplot(all_percent) + 

  # box plots and jitter points, with modified x value
  geom_violin(aes(x=cell_class, y=percent_dat,fill=cell_class), width=0.6, trim=F, color = NA)+
  scale_fill_manual(values=c("#A7D5A0","#93C5DE","#F69999"))+
  geom_boxplot(aes(x=cell_class, y=percent_dat),width=0.1) +
  geom_jitter(aes(x=cell_class, y=percent_dat),alpha =.5, position = position_jitter(width = 0)) +
  #scale_colour_gradient(high="#FAAE62", low="#333333" )+
  # specify x value order
  scale_x_discrete(limits=c('M', 'P', 'K'))+
  
  scale_y_continuous(limits = c(-70, 50),breaks = seq(-70,50, by=10)) +
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
    title = "Cell class comparison of P1 - P4 \n percent change of spiking activity",
    x = "Cell Class",
    y = "Spike rate % change" )
ggsave(filename= "C:/Users/daumail/Documents/LGN_data/single_units/inverted_power_channels/good_single_units_data_4bumps_more/new_peak_alignment_anal/su_peaks_03032020_corrected/orig_peak_values/plots/cell_class_plot_spikerate_percentchange.png")
ggsave(filename= "C:/Users/daumail/Documents/LGN_data/single_units/inverted_power_channels/good_single_units_data_4bumps_more/new_peak_alignment_anal/su_peaks_03032020_corrected/orig_peak_values/plots/cell_class_plot_spikerate_percentchange.svg")




## plot each cell class independently, with simple lines over violins


library(reshape2)
library(ggplot2)

## plot the P cell class
peak1 = org_data.filt$peak1[org_data.filt$layer=="P"]
peak4 = org_data.filt$peak4[org_data.filt$layer=="P"]
p_idx =org_data.filt$channel_idx[org_data.filt$layer=="P"]
p_diff = org_data.filt$diff[org_data.filt$layer == "P"]
p_peakdata = data.frame(p_idx,peak1,peak4, p_diff)

# convert from wide to long
plot_dat <- melt(p_peakdata, id.var=c('p_idx','p_diff'))

# plot
#svg('C:/Users/maier/Documents/LGN_data/single_units/inverted_power_channels/good_single_units_data_4bumps_more/plots/p_class_plot_peak1_peak4_normnegpos.svg')
ggplot(plot_dat) + 

  #scale_colour_gradientn(colours=heat.colors(50, alpha=.7, rev = TRUE))+
  # box plots and jitter points, with modified x value
  geom_violin(aes(x=paste0('dots ', variable), y=value,fill=variable), width=0.6, trim=F, color = NA)+
  scale_fill_manual(values=c("#93C5DE","#B4B4B4"))+
  geom_boxplot(aes(x=paste0('dots ', variable), y=value),width=0.1) +
  # simple lines
  geom_line(aes(x=paste0('dots ', variable), y=value, group=p_idx, color=p_diff)) +
  scale_colour_gradient(high="#93C5DE", low="#333333" )+
  geom_jitter(alpha =.5, aes(x=paste0('dots ',variable), y=value, color=p_diff), position = position_jitter(width = 0)) +
  
  # specify x value order
  scale_x_discrete(limits=c('dots peak1',
                            'dots peak4'))+
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
    title = "P cell class comparison: peak1 vs peak4 spiking activity",
    x = "Peak Number",
    y = "Spike rate (spikes/sec)" )
ggsave(filename= "C:/Users/daumail/Documents/LGN_data/single_units/inverted_power_channels/good_single_units_data_4bumps_more/new_peak_alignment_anal/su_peaks_03032020_corrected/orig_peak_values/plots/p_class_plot_mean_bs_peak1_peak4_violin_line_stacked.png")
ggsave(filename= "C:/Users/daumail/Documents/LGN_data/single_units/inverted_power_channels/good_single_units_data_4bumps_more/new_peak_alignment_anal/su_peaks_03032020_corrected/orig_peak_values/plots/p_class_plot_mean_bs_peak1_peak4_violin_line_stacked.svg")


##make plots with lines of varying width
#library(grid)
library(ggvwline)
library(ggnewscale)
library(reshape2)
library(ggplot2)

## plot the P cell class
ppeak1 = org_data.filt$peak1[org_data.filt$layer=="P"]
ppeak4 = org_data.filt$peak4[org_data.filt$layer=="P"]
p_idx =org_data.filt$channel_idx[org_data.filt$layer=="P"]
p_diff = org_data.filt$diff[org_data.filt$layer == "P"]
p_peakdata = data.frame(p_idx,ppeak1,ppeak4, p_diff)

#percent change, further used for the width
p_percent = -100*(ppeak1 - ppeak4)/ppeak1
w_first = rep(0,length(ppeak1))
w <- c(w_first,abs(p_percent/100))
# convert from wide to long
plot_dat <- melt(p_peakdata, id.var=c('p_idx','p_diff'))
plot_dat$width <-w
plot_dat$colgrad <-c(abs(p_percent/100),abs(p_percent/100))

# plot
#svg('C:/Users/maier/Documents/LGN_data/single_units/inverted_power_channels/good_single_units_data_4bumps_more/plots/p_class_plot_peak1_peak4_normnegpos.svg')
ggplot(mapping=aes(x= variable, y=value))+
  
     #scale_colour_gradientn(colours=heat.colors(50, alpha=.7, rev = TRUE))+
  # box plots and jitter points, with modified x value
   geom_violin(data = plot_dat, aes(fill=variable), width=0.6, trim=F, color = NA)+
   scale_fill_manual(values=c("#F6999A","#B4B4B4"))+
  
   geom_boxplot(data = plot_dat, width=0.1)+ 
  #lines with varying width
   new_scale_fill()+
   geom_vwline(data=plot_dat,aes(x= variable, y=value,group=p_idx, width=width,color = colgrad,fill = colgrad))+ #

   #geom_line(aes(x=paste0('dots ', variable), y=value, group=p_idx, color=p_diff)) +
   scale_fill_gradient(high="#F6999A", low="#333333" )+
  
   geom_jitter(data=plot_dat,alpha =.5, aes(color=colgrad), position = position_jitter(width = 0)) +
   scale_color_gradient(high="#F6999A", low="#333333" )+
   # specify x value order
   scale_x_discrete(limits=c('ppeak1', 'ppeak4'))+
 
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
     title = "P cell class comparison: peak1 vs peak4 spiking activity",
     x = "Peak Number",
     y = "Spike rate (Normalized)" )
ggsave(filename= "/Users/loicdaumail/Documents/Research/adaptation_LGN_project/R/plots/p_class_plot_mean_bs_peak1_peak4_violin_varrythickness.png")
ggsave(filename= "/Users/loicdaumail/Documents/Research/adaptation_LGN_project/R/plots/p_class_plot_mean_bs_peak1_peak4_violin_varrythickness.svg")



## plot the M cell class
mpeak1 = org_data.filt$peak1[org_data.filt$layer=="M"]
mpeak4 = org_data.filt$peak4[org_data.filt$layer=="M"]
m_idx =org_data.filt$channel_idx[org_data.filt$layer=="M"]
m_diff = org_data.filt$diff[org_data.filt$layer == "M"]
m_peakdata = data.frame(m_idx,mpeak1,mpeak4, m_diff)

#percent change, further used for the width
m_percent = -100*(mpeak1 - mpeak4)/mpeak1
w_first = rep(0,length(mpeak1))
w <- c(w_first,abs(m_percent/100))
# convert from wide to long
plot_dat <- melt(m_peakdata, id.var=c('m_idx','m_diff'))
plot_dat$width <-w
plot_dat$colgrad <-c(abs(m_percent/100),abs(m_percent/100))

# plot
#svg('C:/Users/maier/Documents/LGN_data/single_units/inverted_power_channels/good_single_units_data_4bumps_more/plots/p_class_plot_peak1_peak4_normnegpos.svg')
ggplot(mapping=aes(x= variable, y=value))+
  
  #scale_colour_gradientn(colours=heat.colors(50, alpha=.7, rev = TRUE))+
  # box plots and jitter points, with modified x value
  geom_violin(data = plot_dat, aes(fill=variable), width=0.6, trim=F, color = NA)+
  scale_fill_manual(values=c("#93C5DE","#B4B4B4"))+
  
  geom_boxplot(data = plot_dat, width=0.1)+ 
  #lines with varying width
  new_scale_fill()+
  geom_vwline(data=plot_dat,aes(x= variable, y=value,group=m_idx, width=width,color = colgrad,fill = colgrad))+ #
  
  #geom_line(aes(x=paste0('dots ', variable), y=value, group=p_idx, color=p_diff)) +
  scale_fill_gradient(high="#93C5DE", low="#333333" )+
  
  geom_jitter(data=plot_dat,alpha =.5, aes(color=colgrad), position = position_jitter(width = 0)) +
  scale_color_gradient(high="#93C5DE", low="#333333" )+
  # specify x value order
  scale_x_discrete(limits=c('mpeak1', 'mpeak4'))+
  
  scale_y_continuous(limits = c(.6, 1.05),breaks = seq(.6,1.05, by=.1)) +
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
    title = "M cell class comparison: peak1 vs peak4 spiking activity",
    x = "Peak Number",
    y = "Spike rate (Normalized)" )
ggsave(filename= "/Users/loicdaumail/Documents/Research/adaptation_LGN_project/R/plots/m_class_plot_mean_bs_peak1_peak4_violin_varrythickness.png")
ggsave(filename= "/Users/loicdaumail/Documents/Research/adaptation_LGN_project/R/plots/m_class_plot_mean_bs_peak1_peak4_violin_varrythickness.svg")

 
## plot the K cell class
mpeak1 = org_data.filt$peak1[org_data.filt$layer=="K"]
mpeak4 = org_data.filt$peak4[org_data.filt$layer=="K"]
m_idx =org_data.filt$channel_idx[org_data.filt$layer=="K"]
m_diff = org_data.filt$diff[org_data.filt$layer == "K"]
m_peakdata = data.frame(m_idx,mpeak1,mpeak4, m_diff)

#percent change, further used for the width
m_percent = -100*(mpeak1 - mpeak4)/mpeak1
w_first = rep(0,length(mpeak1))
w <- c(w_first,abs(m_percent/100))
# convert from wide to long
plot_dat <- melt(m_peakdata, id.var=c('m_idx','m_diff'))
plot_dat$width <-w
plot_dat$colgrad <-c(abs(m_percent/100),abs(m_percent/100))

# plot 
#svg('C:/Users/maier/Documents/LGN_data/single_units/inverted_power_channels/good_single_units_data_4bumps_more/plots/p_class_plot_peak1_peak4_normnegpos.svg')
ggplot(mapping=aes(x= variable, y=value))+
  
  #scale_colour_gradientn(colours=heat.colors(50, alpha=.7, rev = TRUE))+
  # box plots and jitter points, with modified x value
  geom_violin(data = plot_dat, aes(fill=variable), width=0.6, trim=F, color = NA)+
  scale_fill_manual(values=c("#A8D5A0","#B4B4B4"))+
  
  geom_boxplot(data = plot_dat, width=0.1)+ 
  #lines with varying width
  new_scale_fill()+
  geom_vwline(data=plot_dat,aes(x= variable, y=value,group=m_idx, width=width,color = colgrad,fill = colgrad))+ #
  
  scale_fill_gradient(high="#A8D5A0", low="#333333" )+
  
  geom_jitter(data=plot_dat,alpha =.5, aes(color=colgrad), position = position_jitter(width = 0)) +
  scale_color_gradient(high="#A8D5A0", low="#333333" )+
  # specify x value order
  scale_x_discrete(limits=c('mpeak1', 'mpeak4'))+
  
  scale_y_continuous(limits = c(.35, 1.05),breaks = seq(.3,1, by=.1)) +
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
    title = "k cell class comparison: peak1 vs peak4 spiking activity",
    x = "Peak Number",
    y = "Spike rate (Normalized)" )
ggsave(filename= "/Users/loicdaumail/Documents/Research/adaptation_LGN_project/R/plots/k_class_plot_mean_bs_peak1_peak4_violin_varrythickness.png")
ggsave(filename= "/Users/loicdaumail/Documents/Research/adaptation_LGN_project/R/plots/k_class_plot_mean_bs_peak1_peak4_violin_varrythickness.svg")

####Power plots #####

window1<- t(t(parts_data$parts[,1]))
window2<- t(t(parts_data$parts[,2]))
#create channel index
channel_idx = t(t(rep(1:71, times =1)))

#data frame

org_data = data.frame(channel_idx,window1, window2, file_names_data)


library(tidyr)
#no_na_data <- drop_na(org_data)

row.has.na <- apply(org_data, 1, function(x){any(is.na(x))})
org_data.filt <- org_data[!row.has.na,]
org_data.filt$diff <- ((org_data.filt$window1 - org_data.filt$window2))/(min(org_data.filt$window2 - org_data.filt$window1))

## plot the P cell class
window1 = org_data.filt$window1[org_data.filt$layer=="P"]
window2 = org_data.filt$window2[org_data.filt$layer=="P"]
p_idx =org_data.filt$channel_idx[org_data.filt$layer=="P"]
p_diff = org_data.filt$diff[org_data.filt$layer == "P"]

p_powerdata = data.frame(p_idx,window1,window2,p_diff)

#percent change, further used for the width
p_percent = -100*(window1 - window2)/window1
w_first = rep(0,length(window1))
w <- c(w_first,abs(p_percent/100))
# convert from wide to long
plot_dat <- melt(p_powerdata, id.var=c('p_idx','p_diff'))
plot_dat$width <-w
plot_dat$colgrad <-c(abs(p_percent/100),abs(p_percent/100))

# plot
#svg('C:/Users/maier/Documents/LGN_data/single_units/inverted_power_channels/good_single_units_data_4bumps_more/plots/p_class_plot_peak1_peak4_normnegpos.svg')
ggplot(mapping=aes(x= variable, y=value))+
  
  #scale_colour_gradientn(colours=heat.colors(50, alpha=.7, rev = TRUE))+
  # box plots and jitter points, with modified x value
  geom_violin(data = plot_dat, aes(fill=variable), width=0.6, trim=F, color = NA)+
  scale_fill_manual(values=c("#F6999A","#B4B4B4"))+
  
  geom_boxplot(data = plot_dat, width=0.1)+ 
  #lines with varying width
  new_scale_fill()+
  geom_vwline(data=plot_dat,aes(x= variable, y=value,group=p_idx, width=width,color = colgrad,fill = colgrad))+ #
  
  #geom_line(aes(x=paste0('dots ', variable), y=value, group=p_idx, color=p_diff)) +
  scale_fill_gradient(high="#F6999A", low="#333333" )+
  
  geom_jitter(data=plot_dat,alpha =.5, aes(color=colgrad), position = position_jitter(width = 0)) +
  scale_color_gradient(high="#F6999A", low="#333333" )+
  # specify x value order
  scale_x_discrete(limits=c('window1', 'window2'))+
  
  scale_y_continuous(limits = c(.3, 1),breaks = seq(.3,1, by=.1)) +
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
    title = "P cell class power comparison: window1 vs window 2",
    x = "Window",
    y = "Power (Normalized)" )
ggsave(filename= "/Users/loicdaumail/Documents/Research/adaptation_LGN_project/R/plots/p_class_plot_mean_power_violin_varrythickness.png")
ggsave(filename= "/Users/loicdaumail/Documents/Research/adaptation_LGN_project/R/plots/p_class_plot_mean_power_violin_varrythickness.svg")


## plot the M cell class
window1 = org_data.filt$window1[org_data.filt$layer=="M"]
window2 = org_data.filt$window2[org_data.filt$layer=="M"]
p_idx =org_data.filt$channel_idx[org_data.filt$layer=="M"]
p_diff = org_data.filt$diff[org_data.filt$layer == "M"]

p_powerdata = data.frame(p_idx,window1,window2,p_diff)

#percent change, further used for the width
p_percent = -100*(window1 - window2)/window1
w_first = rep(0,length(window1))
w <- c(w_first,abs(p_percent/100))
# convert from wide to long
plot_dat <- melt(p_powerdata, id.var=c('p_idx','p_diff'))
plot_dat$width <-w
plot_dat$colgrad <-c(abs(p_percent/100),abs(p_percent/100))

# plot
#svg('C:/Users/maier/Documents/LGN_data/single_units/inverted_power_channels/good_single_units_data_4bumps_more/plots/p_class_plot_peak1_peak4_normnegpos.svg')
ggplot(mapping=aes(x= variable, y=value))+
  
  #scale_colour_gradientn(colours=heat.colors(50, alpha=.7, rev = TRUE))+
  # box plots and jitter points, with modified x value
  geom_violin(data = plot_dat, aes(fill=variable), width=0.6, trim=F, color = NA)+
  scale_fill_manual(values=c("#93C5DE","#B4B4B4"))+
  
  geom_boxplot(data = plot_dat, width=0.1)+ 
  #lines with varying width
  new_scale_fill()+
  geom_vwline(data=plot_dat,aes(x= variable, y=value,group=p_idx, width=width,color = colgrad,fill = colgrad))+ #
  
  #geom_line(aes(x=paste0('dots ', variable), y=value, group=p_idx, color=p_diff)) +
  scale_fill_gradient(high="#93C5DE", low="#333333" )+
  
  geom_jitter(data=plot_dat,alpha =.5, aes(color=colgrad), position = position_jitter(width = 0)) +
  scale_color_gradient(high="#93C5DE", low="#333333" )+
  # specify x value order
  scale_x_discrete(limits=c('window1', 'window2'))+
  
  scale_y_continuous(limits = c(.3, 1),breaks = seq(.3,1, by=.1)) +
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
    title = "M cell class power comparison: window1 vs window 2",
    x = "Window",
    y = "Power (Normalized)" )
ggsave(filename= "/Users/loicdaumail/Documents/Research/adaptation_LGN_project/R/plots/m_class_plot_mean_power_violin_varrythickness.png")
ggsave(filename= "/Users/loicdaumail/Documents/Research/adaptation_LGN_project/R/plots/m_class_plot_mean_power_violin_varrythickness.svg")

## plot the K cell class
window1 = org_data.filt$window1[org_data.filt$layer=="K"]
window2 = org_data.filt$window2[org_data.filt$layer=="K"]
p_idx =org_data.filt$channel_idx[org_data.filt$layer=="K"]
p_diff = org_data.filt$diff[org_data.filt$layer == "K"]

p_powerdata = data.frame(p_idx,window1,window2,p_diff)

#percent change, further used for the width
p_percent = -100*(window1 - window2)/window1
w_first = rep(0,length(window1))
w <- c(w_first,abs(p_percent/100))
# convert from wide to long
plot_dat <- melt(p_powerdata, id.var=c('p_idx','p_diff'))
plot_dat$width <-w
plot_dat$colgrad <-c(abs(p_percent/100),abs(p_percent/100))

# plot
#svg('C:/Users/maier/Documents/LGN_data/single_units/inverted_power_channels/good_single_units_data_4bumps_more/plots/p_class_plot_peak1_peak4_normnegpos.svg')
ggplot(mapping=aes(x= variable, y=value))+
  
  #scale_colour_gradientn(colours=heat.colors(50, alpha=.7, rev = TRUE))+
  # box plots and jitter points, with modified x value
  geom_violin(data = plot_dat, aes(fill=variable), width=0.6, trim=F, color = NA)+
  scale_fill_manual(values=c("#A8D5A0","#B4B4B4"))+
  
  geom_boxplot(data = plot_dat, width=0.1)+ 
  #lines with varying width
  new_scale_fill()+
  geom_vwline(data=plot_dat,aes(x= variable, y=value,group=p_idx, width=width,color = colgrad,fill = colgrad))+ #
  
  #geom_line(aes(x=paste0('dots ', variable), y=value, group=p_idx, color=p_diff)) +
  scale_fill_gradient(high="#A8D5A0", low="#333333" )+
  
  geom_jitter(data=plot_dat,alpha =.5, aes(color=colgrad), position = position_jitter(width = 0)) +
  scale_color_gradient(high="#A8D5A0", low="#333333" )+
  # specify x value order
  scale_x_discrete(limits=c('window1', 'window2'))+
  
  scale_y_continuous(limits = c(.1, 1),breaks = seq(.1,1, by=.1)) +
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
    title = "K cell class power comparison: window1 vs window 2",
    x = "Window",
    y = "Power (Normalized)" )
ggsave(filename= "/Users/loicdaumail/Documents/Research/adaptation_LGN_project/R/plots/k_class_plot_mean_power_violin_varrythickness.png")
ggsave(filename= "/Users/loicdaumail/Documents/Research/adaptation_LGN_project/R/plots/k_class_plot_mean_power_violin_varrythickness.svg")
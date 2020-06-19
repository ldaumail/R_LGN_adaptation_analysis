library(R.matlab)

bumps_data <- readMat("all_raw_mean_data_peaks.mat")
peaks =bumps_data[["mean.peaks"]]
#parts_data <- readMat("part1_part2_norm_power.mat")

#file_names contains both the filenames and the layer
file_names_data <- read.csv("filenames_layers.csv", header = TRUE,sep = ',', na.strings="")

#for (i in 1:length(peaks))
peak1<- t(t(peaks[1,]))
peak4<- t(t(peaks[4,]))

percentch = -100*(peak1 - peak4)/peak1
pdiff = peak1-peak4

channel_idx = t(t(rep(1:71, times =1)))

#data frame

org_data = data.frame(channel_idx,peak1, peak4, percentch, pdiff,file_names_data)


library(tidyr)
#no_na_data <- drop_na(org_data)

row.has.na <- apply(org_data, 1, function(x){any(is.na(x))})
org_data.filt <- org_data[!row.has.na,]

mpeak1 = org_data.filt$peak1[org_data.filt$layer=="M"]
mpeak4 = org_data.filt$peak4[org_data.filt$layer=="M"]
mdf = data.frame("Peak1"=mpeak1, "Peak4"=mpeak4)
mdf.long <- mdf %>%
  gather(key = "Peak", value = "Activity", Peak1, Peak4)

ppeak1 = org_data.filt$peak1[org_data.filt$layer=="P"]
ppeak4 = org_data.filt$peak4[org_data.filt$layer=="P"]
pdf = data.frame("Peak1"=ppeak1, "Peak4"=ppeak4)
pdf.long <- pdf %>%
  gather(key = "Peak", value = "Activity", Peak1, Peak4)

kpeak1 = org_data.filt$peak1[org_data.filt$layer=="K"]
kpeak4 = org_data.filt$peak4[org_data.filt$layer=="K"]
kdf = data.frame("Peak1"=kpeak1, "Peak4"=kpeak4)
kdf.long <- kdf %>%
  gather(key = "Peak", value = "Activity", Peak1, Peak4)

t.test(mpeak4, mpeak1, paired = TRUE, alternative = "two.sided")
t.test(ppeak4, ppeak1, paired = TRUE, alternative = "two.sided")
t.test(kpeak4, kpeak1, paired = TRUE, alternative = "two.sided")

#non parametric test
WM =wilcox.test(mdf.long$Activity ~ mdf.long$Peak, paired = TRUE)
WP =wilcox.test(pdf.long$Activity ~ pdf.long$Peak, paired = TRUE)
WK =wilcox.test(kdf.long$Activity ~ kdf.long$Peak, paired = TRUE)
Zm = qnorm(WM$p.value/2)
Zp = qnorm(WP$p.value/2)
Zk = qnorm(WK$p.value/2)

#compute medians (not essential)
library(rstatix)
mdf.long %>%
  group_by(Peak) %>%
  get_summary_stats(Activity, type = "median_iqr")

ks.test(org_data.filt$percentch[org_data.filt$layer=="M"], "pnorm", mean=mean(org_data.filt$percentch[org_data.filt$layer=="M"]), sd=sd(org_data.filt$percentch[org_data.filt$layer=="M"]))
ks.test(org_data.filt$percentch[org_data.filt$layer=="P"], "pnorm", mean=mean(org_data.filt$percentch[org_data.filt$layer=="P"]), sd=sd(org_data.filt$percentch[org_data.filt$layer=="P"]))
ks.test(org_data.filt$percentch[org_data.filt$layer=="K"], "pnorm", mean=mean(org_data.filt$percentch[org_data.filt$layer=="K"]), sd=sd(org_data.filt$percentch[org_data.filt$layer=="K"]))


ks.test(org_data.filt$pdiff[org_data.filt$layer=="M"], "pnorm", mean=mean(org_data.filt$pdiff[org_data.filt$layer=="M"]), sd=sd(org_data.filt$pdiff[org_data.filt$layer=="M"]))
ks.test(org_data.filt$pdiff[org_data.filt$layer=="P"], "pnorm", mean=mean(org_data.filt$pdiff[org_data.filt$layer=="P"]), sd=sd(org_data.filt$pdiff[org_data.filt$layer=="P"]))
ks.test(org_data.filt$pdiff[org_data.filt$layer=="K"], "pnorm", mean=mean(org_data.filt$pdiff[org_data.filt$layer=="K"]), sd=sd(org_data.filt$pdiff[org_data.filt$layer=="K"]))


#same for analysis of the power

file_names_data <- read.csv("filenames_layers.csv", header = TRUE,sep = ',', na.strings="")
parts_data <- readMat("part1_part2_norm_power.mat")
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
mwin1 = org_data.filt$window1[org_data.filt$layer=="M"]
mwin2 = org_data.filt$window2[org_data.filt$layer=="M"]
mdf = data.frame("Window1"=mwin1, "Window2"=mwin2)
mdf.long <- mdf %>%
  gather(key = "Window", value = "Power", Window1, Window2)

pwin1 = org_data.filt$window1[org_data.filt$layer=="P"]
pwin2 = org_data.filt$window2[org_data.filt$layer=="P"]
pdf = data.frame("Window1"=pwin1, "Window2"=pwin2)
pdf.long <- pdf %>%
  gather(key = "Window", value = "Power", Window1, Window2)

kwin1 = org_data.filt$window1[org_data.filt$layer=="K"]
kwin2 = org_data.filt$window2[org_data.filt$layer=="K"]
kdf = data.frame("Window1"=kwin1, "Window2"=kwin2)
kdf.long <- kdf %>%
  gather(key = "Window", value = "Power", Window1, Window2)

#non parametric test
WM =wilcox.test(mdf.long$Power ~ mdf.long$Window, paired = TRUE)
WP =wilcox.test(pdf.long$Power ~ pdf.long$Window, paired = TRUE)
WK =wilcox.test(kdf.long$Power ~ kdf.long$Window, paired = TRUE)
Zm = qnorm(WM$p.value/2)
Zp = qnorm(WP$p.value/2)
Zk = qnorm(WK$p.value/2)

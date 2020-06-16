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
t.test(mpeak4, mpeak1, paired = TRUE, alternative = "two.sided")


ppeak1 = org_data.filt$peak1[org_data.filt$layer=="P"]
ppeak4 = org_data.filt$peak4[org_data.filt$layer=="P"]
t.test(ppeak4, ppeak1, paired = TRUE, alternative = "two.sided")


kpeak1 = org_data.filt$peak1[org_data.filt$layer=="K"]
kpeak4 = org_data.filt$peak4[org_data.filt$layer=="K"]
t.test(kpeak4, kpeak1, paired = TRUE, alternative = "two.sided")

ks.test(org_data.filt$percentch[org_data.filt$layer=="M"], "pnorm", mean=mean(org_data.filt$percentch[org_data.filt$layer=="M"]), sd=sd(org_data.filt$percentch[org_data.filt$layer=="M"]))
ks.test(org_data.filt$percentch[org_data.filt$layer=="P"], "pnorm", mean=mean(org_data.filt$percentch[org_data.filt$layer=="P"]), sd=sd(org_data.filt$percentch[org_data.filt$layer=="P"]))
ks.test(org_data.filt$percentch[org_data.filt$layer=="K"], "pnorm", mean=mean(org_data.filt$percentch[org_data.filt$layer=="K"]), sd=sd(org_data.filt$percentch[org_data.filt$layer=="K"]))


ks.test(org_data.filt$pdiff[org_data.filt$layer=="M"], "pnorm", mean=mean(org_data.filt$pdiff[org_data.filt$layer=="M"]), sd=sd(org_data.filt$pdiff[org_data.filt$layer=="M"]))
ks.test(org_data.filt$pdiff[org_data.filt$layer=="P"], "pnorm", mean=mean(org_data.filt$pdiff[org_data.filt$layer=="P"]), sd=sd(org_data.filt$pdiff[org_data.filt$layer=="P"]))
ks.test(org_data.filt$pdiff[org_data.filt$layer=="K"], "pnorm", mean=mean(org_data.filt$pdiff[org_data.filt$layer=="K"]), sd=sd(org_data.filt$pdiff[org_data.filt$layer=="K"]))
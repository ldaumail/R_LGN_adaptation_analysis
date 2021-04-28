library(R.matlab)

org_data = read.csv("C:/Users/daumail/OneDrive - Vanderbilt/Documents/LGN_data_042021/single_units/adaptation_index/data/AdaptationIndexData.csv", header = TRUE,sep = ',', na.strings="")

library(tidyr)
#no_na_data <- drop_na(org_data)
row.has.na <- apply(org_data, 1, function(x){any(is.na(x))})
org_data.filt <- org_data[!row.has.na,]

#all Idx data mono/bino
IdxMono = org_data.filt$Adaptation.Index[org_data.filt$Condition =="Monocular"]
IdxBino = org_data.filt$Adaptation.Index[org_data.filt$Condition =="Binocular"]

#Idx data mono/bino per cell class
mIdxMono = org_data.filt$Adaptation.Index[org_data.filt$Cell.Class=="M" & org_data.filt$Condition =="Monocular"]
mIdxBino = org_data.filt$Adaptation.Index[org_data.filt$Cell.Class=="M" & org_data.filt$Condition =="Binocular"]

pIdxMono = org_data.filt$Adaptation.Index[org_data.filt$Cell.Class=="P" & org_data.filt$Condition =="Monocular"]
pIdxBino = org_data.filt$Adaptation.Index[org_data.filt$Cell.Class=="P" & org_data.filt$Condition =="Binocular"]

# dependencies
#library(ggplot2)
#library(tibble)
library("ggpubr")
library("colorspace")
source('C:/Users/daumail/OneDrive - Vanderbilt/Documents/R/bootstrap_functions/theme_gar.txt')
source('C:/Users/daumail/OneDrive - Vanderbilt/Documents/R/bootstrap_functions/fun.txt')
#library(beepr)
library(WRS)
#library(WRS2)

##First option: using the WRS package
yuend1<-WRS::yuend(IdxMono,IdxBino,tr=0,alpha=.05)

#here no trimming on M cell class as skewness is <.5
#mbootstraptrim <- WRS::ydbt(mpeak4,mpeak1,tr=.15,alpha=0.05,nboot=1000,side=FALSE,plotit=TRUE,op=1)
mbootstrap <- WRS::ydbt(mIdxMono,mIdxBino,tr=0.15,alpha=0.05,nboot=5000,side=FALSE,plotit=TRUE,op=1)
#trimci(mpeak1,tr=.2,alpha=.05,null.value=0)

#here trimming = 0.15 on P cell class as skewness is >0.5
pbootstraptrim <- WRS::ydbt(pIdxMono,pIdxBino,tr=0,alpha=0.05,nboot=5000,side=FALSE,plotit=TRUE,op=1)

#pbootstrap <- WRS::ydbt(ppeak4,ppeak1,tr=0,alpha=0.05,nboot=5000,side=FALSE,plotit=TRUE,op=1)


## make different sorts of plots

mIdxdiff = mIdxMono - pIdxBino


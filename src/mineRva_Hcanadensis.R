#How to use mineRva functions

setwd("C:/Users/bty/OneDrive - UNCG/GitHub/mineRva")
source(minervaTools.R)

#data<-read.csv("data/180517_Sample_pos_neg_sn3.csv",check.names="FALSE")
data<-read.csv("data/180517_Sample_pos_neg_sn10.csv",check.names="FALSE")

sample_group=sgroup(data,c("MeOH","03125","0625","125","250","500"))
dataset<-dataSet(data,sample_group,rt="rtmed",mz="mzmed")
x_method1<-xM1(dataset,sample_group)
msr<-msr(xm1=x_method1,sample_group)
outliers<-exclude5(xm1=x_method1,msr,T=155) #optional
extract_ions<-ionx(xm1=x_method1,sgrp=sample_group,rtgrp=155,a=0.8) #OPTIONAL: ,exclude=outliers
x_method2<-xM2(x=dataset,sample_group)
relate_ions<-ionr(xm2=x_method2,inx=extract_ions,sample_group,rtgrp=155, a=0.8) #OPTIONAL: ,exclude=outliers

ion_table<-ionTable(extract_ions,relate_ions)
#write.csv(summary, file="T155_summary.csv",row.names = FALSE)
beta_table<-betaTable(extract_ions,relate_ions)
ion_heatmap<-ionHeatmap(extract_ions,relate_ions)
extracted_ion_plots<-inxPlot(xm1=x_method1,inx=extract_ions,sgrp=sample_group)
related_ion_plots<-inrPlot(xm2=x_method2,inr=relate_ions,sgrp=sample_group,frac=0.2) 


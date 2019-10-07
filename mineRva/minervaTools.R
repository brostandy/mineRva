#mineRva: Metabolite Ions Extraction and Visualization data analysis workflow (date 3.19)
#This source is the compilation of tools for mineRva where users can call the .R file with
#source(minervaTools.R)
#
#citation
#If you are use and/or modify this tool for research publication, please cite:
#Article Title: Botanical metabolite ions extraction from full electrospray ionization mass spectrometry using high-dimensional penalized regression
#Rostandy, B. & Gao, X. Metabolomics (2019) 15: 136. https://doi.org/10.1007/s11306-019-1603-5


#Additional citation for mineRva requirements:
citation(package="glmnet")
citation(package="RColorBrewer")
citation(package="gplots")

#group X into sample groups (if applicable)
#INPUT: x=data, groupnames=5 groups and 1 blank, blank is always indicated first 
#OUTPUT: samples in the groups, and number of sample in the groups
sgroup<-function(x,groupnames){ 
  group=list()
  nsamp_in_group=list()
  for (i in 1:length(groupnames)){
    group[[i]]<-grep((paste('^', groupnames[i],sep="")),colnames(x),value=TRUE)
    nsamp_in_group[[i]]<-length(group[[i]])
  }
  return(list(group,nsamp_in_group))
}

#Get required information from data entry and replace zero-entries
#INPUT: x= data, sgrp=sgroup() object, rt=colname for retention time ("string"), 
  #mz=colname for mass-to-charge ("string"), 
  #width=width of ret. time data to be analyzed (in 1 sec increment; 1 second is default)
#OUTPUT: data matrix with rtgroup, rtime, mzmed, rtmed and ion intensity with replacement for zero-entries
dataSet<-function(x,sgrp,width=1,rt,mz)  { #x= data, sgrp=get_group() object, width=integer (unit in sec)
  #rt=colname for retention time ("string"), mz=colname for mass-to-charge ("string")
  comb_samp_grp=Reduce(c,sgrp[[1]])
  data_m<-as.matrix(x[,comb_samp_grp])
  p<-nrow(data_m) 
  for (i in 1:p) {
    xsort<- sort(data_m[i,])
    loc<- which(xsort>0)[1] 
    data_m[i, which(data_m[i,]==0)]= xsort[loc]/1000
    
    rm(xsort)
  } 
  rtmed=(x[,rt])
  w=width #unit in seconds
  rtgroup=round(rtmed/w,0) #to enable sub dataset based on user preference
  rtime=round(rtmed,0)
  mzmed=x[,mz]
  
  return(cbind(rtgroup,rtime,mzmed,rtmed,data_m))
}

#========================Data Tranformation=========================#
#INPUT: sgrp is sgroup()object, x is dataSet() object
#OUTPUT: xM0 returns standardized data matrix
######## xM1 returns log2 ratio of ion intensity to mean intensity of blanks
######## xM2 returns log2 ratio of ion intensity to median intensity of all samples

#Rescale and Transform dataset (1) with sd (standardized dataset)
xM0<-function(x,sgrp){ 
  comb_samp_grp=Reduce(c,sgrp[[1]])
  sd<-apply(x[,comb_samp_grp],1,sd) #get std. dev. for X by row
  mean<-apply(x[,comb_samp_grp],1,mean) #get mean for x by row
  X<-((x[,comb_samp_grp]-mean)/sd) #standardized dataset
  rtgroup=x[,"rtgroup"]
  rtime=x[,"rtime"] #or rtmed=(x[,rt]) and rtime=round(rtmed,0)
  mzmed=x[,"mzmed"]
  rtmed=x[,"rtmed"]
  return(cbind(rtgroup,rtime,mzmed,rtmed,X))
}

#Rescale and Transform dataset (2) against blanks (fold change from mean of blank)
xM1<-function(x,sgrp){ 
  MeOH_mean<-apply(x[,sgrp[[1]][[1]]],1,mean) #calculate mean intensity for blank
  Sample<-Reduce(c,sgrp[[1]][2:length(sgrp[[1]])])
  X<-log2(x[,Sample]/MeOH_mean) #X for Model1
  rtgroup=x[,"rtgroup"]
  rtime=x[,"rtime"]
  mzmed=x[,"mzmed"]
  rtmed=x[,"rtmed"]
  return(cbind(rtgroup,rtime,mzmed,rtmed,X))
}

#Rescale and Transform dataset (3) with median
xM2<-function(x,sgrp){ 
  comb_samp_grp=Reduce(c,sgrp[[1]])
  p_median<-apply(x[,comb_samp_grp],1,median) #get median for all group of each row
  X<-log2(x[,comb_samp_grp]/p_median) #X for Model2
  rtgroup=x[,"rtgroup"]
  rtime=x[,"rtime"]
  mzmed=x[,"mzmed"]
  rtmed=x[,"rtmed"]
  return(cbind(rtgroup,rtime,mzmed,rtmed,X))
}

#Get mean, standard deviation and relative standard deviation for extract
#Also used for rsd cutoff, such as outlier removal
#INPUT: xm1 is xM1() object; it can also be used for any xM() object, but only msr of xM1 can be used for outliers detection
      #and sgrp is sgroup() object
#OUTPUT: mean, standard dev., and relative std. dev. for any transformed dataset
msr<-function(xm1,sgrp){ 
  mean=sd=rsd=list()
  for (k in 2:length(sgrp[[1]])){
    mean[[k]]<-apply(xm1[,sgrp[[1]][[k]]],1,mean)
    sd[[k]]<-apply(xm1[,sgrp[[1]][[k]]],1,sd)
    rsd[[k]]<-sd[[k]]/mean[[k]]
  }
  return(list(mean,sd,rsd))
}

#Get outlier for dataset that have five sample groups (exclude blanks)
#INPUT: xm1 = xM1() object, msr = msr() of xM1() object, rtgrp is ret. time group to be analyzed
#OUTPUT: relative std. dev from each group, outlier from each group, exclude list (following the 2/5 rules)
exclude5<-function(xm1,msr,rtgrp){  
nrsd=length(msr[[3]])-length(is.null(msr[[3]]))
rsd=Reduce(c,msr[[3]])
nrow=length(rsd)/nrsd
rsd_df<-data.frame(matrix(c(xm1[,"rtgroup"],rsd),nrow,nrsd+1))
colnames(rsd_df)<-c("rtgroup","rsd1","rsd2","rsd3","rsd4","rsd5")
y = as.matrix(subset(rsd_df, rtgroup == rtgrp, select = c(rsd1,rsd2,rsd3,rsd4,rsd5)))
rsd.boxplot<-boxplot(y,col=c("lightcoral","orchid3","springgreen3","steelblue3","peru"), 
                     names=c(0.03125, 0.0625, 0.125,0.250, 0.500), plot=FALSE) 
#get the Q1 and Q3
rsd1.Q1<-rsd.boxplot$stats[2,1] ; rsd1.Q3<-rsd.boxplot$stats[4,1]
rsd2.Q1<-rsd.boxplot$stats[2,2] ; rsd2.Q3<-rsd.boxplot$stats[4,2]
rsd3.Q1<-rsd.boxplot$stats[2,3] ; rsd3.Q3<-rsd.boxplot$stats[4,3]
rsd4.Q1<-rsd.boxplot$stats[2,4] ; rsd4.Q3<-rsd.boxplot$stats[4,4]
rsd5.Q1<-rsd.boxplot$stats[2,5] ; rsd5.Q3<-rsd.boxplot$stats[4,5]
#get the cutout regions and return indices
rsd1_out<-which((y[,'rsd1']>(rsd1.Q3+1.5*(rsd1.Q3-rsd1.Q1))) | (y[,'rsd1']<(rsd1.Q1-1.5*(rsd1.Q3-rsd1.Q1))))
rsd2_out<-which((y[,'rsd2']>(rsd2.Q3+1.5*(rsd2.Q3-rsd2.Q1))) | (y[,'rsd2']<(rsd2.Q1-1.5*(rsd2.Q3-rsd2.Q1))))
rsd3_out<-which((y[,'rsd3']>(rsd3.Q3+1.5*(rsd3.Q3-rsd3.Q1))) | (y[,'rsd3']<(rsd3.Q1-1.5*(rsd3.Q3-rsd3.Q1))))
rsd4_out<-which((y[,'rsd4']>(rsd4.Q3+1.5*(rsd4.Q3-rsd4.Q1))) | (y[,'rsd4']<(rsd4.Q1-1.5*(rsd4.Q3-rsd4.Q1))))
rsd5_out<-which((y[,'rsd5']>(rsd5.Q3+1.5*(rsd5.Q3-rsd5.Q1))) | (y[,'rsd5']<(rsd5.Q1-1.5*(rsd5.Q3-rsd5.Q1))))

#give conditional for outlier removal (for 5 sample group + 1 neg. control)
all=Reduce(intersect, list(names(rsd1_out),names(rsd2_out),names(rsd3_out),names(rsd4_out),names(rsd5_out)))
abcd=Reduce(intersect, list(names(rsd1_out),names(rsd2_out),names(rsd3_out),names(rsd4_out))) 
abce=Reduce(intersect, list(names(rsd1_out),names(rsd2_out),names(rsd3_out),names(rsd5_out))) 
abde=Reduce(intersect, list(names(rsd1_out),names(rsd2_out),names(rsd4_out),names(rsd5_out))) 
acde=Reduce(intersect, list(names(rsd1_out),names(rsd3_out),names(rsd4_out),names(rsd5_out))) 
bcde=Reduce(intersect, list(names(rsd2_out),names(rsd3_out),names(rsd4_out),names(rsd5_out))) 
abc=intersect(intersect(names(rsd1_out),names(rsd2_out)),names(rsd3_out)) 
abd=intersect(intersect(names(rsd1_out),names(rsd2_out)),names(rsd4_out)) 
abe=intersect(intersect(names(rsd1_out),names(rsd2_out)),names(rsd5_out)) 
acd=intersect(intersect(names(rsd1_out),names(rsd3_out)),names(rsd4_out)) 
ace=intersect(intersect(names(rsd1_out),names(rsd3_out)),names(rsd5_out))
ade=intersect(intersect(names(rsd1_out),names(rsd4_out)),names(rsd5_out))
bcd=intersect(intersect(names(rsd2_out),names(rsd3_out)),names(rsd4_out))
bce=intersect(intersect(names(rsd2_out),names(rsd3_out)),names(rsd5_out))
bde=intersect(intersect(names(rsd2_out),names(rsd4_out)),names(rsd5_out))
cde=intersect(intersect(names(rsd3_out),names(rsd4_out)),names(rsd5_out)) 
ac=intersect(names(rsd1_out),names(rsd3_out))
ad=intersect(names(rsd1_out),names(rsd4_out))
bc=intersect(names(rsd2_out),names(rsd3_out))
bd=intersect(names(rsd2_out),names(rsd4_out))
be=intersect(names(rsd2_out),names(rsd5_out))
cd=intersect(names(rsd3_out),names(rsd4_out))
ce=intersect(names(rsd3_out),names(rsd5_out))
exclude=as.numeric(unique(c(all,abcd,abce,abde,acde,bcde,abc,abd,abe,acd,ace,ade,bcd,bce,bde,cde,ac,ad,
                            bc,bd,be,cd,ce)))

return(list(boxplot_stats=rsd.boxplot$stats, rsd1_out=rsd1_out, rsd2_out=rsd2_out, 
            rsd3_out=rsd3_out, rsd4_out=rsd4_out, rsd5_out=rsd5_out, exclude=exclude))
}

#Ions Extraction at a given sub dataset of rtgrp==T
#INPUT: xm1 is xM1() object, sgrp is sgroup() object, 
#       exclude is exclude5() object (NULL default), 
#       a = alpha value where 0<=a<=1, rtgrp = ret. time group to be analyzed.
#OUTPUT: list rtmed, mzmed, beta, chosen lambda min, d.f. 
ionx <- function (xm1,sgrp,rtgrp,exclude=NULL,a) { 
  if (!require("glmnet")) {
    install.packages("glmnet", dependencies = TRUE)
    library(glmnet)
  }
  
  Xs<-data.frame(xm1,check.names=FALSE)
  Sample<-Reduce(c,sgrp[[1]][2:length(sgrp[[1]])])
  conc.x= if(is.null(exclude)) subset(Xs, rtgroup == rtgrp, select = Sample) else
    subset(Xs[-(exclude$exclude),], rtgroup == rtgrp, select = Sample)
  count.p=nrow(conc.x)	; mz.p=rt.p=b=numeric(count.p) 
  #other information for dataframe
  mz.p= if(is.null(exclude)) subset(Xs, rtgroup == rtgrp, select = mzmed) else
    subset(Xs[-(exclude$exclude),], rtgroup == rtgrp, select = mzmed)
  rt.p= if(is.null(exclude)) subset(Xs, rtgroup == rtgrp, select = rtmed) else  
    subset(Xs[-(exclude$exclude),], rtgroup == rtgrp, select = rtmed)
  y=rep(c(1,2,3,4,5),c(12,12,12,12,12))
  ##choose lambda, use cross-validation to get lambda.min
  cv.fold=round(count.p*0.05)
  if (cv.fold<3) {
    cv.fold = 3}
  las.cv<-cv.glmnet(t(conc.x),y,nfolds=cv.fold,alpha=a)
  l=las.cv$lambda.min
  ##perform nest penalized regression to get beta
  las.mod<-glmnet(t(conc.x),y,lambda=l,alpha=a)
  b=as.matrix(las.mod$beta)
  df=las.mod$df
  result<-list(mzmed=mz.p,rtmed=rt.p,beta=b) 
  parameter.value<-list(cv.fold=cv.fold,lambda.min=l,df=df)
  return(inx=c(result,parameter.value))
}

#Ion-ion Relatonship Detection at a given sub dataset of rtgrp==T
#INPUT: xm2 is xM2() object, inx is ionx() object, it stands for ion (e)xtracted
#       sgrp= sgroup() object, exclude is exclude5() object (NULL default), 
#       a = alpha value where 0<=a<=1, rtgrp = ret. time group to be analyzed.
#OUTPUT: list rtmed, mzmed, beta, chosen lambda min, d.f. for each ion (e)xtracted.
ionr <- function(xm2,inx,sgrp,rtgrp,exclude=NULL,a){ 
  M1_inx<-which(inx$beta>0) #get rownames of non-zero coefficient from Model1
  M1_inx<-rownames(inx$beta)[M1_inx]
  inr=list()
  for (i in 1:length(M1_inx)){  
    y.index=as.numeric(M1_inx[i])
    Xs<-data.frame(xm2,check.names=FALSE)
    comb_samp_grp=Reduce(c,sgrp[[1]])
    conc.x = if(is.null(exclude)) subset(Xs[-y.index,], rtgroup == rtgrp, select = comb_samp_grp) else
      subset(Xs[-c(y.index,exclude$exclude),], rtgroup == rtgrp, select = comb_samp_grp)
    y=Xs[y.index,comb_samp_grp]
    count.p=nrow(conc.x)
    mz.p=rt.p=b=numeric(count.p)
    
    #other information for dataframe
    mz.p= if(is.null(exclude)) subset(Xs[-y.index], rtgroup == rtgrp, select = mzmed) else
      subset(Xs[-c(y.index,exclude$exclude),], rtgroup == rtgrp, select = mzmed)
    #subset(Xs[-y.index,], rtime == T, select = mzmed)
    rt.p= if(is.null(exclude)) subset(Xs[-y.index,], rtgroup == rtgrp, select = rtmed) else  
      subset(Xs[-c(y.index,exclude$exclude),], rtgroup == rtgrp, select = rtmed)
    #rt.p= subset(Xs[-y.index,], rtime == T, select = rtmed)
    
    cv.fold=round(count.p*0.05)
    if (cv.fold<3) {
      cv.fold=3}
    las.cv<-cv.glmnet(t(conc.x),t(y),nfolds=cv.fold,alpha=a)
    l=las.cv$lambda.min
    
    las.mod<-glmnet(t(conc.x),t(y),lambda=l,alpha=a)
    b=as.matrix(las.mod$beta)
    df=las.mod$df
    
    result<-list(mzmed=mz.p, rtmed=rt.p, beta=b)
    parameter.value<-list(index=y.index,cv.fold=cv.fold,lambda.min=l,df=df)
    inr[[i]]<-c(result,parameter.value)
  }
  
  return(inr)
}


#Summary and Visualization
#Ion Table
#INPUT: inx = ionx() object, inr = ionr() object
#OUTPUT: data frame for ion table
ionTable<-function(inx,inr){
  ion.table<-data.frame()
  sort_M1beta<-sort(inx$beta,decreasing=TRUE,index.return=TRUE)
  M1_inx<-which(sort_M1beta$x>0) #rowname inx (take only non-negative beta)
  M1_inx<-sort_M1beta$ix[M1_inx]
  M1_inx<-rownames(inx$beta)[M1_inx]
  uM1_inx<-which(inx$beta>0)
  uM1_inx<-rownames(inx$beta)[uM1_inx]
  for (i in 1:length(M1_inx)){
    k<-which(uM1_inx==M1_inx[3])
    round_M2beta<-round(inr[[k]]$beta,1)
    M2_inr<-which(round_M2beta!=0)  #which(inr[[k]]$beta!=0)
    M2_inr<-rownames(inr[[k]]$beta)[M2_inr] #get rownames of non-zero coefficients for the i-th inr
    M2_index<-as.numeric(M2_inr)
    M2_mz<-inr[[k]]$mzmed[M2_inr,]
    #M2_beta<-inr[[k]]$beta[M2_inr,]
    M2_rt<-inr[[k]]$rtmed[M2_inr,]
    M1_index<-rep(as.numeric(M1_inx[i]),length(M2_mz))
    M1_mz<-rep(inx$mzmed[M1_inx[i],],length(M2_mz))
    M1_rt<-rep(inx$rtmed[M1_inx[i],],length(M2_mz))
    #M1_beta<-rep(inx$beta[M1_inx[i],],length(M2_beta))
    gs<-data.frame(inx_index=M1_index,inx_mz=M1_mz,inx_rt=M1_rt,
                   inr_index=M2_inr,inr_mz=M2_mz,inr_rt=M2_rt) #inx_beta=M1_beta, ,inr_beta=M2_beta
    #gs<-gs[order(-gs$inr_beta),]#with descending order of beta
    ion.table<-rbind(ion.table,gs)
  }
  return(ion.table)
}

#Beta Table
#INPUT: inx = ionx() object, inr = ionr() object
#OUTPUT: dataframe for beta table
betaTable<-function(inx,inr){
  beta.table=data.frame()
  M1_mz=list()
  beta.table<-cbind(as.numeric(rownames(inx$mzmed)),inx$rtmed,inx$mzmed)
  M1_inx<-which(inx$beta>0) #get rownames of non-zero coefficients from Model1 (take only non-negative beta)
  M1_inx<-rownames(inx$beta)[M1_inx]
  #get beta from Model2
  for (i in 1:length(M1_inx)){
    M1_mz[i]<-round(inx$mzmed[M1_inx[i],],4)
    M2_beta<-round(inr[[i]]$beta,1) 
    #create dummy beta 
    dummy_beta<-matrix(0.9,1,1)
    rownames(dummy_beta)<-M1_inx[i]
    colnames(dummy_beta)<-colnames(M2_beta)
    #append dummy_beta to M2_beta
    M2_beta<-rbind(M2_beta,dummy_beta) 
    beta.table<-cbind(beta.table,M2_beta[match(rownames(beta.table), rownames(M2_beta))])
    names(beta.table)<-c("ori_index","rtmed","mzmed",M1_mz[1:length(M1_mz)])
  }
  return(beta.table)
}

#Ion Heatmap (unsorted)
ionHeatmap<-function(inx,inr){
  if (!require("gplots")) {
    install.packages("gplots", dependencies = TRUE)
    library(gplots)
  }
  if (!require("RColorBrewer")) {
    install.packages("RColorBrewer", dependencies = TRUE)
    library(RColorBrewer)
  }
  beta_table=betaTable(inx,inr)
  
  #colored polarity
  prevRow<-c(-Inf,beta_table$mzmed[-nrow(beta_table)])
  getRowId<-which(beta_table$mzmed<=prevRow) #marks the beginning of another ion-polarity
  getRowNames<-rownames(beta_table)[getRowId:nrow(beta_table)] #list of different polarity ions at time T
  col_ionrPolarity<-c(rep("#d95f02",nrow(beta_table)-length(getRowNames)), 
                      rep("#2c7bb6",length(getRowNames)))
  colCol_df<-data.frame(cbind(mzmed=round(beta_table$mzmed,4),colPolarity=col_ionrPolarity))
  clmn_names<-as.numeric(colnames(beta_table[,4:ncol(beta_table)]))
  colRow_df<-colCol_df[colCol_df$mzmed %in% clmn_names,] 
  
  if (nrow(beta_table)>100){
    print(paste("Data subset for RT=",(min(round(beta_table[,"rtmed"],1)))," - ",(max(round(beta_table[,"rtmed"],1))), " secs", "is too large. Please use manual plot for readability."))
  } else {
    
    cellnote<-t(as.matrix(beta_table[4:ncol(beta_table)]))
    cellnote[cellnote>0]<-1
    cellnote[cellnote<0]<--1
    unsort_hm<-heatmap.2(cellnote, na.col= "lightgray",trace="none", 
                         cellnote=cellnote, notecex=0.6, notecol="black",
                         col = brewer.pal(3, "YlGn"), margins=c(6,6), lhei=c(2.2,6), lwid=c(0.5,10), 
                         main=paste("Unsorted Beta coefficients for RT= ", 
                                    (min(round(beta_table[,"rtmed"],1)))," - ",
                                    (max(round(beta_table[,"rtmed"],1))), " secs", sep=""),
                         dendrogram="none", Rowv=FALSE, Colv=FALSE, 
                         labCol=round(beta_table[,"mzmed"],4), offsetRow=-0.8, offsetCol=-0.5,
                         colRow=as.character(colRow_df$colPolarity), cexRow=0.8,
                         colCol=as.character(colCol_df$colPolarity), cexCol=0.8, 
                         xlab=paste("All m/z at retention time group"),
                         ylab="Extracted m/z", cex.lab=1.5, 
                         key=FALSE) #for small subset (up to 100) 
   
  }
}

#Ion Heatmap (sorted)
ionHeatmap2<-function(inx,inr){
  if (!require("gplots")) {
    install.packages("gplots", dependencies = TRUE)
    library(gplots)
  }
  if (!require("RColorBrewer")) {
    install.packages("RColorBrewer", dependencies = TRUE)
    library(RColorBrewer)
  }
  beta_table=betaTable(inx,inr)
  
  #colored polarity
  prevRow<-c(-Inf,beta_table$mzmed[-nrow(beta_table)])
  getRowId<-which(beta_table$mzmed<=prevRow) #This row in beta_table marks the beginning of another ion-polarity
  getRowNames<-rownames(beta_table)[getRowId:nrow(beta_table)] #list of different polarity ions at time T
  col_ionrPolarity<-c(rep("#d95f02",nrow(beta_table)-length(getRowNames)), 
                      rep("#2c7bb6",length(getRowNames)))
  colCol_df<-data.frame(cbind(mzmed=round(beta_table$mzmed,4),colPolarity=col_ionrPolarity))
  clmn_names<-as.numeric(colnames(beta_table[,4:ncol(beta_table)]))
  colRow_df<-colCol_df[colCol_df$mzmed %in% clmn_names,] 
  
  #Sort beta in the beta_table, then create data frame
  ionxdf<-colSums(na.omit(beta_table[,4:ncol(beta_table)] != 0)) 
  ionrdf<-rowSums(beta_table[,4:ncol(beta_table)] != 0) 
  sortCol<-names(sort(ionxdf,decreasing=FALSE))
  sortRow<-names(sort(ionrdf,decreasing=FALSE))
  betaTable.df1<-beta_table[,c("rtmed","mzmed",sortCol)] #sorted only Col
  betaTable.df2<-beta_table[match(sortRow,rownames(beta_table)),] #sorted only Row
  betaTable.df3<-betaTable.df2[,c("rtmed","mzmed",sortCol)] #sorted both Row and Col
  sortcolRow_df<-colRow_df[match(colnames(betaTable.df3[3:ncol(betaTable.df3)]),colRow_df[,1]),]
  sortcolCol_df<-colCol_df[match(round(betaTable.df3$mzmed,4),colCol_df[,1]),]
  
  if (nrow(beta_table)>100){
    print(paste("Data subset for RT=",(min(round(beta_table[,"rtmed"],1)))," - ",(max(round(beta_table[,"rtmed"],1))), " secs", "is too large. Please use manual plot for readability."))
  } else {
    
    cellnote<-t(as.matrix(betaTable.df3[3:ncol(betaTable.df3)]))
    cellnote[cellnote>0]<-1
    cellnote[cellnote<0]<--1
    sort_hm<-heatmap.2(cellnote, trace="none", 
                       cellnote=cellnote, notecex=0.6, notecol="black",
                       col = brewer.pal(3, "YlGn"), margins=c(6,6), lhei=c(2.2,6), lwid=c(0.5,10), 
                       main=paste("Sorted Beta coefficients for RT= ", 
                                  (min(round(beta_table[,"rtmed"],1)))," - ",
                                  (max(round(beta_table[,"rtmed"],1))), " secs", sep=""),
                       dendrogram="none", Rowv=FALSE, Colv=FALSE, 
                       colRow=as.character(sortcolRow_df$colPolarity), cexRow=0.8,
                       colCol=as.character(sortcolCol_df$colPolarity), cexCol=0.8, 
                       labCol=round(betaTable.df3[,"mzmed"],4), 
                       xlab=paste("All m/z at retention time group"),
                       ylab="Extracted m/z", cex.lab=1.5, 
                       key=FALSE, offsetRow = -0.8, offsetCol = -0.5) #for small subset (up to 100)
  }
}

#univariate model 1 plot concentration vs. intensity
#INPUT: xm1 = xM1() object, inx = ionx() object and sgrp = sgroup() object
ionxPlot<-function(xm1,inx,sgrp){
  kgroup<-length(sgrp[[1]])-1
  nsamp<-Reduce(c,sgrp[[2]][2:length(sgrp[[2]])])
  x.axis<-rep(seq(1,kgroup),c(nsamp))
  Sample<-Reduce(c,sgrp[[1]][2:length(sgrp[[1]])])
  M1_inx<-which(inx$beta!=0) #get rownames of non-zero coefficients from Model1
  M1_inx<-rownames(inx$beta)[M1_inx]
  for (i in 1:length(M1_inx)){ 
    y.axis<-xm1[as.numeric(M1_inx[i]),Sample]
    mz.info<-round(xm1[as.numeric(M1_inx[i]),"mzmed"],4)
    rt.info<-round(xm1[as.numeric(M1_inx[i]),"rtmed"],1)
    plot(x.axis,y.axis,main=paste("m/z=",mz.info," RT =",rt.info,sep=""),xlab="Samples",
         ylab="log2_xM1 ion intensity",pch=16, cex.lab=1.2, cex.main=1.5)
    lines(x.axis,y.axis,lwd=2)
    readline(prompt="Press [enter] to continue")
  }
}

#univariate model 2 plot extracted ion intensity vs. another ion intensity
#INPUT: xm2 = xM2() object, inr = ionr() object, sgrp = sgroup() object 
#       and frac = proportion of top beta to plot
ionrPlot<-function(xm2,inr,sgrp,frac=0.2){
  nionx<-length(inr) 
  comb_samp_grp=Reduce(c,sgrp[[1]]) 
  for (i in 1:nionx){
    M2_inr<-which(inr[[i]]$beta!=0)
    M2_inr<-rownames(inr[[i]]$beta)[M2_inr]
    x.axis<-xm2[inr[[i]]$index,comb_samp_grp] 
    x.info<-round(xm2[inr[[i]]$index,"mzmed"],4)
    nplot=round(frac*(length(M2_inr)),0) #fraction of ions that users expected to be plotted (from hi beta to lo beta)
    if (nplot==0){
      print(paste("No plot for m/z=", x.info)) ; next
    }
    sort_beta<-sort(inr[[i]]$beta,decreasing=TRUE,index.return=TRUE) #to get top beta from large to small
    #top.beta2<-tail(sort(corr_ions[[8]]$beta),nplot) #to get top beta from small to large 
    #what if there are strong negativelly correlated ions? How to get head and tail?
    top.beta.ix<-head(rownames(inr[[i]]$beta)[sort_beta$ix],nplot)
    for (j in 1:length(top.beta.ix)){ #length(top.beta.ix)=nplot
      y.axis<-xm2[as.numeric(top.beta.ix[j]),comb_samp_grp]
      y.info<-round(xm2[as.numeric(top.beta.ix[j]),"mzmed"],4)
      rt.info<-round(xm2[as.numeric(top.beta.ix[j]),"rtmed"],1)
      plot(x.axis,y.axis,main=paste("Ion-Ion Relationship Intensity at RT =", rt.info),
           xlab=paste("log2_xM2 m/z =", x.info,"intensity"),ylab=paste("log2_xM2 m/z =", y.info, "intensity"),
           pch=16, cex.lab=1.2, cex.main=1.5) 
      readline(prompt="Press [enter] to continue")
    }
  }
}

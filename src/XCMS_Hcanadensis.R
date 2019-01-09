###Data Processing for Hydrastis canadensis Extraction with different conc.
###Dataset: 170918_BR16006_62 (Real data)

#Installation Note:
#To install XCMS and CAMERA go to
#http://www.bioconductor.org/packages/3.5/bioc/html/xcms.html
#http://www.bioconductor.org/packages/release/bioc/html/CAMERA.html

library(xcms)
library(CAMERA)
citation(package="xcms")
citation(package="CAMERA")



setwd("C:/Xcalibur/data/Bety/170918_BR16006_62")
#*Not used
#*data1<-list.files("./positive",recursive=T,full=T) #negative
#*data1 #all files are included in the allignment 

#data2_pos/neg indicates only Samples and sample blanks are included
data2_pos<-list.files("./positive/positive_sample",recursive=T,full=T) 
data2_neg<-list.files("./negative/negative_sample",recursive=T,full=T) 

#data3_pos/neg indicates only QC and QC blanks
data3_pos<-list.files("./positive/positive_QC",recursive=T,full=T) 
data3_pos<-list.files("./negative/negative_QC",recursive=T,full=T) 


pos_xset<-xcmsSet(data2_pos, method='centWave',ppm=3, peakwidth=c(2,20),
              polarity="positive", snthresh=3, integrate=2, mzdiff=-1.0000)
#change data2_pos to data3_pos for QC positive mode
#change snthresh=10 for sn10 dataset
  #integrate=2, the descent is done on the real data, 
    #more accurate but prone to noise
  #mzdiff is min difference in m/z for peaks with overlapping r.t.,
    #at a min. of -1 H-atom or -1 isotopic mass; neg value indicates overlap
pos_xset<-group(pos_xset,method="density",bw=3,mzwid=0.015,minfrac=0.5)
pos_xset1<-retcor(pos_xset,family="s",plottype="m",missing=1,smooth="linear")
pos_xset2<-group(pos_xset1,method="density",bw=1.5,mzwid=0.015,minfrac=0.2)
#minfrac at 0.2 will give the max of 12 out of 60 samples for a peak to become valid peaks
pos_xset3<-fillPeaks(pos_xset2, method="chrom",expand.mz=1, expand.rt=1)
pos_xsann<-annotateDiffreport(pos_xset3,sigma=3,perfwhm=0.5,
                          intval=c("intb","maxo","into"), graphMethod="hcs", 
                          calcCiS=TRUE, calcCaS=FALSE, ppm=3, 
                          cor_eic_th=0.75, polarity="positive", 
                          quick=FALSE, sortpval=FALSE)

neg_xset<-xcmsSet(data2_neg, method='centWave',ppm=3, peakwidth=c(2,20),
              polarity="negative", snthresh=3, integrate=2, mzdiff=-1.0000)
#change data2_pos to data3_pos for QC negative mode
#change snthresh=10 for sn10 dataset
neg_xset<-group(neg_xset,method="density",bw=3,mzwid=0.015,minfrac=0.5)
neg_xset1<-retcor(neg_xset,family="s",plottype="m",missing=1,smooth="linear")
neg_xset2<-group(neg_xset1,method="density",bw=1.5,mzwid=0.015,minfrac=0.2)
neg_xset3<-fillPeaks(neg_xset2, method="chrom",expand.mz=1, expand.rt=1)
neg_xsann<-annotateDiffreport(neg_xset3,sigma=3,perfwhm=0.5,
                          intval=c("intb","maxo","into"), graphMethod="hcs", 
                          calcCiS=TRUE, calcCaS=FALSE, ppm=3, 
                          cor_eic_th=0.75, polarity="negative", 
                          quick=FALSE, sortpval=FALSE)

#Create .csv report for XCMS object (Optional)
xset4_data2<-write.csv(c(pos_xsann, neg_xsann),file="filename.csv") #change filename accordingly
xset4_data3<-write.csv(c(pos_xsann,neg_xsann),file="filename.csv") #change filename accordingly
#Other way to create .csv report for XCMS object (Alternative)
xset4<-diffreport(xset3, filebase="filename") #create report with stats but without annotation from CAMERA
xset4<-peakTable(xset3, filebase="filename") #create only peaktable, no stats and no annotation


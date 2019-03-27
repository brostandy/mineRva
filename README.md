# mineRva
metabolite ions extRaction and visualization (mineRva) is a statistical analysis model to extract ions from LC-MS1 experiment of a complex mixture analysis. The tool is written in R environment. 

## Introduction
Metabolite ions from LC-MS1 experiment usually is the initial step to model or generate hypothesis for metabolomics study. However,
MS1 data is usually too limited for the purpose of metabolite identification. Typically, MS2 or MSn data are acquired to enable compound
library search. MS2 or MSn data are acquired in LC-MS/MS experiment either by Data Dependent Analysis or some research had also
introduced Data Independent Analysis, where All Ion Fragmentation (AIF) process are aligned in the experiment. These MS/MS acquisition 
lead to higher volume of data to be analyzed. Therefore, in our attempt to cut analytical workload on mass spectrometric data, we utilized 
only MS1 data to better direct analysts in further MS/MS analysis without overwhelming analysts with processing MS2 data.

mineRva is a tool created in R-languange to extract metabolite ions with a particular experimental design for complex mixture analysis via LC-ESI-MS. This tool enables analysts to extract putative ions that belongs to the biological sample, thus directing analysts for targeted MS/MS experiment.

## Instructions
1. Download the minervaTools.R source code by cloning your own copy of the source repository with the commands:

   $> git clone https://github.com/brostandy/mineRva.git<br/>
   $> cd mineRva

2. Usage, in R console, type:
   > source("minervaTools.R")

## mineRva Requirement Credits
mineRva tool consists of other R packages:

R Core Team (2018). R: A language and environment for statistical computing. R Foundation for
  Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.
  
Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010). Regularization Paths for Generalized
  Linear Models via Coordinate Descent. Journal of Statistical Software, 33(1), 1-22. URL
  http://www.jstatsoft.org/v33/i01/.
  
Gregory R. Warnes, Ben Bolker, Lodewijk Bonebakker, Robert Gentleman, Wolfgang Huber Andy Liaw,
  Thomas Lumley, Martin Maechler, Arni Magnusson, Steffen Moeller, Marc Schwartz and Bill
  Venables (2016). gplots: Various R Programming Tools for Plotting Data. R package version
  3.0.1. https://CRAN.R-project.org/package=gplots
  
Erich Neuwirth (2014). RColorBrewer: ColorBrewer Palettes. R package version 1.1-2.
  https://CRAN.R-project.org/package=RColorBrewer
  



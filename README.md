# mineRva
metabolite ions extRaction and visualization (mineRva) is to extract ions from LC-MS1 experiment of a complex mixture analysis

## Introduction
Metabolite ions from LC-MS1 experiment usually is the initial step to model or generate hypothesis for metabolomics study. However,
MS1 data is usually too limited for the purpose of metabolite identification. Typically, MS2 or MSn data are acquired to enable compound
library search. MS2 or MSn data are acquired in LC-MS/MS experiment either by Data Dependent Analysis or some research had also
introduced Data Independent Analysis, where All Ion Fragmentation (AIF) process are aligned in the experiment. These MS/MS acquisition 
lead to higher volume of data to be analyzed. Therefore, in our attempt to cut analytical workload on mass spectrometric data, we utilized 
only MS1 data to better direct analysts in further MS/MS analysis without overwhelming analysts with processing MS2 data.

mineRva is a tool created in R-languange to extract metabolite ions with a particular experimental design for complex mixture analysis with
LC-ESI-MS. This tool enables analysts to extract putative ions that belongs to the biological sample, thus directing analysts for targeted
MS/MS experiment.


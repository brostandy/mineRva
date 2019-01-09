mineRva Functions and Workflow
================
Bety Rostandy
January 8, 2019

### What is mineRva?

mineRva is an R package to extract metabolite ions for a complex mixture of LC-MS analysis after preprocessing steps (i.e peaks filling, peaks allignment, retention time correction, and so on). The analytical statistics model is build for LC-MS complex (i.e. biological) mixture dataset that must contain retention time, mass-to-charge ratio, and measured ions' intensity from samples. This package is build to analyze specific design of experiment in measuring biomolecules via LC-MS system.

### Metabolite Ions Extraction and Corelation Functions

Users feed in preprocessed LC-MS data, where ions' peak had been deconvoluted to a table. Then, to enable extraction of dataset from data, we are going to define sample groups with `sgrp()`. Users give groupnames based on the column names that define the sample groups of preprocessed data. The first group is always the negative control or blanks.

After that, the dataset is extracted and zero-entries are treated with `dataSet()`.

This dataset was extracted and zero-entries is treated the same for all model during this step, however, each model is rescaled and transformed accordingly. For Model 1, the regression use concentration of the extract as response variable, and ion intensity is normalized to be the ratio of samples vs. blanks. Model 1 follows such equation:

$$log2\\frac{I\_{i,M\_p}^k}{\\frac{1}{27}\\sum\_{i=1}^{27}I^0\_{i,M\_p}}$$
 where, *k* = group of extract concentration from 1 to 5 ;

*i* = number of sample in the group of extract concentration (each group has 12 samples, except blank has 27 samples) ;

*M*<sub>*p*</sub> = ion observed at a time == T, where T is the time round up to nearest second(s) ;

*I* = ion intensity observed.

Notice that *M*<sub>*p*</sub> is ion observed at a time == T; this means that we subset out the data into different time point (rounding up to second).

Therefore, next step is to rescale and transform the dataset for model 1 with `xM1()`

Then, users can choose to perform RSD cutoff to determine outlier (only for Model 1, as Model 2 uses Model 1 result). If users choose to perform RSD cutoff an exclude list will be returned. This exclude list follows certain conditions that is set for five sample groups experiment (i.e. six *groupnames* with blank).The outlier for RSD cutoff is determined at such: *X* &gt; *Q*3 + 1.5 \* *I**Q**R* or *X* &lt; *Q*1 − 1.5 \* *I**Q**R*. In order to get RSD, mean and standard deviation calculations must be performed. Therefore, `msr()` performs calculations for mean(*m*), standard deviation (*s*) and relative standard deviatio (*r*).

The exclude list is returned after running the function `exclude5()`. (Note: as of current development, ions exclusion can only be done if users have five sample groups and one negative control.)

Finally, Model 1 can be performed (with or without outlier removal) with `ionx()`. Model 1 has the objective of extracting significant ions from LC-MS data, therefore the function `ionx()` stands for ion (e)xtraction. $More of the statistical methodology for the model will be further discussed in publication.

After fitting model 1, where metabolite ions are extracted at a round-up retention time (*T*) given by the user, these extracted ions of model 1 then fitted through model 2 to find its correlated ions. Ion intensity is normalized for model 2 before correlation. Each ion follow model 2 equation below:

$$log2\\frac{I\_{i,M\_p}^{k\_2}}{median\_{p}}$$
 where, *k*<sub>2</sub> = groups for model 2 from 1 to 6, includes negative control (i.e. blank) ;

*i* = number of sample in the group of extract concentration (each group has 12 samples, except blank has 27 samples) ;

*M*<sub>*p*</sub> = ion observed at a time == T, where T is the time round up to nearest second(s) ;

*I* = ion intensity observed.

*m**e**d**i**a**n*<sub>*p*</sub> = median intensity of ion observed from 1 to *p*.

Hence, the function `xM2()` performs data scaling and transformation for model 2.

Model 2 function, `ionr()`, is to perform correlation of ions that were extracted from Model 1. For each selected ion from Model 1, we against it with all other ions at the given retention time. Ions that are correlated to the selected ions were pulled out.

User can call model 1 and model 2 functions iteratively as below:

``` r
#set working directory
setwd("C:/Users/bty/OneDrive - UNCG/Projects/Goldenseal/170927_Hcanadensis with diff conc")

#data<-read.csv("Data/180517_Sample_pos_neg_sn3.csv",check.names="FALSE")
data<-read.csv("Data/180517_Sample_pos_neg_sn10.csv",check.names="FALSE")

sample_group=sgroup(data,c("MeOH","03125","0625","125","250","500"))
dataset<-dataSet(data,sample_group,rt="rtmed",mz="mzmed")
x_method1<-xM1(dataset,sample_group)
msr<-msr(x_method1,sample_group)
h<-exclude5(x_method1,msr,T=155)
extract_ions<-ionx(xm=x_method1,sample_group,T=155,a=0.8) #exclude=h
x_method2<-xM2(dataset,sample_group)
corr_ions<-ionr(xm=x_method2,extract_ions,sample_group,T=155,a=0.8)
```

After performing model 1 and model 2, a summary of an *i*-th ion extracted and ion correlated to it can be produced with `ionTable()`.

We then be able to see model 1 and model 2 results in a dataframe where:

``` r
ion_table<-ionTable(extract_ions,corr_ions)
ion_table
##    inx_index    inx_mz   inx_rt inr_index    inr_mz   inr_rt
## 1        637  95.04397 155.2691      1306 102.03395 155.4705
## 2        637  95.04397 155.2691      4858 158.96146 154.9756
## 3        637  95.04397 155.2691      6877 208.03949 155.4236
## 4        637  95.04397 155.2691     13081 785.28284 154.8956
## 5        862  98.01370 155.1858      1306 102.03395 155.4705
## 6        862  98.01370 155.1858      5759 177.00717 155.3517
## 7        862  98.01370 155.1858     14986 190.92879 154.5651
## 8       1306 102.03395 155.4705       637  95.04397 155.2691
## 9       1306 102.03395 155.4705       862  98.01370 155.1858
## 10      1306 102.03395 155.4705      1627 109.02250 155.0114
## 11      1306 102.03395 155.4705      2138 112.01802 155.2477
## 12      1306 102.03395 155.4705      6877 208.03949 155.4236
## 13      1306 102.03395 155.4705     10795 401.14241 155.4359
## 14      1306 102.03395 155.4705     13081 785.28284 154.8956
## 15      1306 102.03395 155.4705     13108 800.27447 155.4060
## 16      4026 143.99697 155.2523     13108 800.27447 155.4060
## 17      4725 156.99083 155.4705      1306 102.03395 155.4705
## 18      4725 156.99083 155.4705      2138 112.01802 155.2477
## 19      4725 156.99083 155.4705      3892 143.03974 155.4589
## 20      4725 156.99083 155.4705      5759 177.00717 155.3517
## 21      4725 156.99083 155.4705     10066 358.16413 154.8946
## 22      4725 156.99083 155.4705     10774 400.13906 155.4998
## 23      4725 156.99083 155.4705     10795 401.14241 155.4359
## 24      4858 158.96146 154.9756       637  95.04397 155.2691
## 25      4858 158.96146 154.9756      6877 208.03949 155.4236
## 26      5759 177.00717 155.3517       862  98.01370 155.1858
## 27      5759 177.00717 155.3517      1627 109.02250 155.0114
## 28      5759 177.00717 155.3517      3497 138.99543 154.9399
## 29      5759 177.00717 155.3517      4026 143.99697 155.2523
## 30      5759 177.00717 155.3517      4725 156.99083 155.4705
## 31      5759 177.00717 155.3517      4858 158.96146 154.9756
## 32      5759 177.00717 155.3517      5853 180.98994 154.6402
## 33      5759 177.00717 155.3517      6877 208.03949 155.4236
## 34      5759 177.00717 155.3517     10036 356.18531 155.4821
## 35      5759 177.00717 155.3517     10066 358.16413 154.8946
## 36      5759 177.00717 155.3517     10774 400.13906 155.4998
## 37      5759 177.00717 155.3517     14986 190.92879 154.5651
## 38      6877 208.03949 155.4236       637  95.04397 155.2691
## 39      6877 208.03949 155.4236      1306 102.03395 155.4705
## 40      6877 208.03949 155.4236      2138 112.01802 155.2477
## 41      6877 208.03949 155.4236      4858 158.96146 154.9756
## 42      6877 208.03949 155.4236      5113 161.96972 155.3870
## 43      6877 208.03949 155.4236      5759 177.00717 155.3517
## 44      6877 208.03949 155.4236     13108 800.27447 155.4060
## 45     10066 358.16413 154.8946       637  95.04397 155.2691
## 46     10066 358.16413 154.8946      3497 138.99543 154.9399
## 47     10066 358.16413 154.8946      3892 143.03974 155.4589
## 48     10066 358.16413 154.8946      5759 177.00717 155.3517
## 49     10066 358.16413 154.8946     10774 400.13906 155.4998
## 50     10066 358.16413 154.8946     10795 401.14241 155.4359
## 51     10066 358.16413 154.8946     13081 785.28284 154.8956
## 52     10066 358.16413 154.8946     13106 799.27063 155.4821
## 53     10066 358.16413 154.8946     14986 190.92879 154.5651
## 54     10774 400.13906 155.4998      1627 109.02250 155.0114
## 55     10774 400.13906 155.4998      2138 112.01802 155.2477
## 56     10774 400.13906 155.4998      3497 138.99543 154.9399
## 57     10774 400.13906 155.4998      3892 143.03974 155.4589
## 58     10774 400.13906 155.4998      4026 143.99697 155.2523
## 59     10774 400.13906 155.4998      4858 158.96146 154.9756
## 60     10774 400.13906 155.4998      5113 161.96972 155.3870
## 61     10774 400.13906 155.4998      5759 177.00717 155.3517
## 62     10774 400.13906 155.4998      6877 208.03949 155.4236
## 63     10774 400.13906 155.4998     10066 358.16413 154.8946
## 64     10774 400.13906 155.4998     10795 401.14241 155.4359
## 65     10774 400.13906 155.4998     13081 785.28284 154.8956
## 66     10774 400.13906 155.4998     13106 799.27063 155.4821
## 67     10774 400.13906 155.4998     13108 800.27447 155.4060
## 68     10795 401.14241 155.4359      1306 102.03395 155.4705
## 69     10795 401.14241 155.4359      1627 109.02250 155.0114
## 70     10795 401.14241 155.4359      2138 112.01802 155.2477
## 71     10795 401.14241 155.4359      3497 138.99543 154.9399
## 72     10795 401.14241 155.4359      3892 143.03974 155.4589
## 73     10795 401.14241 155.4359      4026 143.99697 155.2523
## 74     10795 401.14241 155.4359      4858 158.96146 154.9756
## 75     10795 401.14241 155.4359      6877 208.03949 155.4236
## 76     10795 401.14241 155.4359     10066 358.16413 154.8946
## 77     10795 401.14241 155.4359     10774 400.13906 155.4998
## 78     10795 401.14241 155.4359     14986 190.92879 154.5651
## 79     13108 800.27447 155.4060       862  98.01370 155.1858
## 80     13108 800.27447 155.4060      1306 102.03395 155.4705
## 81     13108 800.27447 155.4060      1627 109.02250 155.0114
## 82     13108 800.27447 155.4060      2138 112.01802 155.2477
## 83     13108 800.27447 155.4060      3497 138.99543 154.9399
## 84     13108 800.27447 155.4060      3892 143.03974 155.4589
## 85     13108 800.27447 155.4060      4026 143.99697 155.2523
## 86     13108 800.27447 155.4060      6877 208.03949 155.4236
## 87     13108 800.27447 155.4060     10036 356.18531 155.4821
## 88     13108 800.27447 155.4060     10774 400.13906 155.4998
## 89     13108 800.27447 155.4060     10795 401.14241 155.4359
## 90     13108 800.27447 155.4060     13081 785.28284 154.8956
## 91     13108 800.27447 155.4060     13106 799.27063 155.4821
## 92     13108 800.27447 155.4060     14986 190.92879 154.5651
## 93     14986 190.92879 154.5651       862  98.01370 155.1858
## 94     14986 190.92879 154.5651      4858 158.96146 154.9756
## 95     14986 190.92879 154.5651      5113 161.96972 155.3870
## 96     14986 190.92879 154.5651     10036 356.18531 155.4821
## 97     14986 190.92879 154.5651     10066 358.16413 154.8946
```

Another type of dataframe can be generated for beta coefficients in model 2, we called this function `betaTable()`. This table collects beta coefficients of all correlated ions to every extracted extracted ions.

We can see the beta table from here. Since for illustration we used T at 155 second, below is the beta table for T==155.

``` r
beta_table<-betaTable(extract_ions,corr_ions)
beta_table
##          rtmed     mzmed 95.044 98.0137 102.034 143.0397 143.997 156.9908
## 637   155.2691  95.04397   0.99    0.00    0.05     0.00    0.00     0.00
## 862   155.1858  98.01370   0.00    0.99   -0.02     0.00    0.00     0.00
## 1306  155.4705 102.03395   0.01   -0.12    0.99     0.00    0.00     0.01
## 1627  155.0114 109.02250   0.00    0.00    0.08     0.00    0.00     0.00
## 2138  155.2477 112.01802   0.00    0.00    0.15     0.00    0.00     0.16
## 3497  154.9399 138.99543   0.00    0.00    0.00     0.00    0.00     0.00
## 3892  155.4589 143.03974   0.00    0.00    0.00     0.99    0.00     0.03
## 4026  155.2523 143.99697   0.00    0.00    0.00     0.00    0.99     0.00
## 4725  155.4705 156.99083   0.00    0.00    0.00     0.00    0.00     0.99
## 4858  154.9756 158.96146  -0.06    0.00    0.00     0.00    0.00     0.00
## 5113  155.3870 161.96972   0.00    0.00    0.00     0.00    0.00     0.00
## 5759  155.3517 177.00717   0.00    0.09    0.00     0.00    0.00     0.07
## 5853  154.6402 180.98994   0.00    0.00    0.00     0.00    0.00     0.00
## 6877  155.4236 208.03949   0.22    0.00    0.11     0.00    0.00     0.00
## 10036 155.4821 356.18531   0.00    0.00    0.00     0.00    0.00     0.00
## 10066 154.8946 358.16413   0.00    0.00    0.00     0.00    0.00     0.00
## 10774 155.4998 400.13906   0.00    0.00    0.00     0.00    0.00    -0.02
## 10795 155.4359 401.14241   0.00    0.00   -0.04     0.00    0.00     0.00
## 13081 154.8956 785.28284   0.00    0.00   -0.01     0.00    0.00     0.00
## 13106 155.4821 799.27063   0.00    0.00    0.00     0.00    0.00     0.00
## 13108 155.4060 800.27447   0.00    0.00    0.00     0.00    0.00     0.00
## 14986 154.5651 190.92879   0.00    0.06    0.00     0.00    0.00     0.00
##       158.9615 161.9697 177.0072 208.0395 358.1641 400.1391 401.1424
## 637      -0.09     0.00     0.00     0.69     0.30     0.00     0.00
## 862       0.00     0.00     0.13     0.00     0.00     0.00     0.00
## 1306      0.00     0.00     0.00     0.12     0.00     0.00    -0.02
## 1627      0.00     0.00     0.21     0.00     0.00     0.15    -0.16
## 2138      0.00     0.00     0.00     0.03     0.00     0.00    -0.03
## 3497      0.00     0.00     0.02     0.00    -0.02     0.04    -0.03
## 3892      0.00     0.00     0.00     0.00     0.02    -0.12     0.08
## 4026      0.00     0.00     0.27     0.00     0.00    -0.06     0.03
## 4725      0.00     0.00     0.17     0.00     0.00     0.00     0.00
## 4858      0.99     0.00     0.03    -0.12     0.00    -0.11     0.10
## 5113      0.00     0.99     0.00    -0.05     0.00    -0.04     0.00
## 5759      0.00     0.00     0.99    -0.07    -0.05    -0.05     0.00
## 5853      0.00     0.00     0.02     0.00     0.00     0.00     0.00
## 6877     -0.04     0.00    -0.13     0.99     0.00     0.09    -0.12
## 10036     0.00     0.00     0.15     0.00     0.00     0.00     0.00
## 10066     0.00     0.00    -0.01     0.00     0.99     0.00     0.00
## 10774     0.00     0.00    -0.01     0.00     0.30     0.99     1.03
## 10795     0.00     0.00     0.00     0.00     0.17     0.89     0.99
## 13081     0.00     0.00     0.00     0.00     0.08     0.01     0.00
## 13106     0.00     0.00     0.00     0.00     0.08     0.01     0.00
## 13108     0.00     0.00     0.00    -0.01     0.00     0.00     0.00
## 14986     0.00     0.00    -0.06     0.00    -0.23     0.00     0.00
##       800.2745 190.9288
## 637       0.00     0.00
## 862      -0.12     0.02
## 1306     -0.22     0.00
## 1627      1.28     0.00
## 2138     -0.59     0.00
## 3497      0.34     0.00
## 3892     -0.77     0.00
## 4026     -0.47     0.00
## 4725      0.00     0.00
## 4858      0.00    -0.04
## 5113      0.00    -0.06
## 5759      0.00     0.00
## 5853      0.00     0.00
## 6877     -0.80     0.00
## 10036     1.28     0.08
## 10066     0.00    -0.02
## 10774     0.69     0.00
## 10795     0.23     0.00
## 13081     0.20     0.00
## 13106     0.13     0.00
## 13108     0.99     0.00
## 14986     0.03     0.99
```

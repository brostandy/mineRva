mineRva Functions and Workflow
================
Bety Rostandy
January 8, 2019

### What is mineRva?

mineRva is an R package to extract metabolite ions for a complex mixture of LC-MS analysis after preprocessing steps (i.e peaks filling, peaks allignment, retention time correction, and so on). The analytical statistics model is build for LC-MS complex (i.e. biological) mixture dataset that must contain retention time, mass-to-charge ratio, and measured ions' intensity from samples. This package is build to analyze specific design of experiment in measuring biomolecules via LC-MS system.

### Metabolite Ions Extraction and Corelation Functions

Users feed in preprocessed LC-MS data, where ions' peak had been deconvoluted to a table. Then, to enable extraction of dataset from data, we are going to define sample groups with `sgrp()`. Users give groupnames based on the column names that define the sample groups of preprocessed data. The first group is always the negative control or blanks.

After that, the dataset is extracted and zero-entries are treated with `dataSet()`. This dataset is extracted and its zero-entries are treated the same for all model, however, each model was rescaled and transformed accordingly.

Therefore to rescale and transform the dataset for model 1, we use `xM1()`.

Then, users can choose to perform RSD cutoff to determine outlier (only for Model 1, as Model 2 uses Model 1 result). If users choose to perform RSD cutoff an exclude list will be returned. This exclude list follows certain conditions that is set for five sample groups experiment (i.e. six *groupnames* with blank).The outlier for RSD cutoff is determined at such: *X* &gt; *Q*3 + 1.5 \* *I**Q**R* or *X* &lt; *Q*1 − 1.5 \* *I**Q**R*. In order to get RSD, mean and standard deviation calculations must be performed. Therefore, `msr()` performs calculations for mean(*m*), standard deviation (*s*) and relative standard deviatio (*r*).

The exclude list is returned after running the function `exclude5()`. (Note: as of current development, ions exclusion can only be done if users have five sample groups and one negative control.)

Finally, Model 1 can be performed (with or without outlier removal) with `ionx()`. Model 1 has the objective of extracting significant ions from LC-MS data, therefore the function `ionx()` stands for ion (e)xtraction.

After fitting model 1, where metabolite ions are extracted at a round-up retention time (*T*) given by the user, these extracted ions of model 1 then fitted through model 2 to find its correlated ions. Ion intensity is normalized for model 2 before correlation with xM2()\`.

Model 2 function, `ionr()`, is to perform correlation of ions that were extracted from Model 1. For each selected ion from Model 1, we against it with all other ions at the given retention time. Ions that are correlated to the selected ions are pulled out.

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
extract_ions<-ionx(xm=x_method1,sample_group,exclude=h,T=155,a=0.8) #
x_method2<-xM2(dataset,sample_group)
corr_ions<-ionr(xm=x_method2,extract_ions,sample_group,T=155,a=0.8)
```

After performing model 1 and model 2, an ion table of an *i*-th ion extracted and ion correlated to it can be produced with `ionTable()`.

We then are able to see model 1 and model 2 results in a dataframe where:

``` r
ion_table<-ionTable(extract_ions,corr_ions)
ion_table
##    inx_index   inx_mz   inx_rt inr_index    inr_mz   inr_rt
## 1       3497 138.9954 154.9399      1627 109.02250 155.0114
## 2       3497 138.9954 154.9399      5113 161.96972 155.3870
## 3       3497 138.9954 154.9399      5853 180.98994 154.6402
## 4       3497 138.9954 154.9399     10036 356.18531 155.4821
## 5       5853 180.9899 154.6402      1627 109.02250 155.0114
## 6       5853 180.9899 154.6402      3497 138.99543 154.9399
## 7      10036 356.1853 155.4821       637  95.04397 155.2691
## 8      10036 356.1853 155.4821      2138 112.01802 155.2477
## 9      10036 356.1853 155.4821      3497 138.99543 154.9399
## 10     10036 356.1853 155.4821      3892 143.03974 155.4589
## 11     10036 356.1853 155.4821      4026 143.99697 155.2523
## 12     10036 356.1853 155.4821      4725 156.99083 155.4705
## 13     10036 356.1853 155.4821      4858 158.96146 154.9756
## 14     10036 356.1853 155.4821      5759 177.00717 155.3517
## 15     10036 356.1853 155.4821      5853 180.98994 154.6402
## 16     10036 356.1853 155.4821      6877 208.03949 155.4236
## 17     10036 356.1853 155.4821     10066 358.16413 154.8946
## 18     10036 356.1853 155.4821     10774 400.13906 155.4998
## 19     10036 356.1853 155.4821     10795 401.14241 155.4359
## 20     10036 356.1853 155.4821     13081 785.28284 154.8956
## 21     10036 356.1853 155.4821     13106 799.27063 155.4821
## 22     10036 356.1853 155.4821     13108 800.27447 155.4060
## 23     10036 356.1853 155.4821     14986 190.92879 154.5651
## 24     10066 358.1641 154.8946     10774 400.13906 155.4998
## 25     10066 358.1641 154.8946     10795 401.14241 155.4359
## 26     10066 358.1641 154.8946     13081 785.28284 154.8956
## 27     10066 358.1641 154.8946     13106 799.27063 155.4821
## 28     10774 400.1391 155.4998      1627 109.02250 155.0114
## 29     10774 400.1391 155.4998      2138 112.01802 155.2477
## 30     10774 400.1391 155.4998      3497 138.99543 154.9399
## 31     10774 400.1391 155.4998      3892 143.03974 155.4589
## 32     10774 400.1391 155.4998      4026 143.99697 155.2523
## 33     10774 400.1391 155.4998      4858 158.96146 154.9756
## 34     10774 400.1391 155.4998      5113 161.96972 155.3870
## 35     10774 400.1391 155.4998      5759 177.00717 155.3517
## 36     10774 400.1391 155.4998      6877 208.03949 155.4236
## 37     10774 400.1391 155.4998     10066 358.16413 154.8946
## 38     10774 400.1391 155.4998     10795 401.14241 155.4359
## 39     10774 400.1391 155.4998     13081 785.28284 154.8956
## 40     10774 400.1391 155.4998     13106 799.27063 155.4821
## 41     10774 400.1391 155.4998     13108 800.27447 155.4060
## 42     10774 400.1391 155.4998     14986 190.92879 154.5651
## 43     10795 401.1424 155.4359      1627 109.02250 155.0114
## 44     10795 401.1424 155.4359      2138 112.01802 155.2477
## 45     10795 401.1424 155.4359      3497 138.99543 154.9399
## 46     10795 401.1424 155.4359      3892 143.03974 155.4589
## 47     10795 401.1424 155.4359      4026 143.99697 155.2523
## 48     10795 401.1424 155.4359      4725 156.99083 155.4705
## 49     10795 401.1424 155.4359      4858 158.96146 154.9756
## 50     10795 401.1424 155.4359      5113 161.96972 155.3870
## 51     10795 401.1424 155.4359      5759 177.00717 155.3517
## 52     10795 401.1424 155.4359      5853 180.98994 154.6402
## 53     10795 401.1424 155.4359      6877 208.03949 155.4236
## 54     10795 401.1424 155.4359     10036 356.18531 155.4821
## 55     10795 401.1424 155.4359     10066 358.16413 154.8946
## 56     10795 401.1424 155.4359     10774 400.13906 155.4998
## 57     10795 401.1424 155.4359     14986 190.92879 154.5651
## 58     13108 800.2745 155.4060      1627 109.02250 155.0114
## 59     13108 800.2745 155.4060      2138 112.01802 155.2477
## 60     13108 800.2745 155.4060      3497 138.99543 154.9399
## 61     13108 800.2745 155.4060      3892 143.03974 155.4589
## 62     13108 800.2745 155.4060      4026 143.99697 155.2523
## 63     13108 800.2745 155.4060      6877 208.03949 155.4236
## 64     13108 800.2745 155.4060     10036 356.18531 155.4821
## 65     13108 800.2745 155.4060     10774 400.13906 155.4998
## 66     13108 800.2745 155.4060     10795 401.14241 155.4359
## 67     13108 800.2745 155.4060     13081 785.28284 154.8956
## 68     13108 800.2745 155.4060     13106 799.27063 155.4821
## 69     14986 190.9288 154.5651      5113 161.96972 155.3870
## 70     14986 190.9288 154.5651     10036 356.18531 155.4821
## 71     14986 190.9288 154.5651     10066 358.16413 154.8946
```

Another type of dataframe can be generated for beta coefficients in model 2, we called this function `betaTable()`. This table collects beta coefficients of all correlated ions to every extracted extracted ions.

We can see the beta table from here. Since for illustration we used T at 155 second, below is the beta table for T==155.

``` r
beta_table<-betaTable(extract_ions,corr_ions)
beta_table
##          rtmed     mzmed 138.9954 156.9908 161.9697 180.9899 356.1853
## 637   155.2691  95.04397     0.00     0.00     0.00     0.00    -0.21
## 1306  155.4705 102.03395     0.00     0.00     0.00     0.00     0.00
## 1627  155.0114 109.02250     0.07     0.00     0.00     0.03     0.00
## 2138  155.2477 112.01802     0.00     0.00     0.00     0.00    -0.04
## 3497  154.9399 138.99543     0.99     0.00     0.00     0.08     0.13
## 3892  155.4589 143.03974     0.00     0.00     0.00     0.00     0.23
## 4026  155.2523 143.99697     0.00     0.00     0.00     0.00    -0.26
## 4725  155.4705 156.99083     0.00     0.99     0.00     0.00     0.02
## 4858  154.9756 158.96146     0.00     0.00     0.00     0.00     0.21
## 5113  155.3870 161.96972    -0.09     0.00     0.99     0.00     0.00
## 5759  155.3517 177.00717     0.00     0.00     0.00     0.00     0.28
## 5853  154.6402 180.98994     0.11     0.00     0.00     0.99     0.21
## 6877  155.4236 208.03949     0.00     0.00     0.00     0.00    -0.01
## 10036 155.4821 356.18531     0.05     0.00     0.00     0.00     0.99
## 10066 154.8946 358.16413     0.00     0.00     0.00     0.00     0.02
## 10774 155.4998 400.13906     0.00     0.00     0.00     0.00    -0.14
## 10795 155.4359 401.14241     0.00     0.00     0.00     0.00    -0.12
## 13081 154.8956 785.28284     0.00     0.00     0.00     0.00     0.07
## 13106 155.4821 799.27063     0.00     0.00     0.00     0.00     0.02
## 13108 155.4060 800.27447     0.00     0.00     0.00     0.00     0.10
## 14986 154.5651 190.92879     0.00     0.00     0.00     0.00     0.46
##       358.1641 400.1391 401.1424 800.2745 190.9288
## 637       0.00     0.00     0.00     0.00     0.00
## 1306      0.00     0.00     0.00     0.00     0.00
## 1627      0.00     0.16    -0.21     0.92     0.00
## 2138      0.00     0.02    -0.06    -0.28     0.00
## 3497      0.00     0.05    -0.07     0.22     0.00
## 3892      0.00    -0.14     0.20    -0.18     0.00
## 4026      0.00    -0.06     0.06    -0.14     0.00
## 4725      0.00     0.00    -0.01     0.00     0.00
## 4858      0.00    -0.11     0.13     0.00     0.00
## 5113      0.00    -0.04     0.03     0.00    -0.01
## 5759      0.00    -0.06     0.03     0.00     0.00
## 5853      0.00     0.00     0.02     0.00     0.00
## 6877      0.00     0.10    -0.14    -0.65     0.00
## 10036     0.00     0.00    -0.02     1.25     0.04
## 10066     0.99     0.00     0.00     0.00     0.00
## 10774     0.30     0.99     1.05     0.67     0.00
## 10795     0.17     0.89     0.99     0.19     0.00
## 13081     0.06     0.01     0.00     0.21     0.00
## 13106     0.07     0.01     0.00     0.15     0.00
## 13108     0.00     0.00     0.00     0.99     0.00
## 14986     0.00     0.01    -0.02     0.00     0.99
```

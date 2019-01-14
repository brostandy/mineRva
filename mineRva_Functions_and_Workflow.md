mineRva Functions and Workflow
================
Bety Rostandy
January 13, 2019

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
setwd("C:/Users/bty/OneDrive - UNCG/GitHub/mineRva")

#data<-read.csv("data/180517_Sample_pos_neg_sn3.csv",check.names="FALSE")
data<-read.csv("data/180517_Sample_pos_neg_sn10.csv",check.names="FALSE")

sample_group=sgroup(data,c("MeOH","03125","0625","125","250","500"))
dataset<-dataSet(data,sample_group,rt="rtmed",mz="mzmed")
x_method1<-xM1(dataset,sample_group)
msr<-msr(x_method1,sample_group)
outliers<-exclude5(x_method1,msr,T=155)
extract_ions<-ionx(xm=x_method1,sample_group,exclude=outliers,T=155,a=0.8) 
x_method2<-xM2(dataset,sample_group)
corr_ions<-ionr(xm=x_method2,extract_ions,sample_group,T=155,a=0.8)
```

After performing model 1 and model 2, an ion table of an *i*-th ion extracted and ion correlated to it can be produced with `ionTable()`.

We then are able to see model 1 and model 2 results in a dataframe where:

``` r
ion_table<-ionTable(extract_ions,corr_ions)
ion_table
##     inx_index    inx_mz   inx_rt inr_index    inr_mz   inr_rt
## 1         637  95.04397 155.2691      1306 102.03395 155.4705
## 2         637  95.04397 155.2691      4858 158.96146 154.9756
## 3         637  95.04397 155.2691      6877 208.03949 155.4236
## 4         637  95.04397 155.2691     13081 785.28284 154.8956
## 5        1306 102.03395 155.4705       637  95.04397 155.2691
## 6        1306 102.03395 155.4705       862  98.01370 155.1858
## 7        1306 102.03395 155.4705      1627 109.02250 155.0114
## 8        1306 102.03395 155.4705      2138 112.01802 155.2477
## 9        1306 102.03395 155.4705      6877 208.03949 155.4236
## 10       1306 102.03395 155.4705     10795 401.14241 155.4359
## 11       1306 102.03395 155.4705     13081 785.28284 154.8956
## 12       1306 102.03395 155.4705     13108 800.27447 155.4060
## 13       1627 109.02250 155.0114       637  95.04397 155.2691
## 14       1627 109.02250 155.0114      1306 102.03395 155.4705
## 15       1627 109.02250 155.0114      2138 112.01802 155.2477
## 16       1627 109.02250 155.0114      3497 138.99543 154.9399
## 17       1627 109.02250 155.0114      3892 143.03974 155.4589
## 18       1627 109.02250 155.0114      4026 143.99697 155.2523
## 19       1627 109.02250 155.0114      4725 156.99083 155.4705
## 20       1627 109.02250 155.0114      4858 158.96146 154.9756
## 21       1627 109.02250 155.0114      5113 161.96972 155.3870
## 22       1627 109.02250 155.0114      5759 177.00717 155.3517
## 23       1627 109.02250 155.0114      5853 180.98994 154.6402
## 24       1627 109.02250 155.0114      6877 208.03949 155.4236
## 25       1627 109.02250 155.0114     10036 356.18531 155.4821
## 26       1627 109.02250 155.0114     10066 358.16413 154.8946
## 27       1627 109.02250 155.0114     10774 400.13906 155.4998
## 28       1627 109.02250 155.0114     10795 401.14241 155.4359
## 29       1627 109.02250 155.0114     13081 785.28284 154.8956
## 30       1627 109.02250 155.0114     13106 799.27063 155.4821
## 31       1627 109.02250 155.0114     13108 800.27447 155.4060
## 32       1627 109.02250 155.0114     14986 190.92879 154.5651
## 33       3497 138.99543 154.9399      1627 109.02250 155.0114
## 34       3497 138.99543 154.9399      5113 161.96972 155.3870
## 35       3497 138.99543 154.9399      5853 180.98994 154.6402
## 36       3497 138.99543 154.9399     10036 356.18531 155.4821
## 37       4725 156.99083 155.4705      1306 102.03395 155.4705
## 38       4725 156.99083 155.4705      2138 112.01802 155.2477
## 39       4725 156.99083 155.4705      3892 143.03974 155.4589
## 40       4725 156.99083 155.4705      5759 177.00717 155.3517
## 41       4725 156.99083 155.4705     10066 358.16413 154.8946
## 42       4725 156.99083 155.4705     10774 400.13906 155.4998
## 43       4725 156.99083 155.4705     10795 401.14241 155.4359
## 44       4858 158.96146 154.9756       637  95.04397 155.2691
## 45       4858 158.96146 154.9756      6877 208.03949 155.4236
## 46       5113 161.96972 155.3870      3497 138.99543 154.9399
## 47       5113 161.96972 155.3870      3892 143.03974 155.4589
## 48       5759 177.00717 155.3517       862  98.01370 155.1858
## 49       5759 177.00717 155.3517      1627 109.02250 155.0114
## 50       5759 177.00717 155.3517      3497 138.99543 154.9399
## 51       5759 177.00717 155.3517      4026 143.99697 155.2523
## 52       5759 177.00717 155.3517      4725 156.99083 155.4705
## 53       5759 177.00717 155.3517      4858 158.96146 154.9756
## 54       5759 177.00717 155.3517      5853 180.98994 154.6402
## 55       5759 177.00717 155.3517      6877 208.03949 155.4236
## 56       5759 177.00717 155.3517     10036 356.18531 155.4821
## 57       5759 177.00717 155.3517     10066 358.16413 154.8946
## 58       5759 177.00717 155.3517     10774 400.13906 155.4998
## 59       5759 177.00717 155.3517     13081 785.28284 154.8956
## 60       5759 177.00717 155.3517     14986 190.92879 154.5651
## 61       5853 180.98994 154.6402      1627 109.02250 155.0114
## 62       5853 180.98994 154.6402      3497 138.99543 154.9399
## 63       6877 208.03949 155.4236       637  95.04397 155.2691
## 64       6877 208.03949 155.4236      1306 102.03395 155.4705
## 65       6877 208.03949 155.4236      2138 112.01802 155.2477
## 66       6877 208.03949 155.4236      4858 158.96146 154.9756
## 67       6877 208.03949 155.4236      5113 161.96972 155.3870
## 68       6877 208.03949 155.4236      5759 177.00717 155.3517
## 69       6877 208.03949 155.4236     13108 800.27447 155.4060
## 70      10066 358.16413 154.8946     10774 400.13906 155.4998
## 71      10066 358.16413 154.8946     10795 401.14241 155.4359
## 72      10066 358.16413 154.8946     13081 785.28284 154.8956
## 73      10066 358.16413 154.8946     13106 799.27063 155.4821
## 74      10066 358.16413 154.8946     14986 190.92879 154.5651
## 75      10774 400.13906 155.4998      1627 109.02250 155.0114
## 76      10774 400.13906 155.4998      3497 138.99543 154.9399
## 77      10774 400.13906 155.4998      3892 143.03974 155.4589
## 78      10774 400.13906 155.4998      4026 143.99697 155.2523
## 79      10774 400.13906 155.4998      4858 158.96146 154.9756
## 80      10774 400.13906 155.4998      5113 161.96972 155.3870
## 81      10774 400.13906 155.4998      5759 177.00717 155.3517
## 82      10774 400.13906 155.4998      6877 208.03949 155.4236
## 83      10774 400.13906 155.4998     10066 358.16413 154.8946
## 84      10774 400.13906 155.4998     10795 401.14241 155.4359
## 85      10774 400.13906 155.4998     13081 785.28284 154.8956
## 86      10774 400.13906 155.4998     13106 799.27063 155.4821
## 87      10774 400.13906 155.4998     13108 800.27447 155.4060
## 88      10795 401.14241 155.4359      1306 102.03395 155.4705
## 89      10795 401.14241 155.4359      1627 109.02250 155.0114
## 90      10795 401.14241 155.4359      2138 112.01802 155.2477
## 91      10795 401.14241 155.4359      3497 138.99543 154.9399
## 92      10795 401.14241 155.4359      3892 143.03974 155.4589
## 93      10795 401.14241 155.4359      4026 143.99697 155.2523
## 94      10795 401.14241 155.4359      4858 158.96146 154.9756
## 95      10795 401.14241 155.4359      5113 161.96972 155.3870
## 96      10795 401.14241 155.4359      6877 208.03949 155.4236
## 97      10795 401.14241 155.4359     10774 400.13906 155.4998
## 98      10795 401.14241 155.4359     14986 190.92879 154.5651
## 99      13081 785.28284 154.8956      5759 177.00717 155.3517
## 100     13081 785.28284 154.8956     10036 356.18531 155.4821
## 101     13081 785.28284 154.8956     10066 358.16413 154.8946
## 102     13081 785.28284 154.8956     10774 400.13906 155.4998
## 103     13081 785.28284 154.8956     10795 401.14241 155.4359
## 104     13081 785.28284 154.8956     13106 799.27063 155.4821
## 105     13081 785.28284 154.8956     13108 800.27447 155.4060
## 106     13108 800.27447 155.4060      1627 109.02250 155.0114
## 107     13108 800.27447 155.4060      2138 112.01802 155.2477
## 108     13108 800.27447 155.4060      3497 138.99543 154.9399
## 109     13108 800.27447 155.4060      3892 143.03974 155.4589
## 110     13108 800.27447 155.4060      4026 143.99697 155.2523
## 111     13108 800.27447 155.4060      6877 208.03949 155.4236
## 112     13108 800.27447 155.4060     10036 356.18531 155.4821
## 113     13108 800.27447 155.4060     10774 400.13906 155.4998
## 114     13108 800.27447 155.4060     10795 401.14241 155.4359
## 115     13108 800.27447 155.4060     13081 785.28284 154.8956
## 116     13108 800.27447 155.4060     13106 799.27063 155.4821
```

Another type of dataframe can be generated for beta coefficients in model 2, we called this function `betaTable()`. This table collects beta coefficients of all correlated ions to every extracted extracted ions.

We can see the beta table from here. Since for illustration we used T at 155 second, below is the beta table for T==155.

``` r
beta_table<-betaTable(extract_ions,corr_ions)
beta_table
##          rtmed     mzmed 95.044 102.034 109.0225 138.9954 143.0397 143.997
## 637   155.2691  95.04397   0.99    0.05     0.01     0.00     0.00    0.00
## 1306  155.4705 102.03395   0.01    0.99     0.34     0.00     0.00    0.00
## 1627  155.0114 109.02250   0.00    0.08     0.99     0.08     0.00    0.00
## 2138  155.2477 112.01802   0.00    0.15     0.09     0.00     0.00    0.00
## 3497  154.9399 138.99543   0.00    0.00     0.04     0.99     0.00    0.00
## 3892  155.4589 143.03974   0.00    0.00     0.58     0.00     0.99    0.00
## 4026  155.2523 143.99697   0.00    0.00    -0.08     0.00     0.00    0.99
## 4725  155.4705 156.99083   0.00    0.00    -0.08     0.00     0.00    0.00
## 4858  154.9756 158.96146  -0.06    0.00     0.16     0.00     0.00    0.00
## 5113  155.3870 161.96972   0.00    0.00    -0.15    -0.12     0.00    0.00
## 5759  155.3517 177.00717   0.00    0.00     0.20     0.00     0.00    0.00
## 5853  154.6402 180.98994   0.00    0.00     0.24     0.13     0.00    0.00
## 6877  155.4236 208.03949   0.22    0.11    -0.07     0.00     0.00    0.00
## 10036 155.4821 356.18531   0.00    0.00    -0.10     0.06     0.00    0.00
## 10066 154.8946 358.16413   0.00    0.00     0.01     0.00     0.00    0.00
## 10774 155.4998 400.13906   0.00    0.00     0.60     0.00     0.00    0.00
## 10795 155.4359 401.14241   0.00   -0.04    -0.64     0.00     0.00    0.00
## 13081 154.8956 785.28284   0.00   -0.01    -0.01     0.00     0.00    0.00
## 13106 155.4821 799.27063   0.00    0.00    -0.01     0.00     0.00    0.00
## 13108 155.4060 800.27447   0.00    0.00     0.06     0.00     0.00    0.00
## 14986 154.5651 190.92879   0.00    0.00     0.06     0.00     0.00    0.00
##       156.9908 158.9615 161.9697 177.0072 180.9899 208.0395 358.1641
## 637       0.00    -0.15     0.00     0.00     0.00     0.68     0.00
## 1306      0.01     0.00     0.00     0.00     0.00     0.11     0.00
## 1627      0.00     0.00     0.00     0.23     0.00     0.00     0.00
## 2138      0.16     0.00     0.00     0.00     0.00     0.01     0.00
## 3497      0.00     0.00    -0.04     0.03     0.05     0.00     0.00
## 3892      0.03     0.00     0.04     0.00     0.00     0.00     0.00
## 4026      0.00     0.00     0.00     0.31     0.00     0.00     0.00
## 4725      0.99     0.00     0.00     0.20     0.00     0.00     0.00
## 4858      0.00     0.99     0.00     0.03     0.00    -0.11     0.00
## 5113      0.00     0.00     0.99     0.00     0.00     0.00     0.00
## 5759      0.07     0.00     0.00     0.99     0.00    -0.05     0.00
## 5853      0.00     0.00     0.00     0.06     0.99     0.00     0.00
## 6877      0.00    -0.08     0.00    -0.16     0.00     0.99     0.00
## 10036     0.00     0.00     0.00     0.14     0.00     0.00     0.00
## 10066     0.00     0.00     0.00    -0.02     0.00     0.00     0.99
## 10774    -0.02     0.00     0.00    -0.03     0.00     0.00     0.31
## 10795     0.00     0.00     0.00     0.00     0.00     0.00     0.17
## 13081     0.00     0.00     0.00     0.01     0.00     0.00     0.07
## 13106     0.00     0.00     0.00     0.00     0.00     0.00     0.07
## 13108     0.00     0.00     0.00     0.00     0.00    -0.01     0.00
## 14986     0.00     0.00     0.00    -0.11     0.00     0.00    -0.04
##       400.1391 401.1424 785.2828 800.2745 190.9288
## 637       0.00     0.00     0.00     0.00     0.00
## 1306      0.00    -0.01     0.00     0.00     0.00
## 1627      0.13    -0.18     0.00     0.86     0.00
## 2138      0.00    -0.04     0.00    -0.21     0.00
## 3497      0.03    -0.05     0.00     0.19     0.00
## 3892     -0.11     0.14     0.00    -0.06     0.00
## 4026     -0.07     0.06     0.00    -0.08     0.00
## 4725      0.00     0.00     0.00     0.00     0.00
## 4858     -0.11     0.11     0.00     0.00     0.00
## 5113     -0.04     0.01     0.00     0.00     0.00
## 5759     -0.04     0.00     0.26     0.00     0.00
## 5853      0.00     0.00     0.00     0.00     0.00
## 6877      0.08    -0.13     0.00    -0.61     0.00
## 10036     0.00     0.00     0.69     1.25     0.00
## 10066     0.00     0.00     0.19     0.00     0.00
## 10774     0.99     1.04     0.56     0.66     0.00
## 10795     0.88     0.99     0.09     0.18     0.00
## 13081     0.01     0.00     0.99     0.22     0.00
## 13106     0.01     0.00     0.24     0.15     0.00
## 13108     0.00     0.00     0.24     0.99     0.00
## 14986     0.00    -0.02     0.00     0.00     0.99
```

The `betaTable()` result can also be viewed with `ionHeatmap()`. The extracted ions from `ionx()` can be viewed with `ionxPlot()` and the correlated ions from `ionr()` can be viewed with `ionrPlot()`.

Therefore,

``` r
ion_heatmap<-ionHeatmap(extract_ions,corr_ions)
```

![](mineRva_Functions_and_Workflow_files/figure-markdown_github/unnamed-chunk-17-1.png)

``` r
extracted_ion_plots<-ionxPlot(xm1=x_method1,inx=extract_ions,sgrp=sample_group)
```

![](mineRva_Functions_and_Workflow_files/figure-markdown_github/unnamed-chunk-17-2.png)

    ## Press [enter] to continue

![](mineRva_Functions_and_Workflow_files/figure-markdown_github/unnamed-chunk-17-3.png)

    ## Press [enter] to continue

![](mineRva_Functions_and_Workflow_files/figure-markdown_github/unnamed-chunk-17-4.png)

    ## Press [enter] to continue

![](mineRva_Functions_and_Workflow_files/figure-markdown_github/unnamed-chunk-17-5.png)

    ## Press [enter] to continue

![](mineRva_Functions_and_Workflow_files/figure-markdown_github/unnamed-chunk-17-6.png)

    ## Press [enter] to continue

![](mineRva_Functions_and_Workflow_files/figure-markdown_github/unnamed-chunk-17-7.png)

    ## Press [enter] to continue

![](mineRva_Functions_and_Workflow_files/figure-markdown_github/unnamed-chunk-17-8.png)

    ## Press [enter] to continue

![](mineRva_Functions_and_Workflow_files/figure-markdown_github/unnamed-chunk-17-9.png)

    ## Press [enter] to continue

![](mineRva_Functions_and_Workflow_files/figure-markdown_github/unnamed-chunk-17-10.png)

    ## Press [enter] to continue

![](mineRva_Functions_and_Workflow_files/figure-markdown_github/unnamed-chunk-17-11.png)

    ## Press [enter] to continue

![](mineRva_Functions_and_Workflow_files/figure-markdown_github/unnamed-chunk-17-12.png)

    ## Press [enter] to continue

![](mineRva_Functions_and_Workflow_files/figure-markdown_github/unnamed-chunk-17-13.png)

    ## Press [enter] to continue

![](mineRva_Functions_and_Workflow_files/figure-markdown_github/unnamed-chunk-17-14.png)

    ## Press [enter] to continue

![](mineRva_Functions_and_Workflow_files/figure-markdown_github/unnamed-chunk-17-15.png)

    ## Press [enter] to continue

![](mineRva_Functions_and_Workflow_files/figure-markdown_github/unnamed-chunk-17-16.png)

    ## Press [enter] to continue

![](mineRva_Functions_and_Workflow_files/figure-markdown_github/unnamed-chunk-17-17.png)

    ## Press [enter] to continue

![](mineRva_Functions_and_Workflow_files/figure-markdown_github/unnamed-chunk-17-18.png)

    ## Press [enter] to continue

![](mineRva_Functions_and_Workflow_files/figure-markdown_github/unnamed-chunk-17-19.png)

    ## Press [enter] to continue

``` r
correlated_ion_plots<-ionrPlot(xm2=x_method2,inr=corr_ions,sgrp=sample_group,frac=0.2) 
```

![](mineRva_Functions_and_Workflow_files/figure-markdown_github/unnamed-chunk-17-20.png)

    ## Press [enter] to continue

![](mineRva_Functions_and_Workflow_files/figure-markdown_github/unnamed-chunk-17-21.png)

    ## Press [enter] to continue

![](mineRva_Functions_and_Workflow_files/figure-markdown_github/unnamed-chunk-17-22.png)

    ## Press [enter] to continue

![](mineRva_Functions_and_Workflow_files/figure-markdown_github/unnamed-chunk-17-23.png)

    ## Press [enter] to continue

![](mineRva_Functions_and_Workflow_files/figure-markdown_github/unnamed-chunk-17-24.png)

    ## Press [enter] to continue

![](mineRva_Functions_and_Workflow_files/figure-markdown_github/unnamed-chunk-17-25.png)

    ## Press [enter] to continue

![](mineRva_Functions_and_Workflow_files/figure-markdown_github/unnamed-chunk-17-26.png)

    ## Press [enter] to continue

![](mineRva_Functions_and_Workflow_files/figure-markdown_github/unnamed-chunk-17-27.png)

    ## Press [enter] to continue
    ## [1] "No plot for m/z= 143.0397"
    ## [1] "No plot for m/z= 143.997"

![](mineRva_Functions_and_Workflow_files/figure-markdown_github/unnamed-chunk-17-28.png)

    ## Press [enter] to continue
    ## [1] "No plot for m/z= 158.9615"
    ## [1] "No plot for m/z= 161.9697"

![](mineRva_Functions_and_Workflow_files/figure-markdown_github/unnamed-chunk-17-29.png)

    ## Press [enter] to continue

![](mineRva_Functions_and_Workflow_files/figure-markdown_github/unnamed-chunk-17-30.png)

    ## Press [enter] to continue

![](mineRva_Functions_and_Workflow_files/figure-markdown_github/unnamed-chunk-17-31.png)

    ## Press [enter] to continue
    ## [1] "No plot for m/z= 180.9899"

![](mineRva_Functions_and_Workflow_files/figure-markdown_github/unnamed-chunk-17-32.png)

    ## Press [enter] to continue

![](mineRva_Functions_and_Workflow_files/figure-markdown_github/unnamed-chunk-17-33.png)

    ## Press [enter] to continue

![](mineRva_Functions_and_Workflow_files/figure-markdown_github/unnamed-chunk-17-34.png)

    ## Press [enter] to continue

![](mineRva_Functions_and_Workflow_files/figure-markdown_github/unnamed-chunk-17-35.png)

    ## Press [enter] to continue

![](mineRva_Functions_and_Workflow_files/figure-markdown_github/unnamed-chunk-17-36.png)

    ## Press [enter] to continue

![](mineRva_Functions_and_Workflow_files/figure-markdown_github/unnamed-chunk-17-37.png)

    ## Press [enter] to continue

![](mineRva_Functions_and_Workflow_files/figure-markdown_github/unnamed-chunk-17-38.png)

    ## Press [enter] to continue

![](mineRva_Functions_and_Workflow_files/figure-markdown_github/unnamed-chunk-17-39.png)

    ## Press [enter] to continue

![](mineRva_Functions_and_Workflow_files/figure-markdown_github/unnamed-chunk-17-40.png)

    ## Press [enter] to continue

![](mineRva_Functions_and_Workflow_files/figure-markdown_github/unnamed-chunk-17-41.png)

    ## Press [enter] to continue
    ## [1] "No plot for m/z= 190.9288"

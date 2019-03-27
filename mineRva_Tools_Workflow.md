mineRva Tools\_Workflow
================
Bety Rostandy
March 26, 2019

### What is mineRva?

mineRva is a tool built in R environment to extract metabolite ions for a complex mixture of LC-MS analysis after preprocessing steps (i.e peaks filling, peaks allignment, retention time correction, and so on). The analytical statistics model is built for LC-MS complex (i.e. biological) mixture dataset that must contain retention time, mass-to-charge ratio, and measured ions' intensity from samples. This tool is built to analyze a specific design of experiment in measuring biomolecules via LC-MS system.

### Metabolite Ions Extraction and Corelation Functions

Users feed in preprocessed LC-MS data, where ions' peak had been deconvoluted to a table. Then, to enable extraction of dataset from data, we are going to define sample groups with `sgroup()`. Users give groupnames based on the column names that define the sample groups of preprocessed data. The first group is always the negative control or blanks.

After that, the dataset is extracted and zero-entries are treated with `dataSet()`. This dataset is extracted and its zero-entries are treated the same for all model, however, each model was rescaled and transformed accordingly.

Therefore to rescale and transform the dataset for model 1, we use `xM1()`.

Then, users can choose to perform RSD cutoff to determine outlier (only for Model 1, as Model 2 uses Model 1 result). If users choose to perform RSD cutoff an exclude list will be returned. This exclude list follows certain conditions that is set for five sample groups experiment (i.e. six *groupnames* with blank).The outlier for RSD cutoff is determined at such: *X* &gt; *Q*3 + 1.5 \* *I**Q**R* or *X* &lt; *Q*1 − 1.5 \* *I**Q**R*. In order to get RSD, mean and standard deviation calculations must be performed. Therefore, `msr()` performs calculations for mean(*m*), standard deviation (*s*) and relative standard deviation (*r*).

The exclude list is returned after running the function `exclude5()`. (Note: as of current development, ions exclusion can only be done if users have five sample groups and one negative control.)

Finally, Model 1 can be performed (with or without outlier removal) with `ionx()`. Model 1 has the objective of extracting "important" ions from LC-MS data.

After fitting model 1, where metabolite ions are extracted at a retention time group (*rtgrp*) given by the user, these extracted ions of model 1 then fitted through model 2 to find its related ions. Ion intensity is normalized for model 2 beforehand with xM2()\`.

Model 2 function, `ionr()`, is to perform ion-ion relationship detection that were extracted from Model 1. For each extracted ion from Model 1, we against it with all other ions at the given retention time group.

User can call model 1 and model 2 functions consecutively as below:

After performing model 1 and model 2, an ion table of an *inx*-th of ion extracted and ion related to it can be produced with `ionTable()`.

We then are able to see model 1 and model 2 results in a dataframe where:

    ##    inx_index   inx_mz   inx_rt inr_index   inr_mz   inr_rt
    ## 1      16620 297.0774 423.1332      7907 259.2827 422.6997
    ## 2      16620 297.0774 423.1332      9630 338.1935 422.7710
    ## 3      16620 297.0774 423.1332     11372 445.2171 422.5748
    ## 4      16620 297.0774 423.1332     11387 446.2211 422.8140
    ## 5       8802 299.0916 423.2699      7907 259.2827 422.6997
    ## 6       8802 299.0916 423.2699      9630 338.1935 422.7710
    ## 7       8802 299.0916 423.2699     11372 445.2171 422.5748
    ## 8       8802 299.0916 423.2699     11387 446.2211 422.8140
    ## 9       9610 337.1907 422.5514      7907 259.2827 422.6997
    ## 10      9610 337.1907 422.5514      9630 338.1935 422.7710
    ## 11      9610 337.1907 422.5514     11372 445.2171 422.5748
    ## 12      9610 337.1907 422.5514     11387 446.2211 422.8140
    ## 13      8826 300.0949 423.3019      7907 259.2827 422.6997
    ## 14      8826 300.0949 423.3019      9630 338.1935 422.7710
    ## 15      8826 300.0949 423.3019     11372 445.2171 422.5748
    ## 16      8826 300.0949 423.3019     11387 446.2211 422.8140
    ## 17     16649 298.0808 423.1897      7907 259.2827 422.6997
    ## 18     16649 298.0808 423.1897      9630 338.1935 422.7710
    ## 19     16649 298.0808 423.1897     11372 445.2171 422.5748
    ## 20     16649 298.0808 423.1897     11387 446.2211 422.8140
    ## 21      9232 318.3005 423.4211      7907 259.2827 422.6997
    ## 22      9232 318.3005 423.4211      9630 338.1935 422.7710
    ## 23      9232 318.3005 423.4211     11372 445.2171 422.5748
    ## 24      9232 318.3005 423.4211     11387 446.2211 422.8140
    ## 25     13013 754.2136 423.1698      7907 259.2827 422.6997
    ## 26     13013 754.2136 423.1698      9630 338.1935 422.7710
    ## 27     13013 754.2136 423.1698     11372 445.2171 422.5748
    ## 28     13013 754.2136 423.1698     11387 446.2211 422.8140
    ## 29     11372 445.2171 422.5748      7907 259.2827 422.6997
    ## 30     11372 445.2171 422.5748      9630 338.1935 422.7710
    ## 31     11372 445.2171 422.5748     11372 445.2171 422.5748
    ## 32     11372 445.2171 422.5748     11387 446.2211 422.8140
    ## 33      5770 177.0885 422.5229      7907 259.2827 422.6997
    ## 34      5770 177.0885 422.5229      9630 338.1935 422.7710
    ## 35      5770 177.0885 422.5229     11372 445.2171 422.5748
    ## 36      5770 177.0885 422.5229     11387 446.2211 422.8140

Another type of dataframe can be generated for beta coefficients in model 2, we called this function `betaTable()`. This table collects beta coefficients of all related ions to every extracted ions.

We can see the beta table from here. Since for illustration we used retention time group of 423 second, below is the beta table for RT = 422.5 - 423.5.

    ##       ori_index    rtmed     mzmed 177.0885 299.0916 300.0949 318.3005
    ## 2183       2183 423.4343 112.51974      0.0      0.0      0.0      0.0
    ## 2226       2226 422.6815 112.96691      0.0      0.0      0.0      0.0
    ## 2621       2621 422.9557 117.95975      0.0      0.0      0.0      0.0
    ## 3657       3657 423.3983 140.96185      0.0      0.0      0.0      0.0
    ## 3943       3943 423.0972 143.95911      0.0      0.0      0.0      0.0
    ## 5369       5369 422.7492 168.01635      0.0      0.0      0.0      0.0
    ## 5770       5770 422.5229 177.08849      0.9      0.0      0.0      0.0
    ## 7434       7434 422.9004 230.12009      0.0      0.0      0.0      0.0
    ## 7662       7662 422.8300 245.12535      0.0      0.0      0.0      0.0
    ## 7672       7672 422.7905 245.62740      0.0      0.0      0.0      0.0
    ## 7888       7888 422.6081 258.27934      0.0      0.0      0.0      0.0
    ## 7907       7907 422.6997 259.28271      0.0      0.0      0.0      0.1
    ## 8050       8050 423.2421 269.16500     -0.2      0.0      0.0      0.0
    ## 8304       8304 422.8505 276.14375      0.0      0.0      0.0      0.0
    ## 8445       8445 423.1034 283.15151      0.0      0.0      0.0      0.0
    ## 8461       8461 422.8483 283.65338      0.0      0.0      0.0      0.0
    ## 8472       8472 422.8505 284.15474      0.0      0.0      0.0      0.0
    ## 8625       8625 422.6625 292.15672      0.0      0.0      0.0      0.0
    ## 8802       8802 423.2699 299.09157      0.0      0.9      0.6      0.0
    ## 8808       8808 422.5365 299.16468      0.0      0.0      0.0      0.0
    ## 8826       8826 423.3019 300.09494      0.0      0.1      0.9      0.0
    ## 8976       8976 423.2891 308.18544      0.0      0.0      0.0      0.0
    ## 9232       9232 423.4211 318.30048      0.0      0.0      0.0      0.9
    ## 9259       9259 423.4343 319.30392      0.0      0.0      0.0      0.5
    ## 9278       9278 422.8300 321.17771      0.0      0.0      0.0      0.0
    ## 9281       9281 422.8140 321.67944      0.0      0.0      0.0      0.0
    ## 9289       9289 422.9738 322.18015      0.0      0.0      0.0      0.0
    ## 9462       9462 422.8300 330.18298      0.0      0.0      0.0      0.0
    ## 9469       9469 422.8483 330.68472      0.0      0.0      0.0      0.0
    ## 9600       9600 422.7244 337.10720      0.0      0.0      0.0      0.0
    ## 9610       9610 422.5514 337.19068      0.0      0.0      0.0      0.0
    ## 9618       9618 422.6290 337.69242      0.0      0.0      0.0      0.0
    ## 9630       9630 422.7710 338.19355      0.0      0.0      0.0      0.0
    ## 9800       9800 423.0138 345.26369      0.0      0.0      0.0      0.1
    ## 9815       9815 422.7905 346.19612      0.0      0.0      0.0      0.0
    ## 9831       9831 422.5514 346.69783      0.0      0.0      0.0      0.0
    ## 9950       9950 422.5982 353.70557      0.0      0.0      0.0      0.0
    ## 10149     10149 423.3510 362.29020      0.0      0.0      0.0      0.0
    ## 10238     10238 423.3771 367.24576      0.0      0.0      0.0      0.0
    ## 10277     10277 423.0501 369.16459      0.0      0.0      0.0      0.0
    ## 10429     10429 423.0988 376.30593      0.0      0.0      0.0      0.0
    ## 10573     10573 422.6085 387.19313      0.0      0.0      0.0      0.0
    ## 11372     11372 422.5748 445.21707      0.0      0.0      0.0      0.0
    ## 11387     11387 422.8140 446.22105      0.1      0.0      0.0      0.0
    ## 11745     11745 422.6826 487.26833      0.0      0.0      0.0      0.0
    ## 12149     12149 422.5875 542.39047      0.0      0.0      0.0      0.0
    ## 12795     12795 422.6085 673.37429      0.0      0.0      0.0      0.0
    ## 12800     12800 422.5371 674.37728      0.0      0.0      0.0      0.0
    ## 12808     12808 422.5982 675.37979      0.0      0.0      0.0      0.0
    ## 12909     12909 422.5514 705.40064      0.0      0.0      0.0      0.0
    ## 12912     12912 422.5229 706.40372      0.0      0.0      0.0      0.0
    ## 12915     12915 422.5982 707.40628      0.0      0.0      0.0      0.0
    ## 12917     12917 422.7671 708.40903      0.0      0.0      0.0      0.0
    ## 12963     12963 422.5217 727.38257      0.0      0.0      0.0      0.0
    ## 12966     12966 422.5748 728.38575      0.0      0.0      0.0      0.0
    ## 13013     13013 423.1698 754.21363      0.0      0.0      0.0      0.0
    ## 13016     13016 423.2487 755.21696      0.0      0.1      0.0      0.0
    ## 13191     13191 423.4472 930.39890      0.0      0.1      0.1      0.0
    ## 13407     13407 423.1755  98.93943      0.0      0.0      0.0      0.0
    ## 13437     13437 422.7354  99.92609      0.0      0.0      0.0     -0.1
    ## 13861     13861 422.7129 135.89510      0.0      0.0      0.0      0.0
    ## 13944     13944 422.6523 142.92929      0.0      0.0      0.0      0.0
    ## 15099     15099 423.4712 194.92461      0.0      0.0      0.0      0.0
    ## 16620     16620 423.1332 297.07739      0.0      0.1      0.0      0.0
    ## 16649     16649 423.1897 298.08079      0.0      0.4      0.0      0.0
    ## 16658     16658 422.7984 298.15696      0.0      0.0      0.0      0.0
    ## 18022     18022 423.3396 449.00633      0.0      0.0      0.0      0.0
    ##       337.1907 445.2171 754.2136 297.0774 298.0808
    ## 2183       0.0      0.0      0.0     -0.2      0.1
    ## 2226       0.0      0.0      0.0      0.0      0.0
    ## 2621       0.0      0.1      0.0      0.0      0.0
    ## 3657       0.0      0.0      0.0      0.0      0.0
    ## 3943       0.0      0.0      0.0      0.0      0.0
    ## 5369       0.0      0.0      0.0      0.0      0.0
    ## 5770       0.0      0.4      0.0      0.0      0.0
    ## 7434       0.0      0.0      0.0     -0.4      0.0
    ## 7662       0.0      0.0      0.0      0.0      0.0
    ## 7672       0.0      0.0      0.0      0.0      0.0
    ## 7888       0.0      0.0      0.0      0.0      0.0
    ## 7907      -0.1      0.1      0.0      0.0      0.0
    ## 8050       0.0      0.0      0.0      0.0      0.0
    ## 8304       0.0      0.0      0.0      0.0      0.0
    ## 8445       0.0      0.0      0.0      0.0      0.0
    ## 8461       0.0      0.0      0.0      0.0      0.0
    ## 8472       0.0      0.1      0.0      0.0      0.0
    ## 8625       0.0      0.0      0.0      0.0      0.0
    ## 8802       0.0      0.0      0.0      0.1      0.3
    ## 8808       0.0      0.0      0.0      0.0      0.0
    ## 8826       0.0      0.0      0.0      0.0      0.0
    ## 8976       0.0      0.0      0.0      0.0      0.0
    ## 9232       0.0      0.1      0.0      0.0      0.0
    ## 9259       0.0      0.0      0.0      0.0      0.0
    ## 9278       0.0      0.0      0.0      0.0      0.0
    ## 9281       0.0      0.0      0.0      0.0      0.0
    ## 9289       0.0      0.0      0.0      0.0      0.0
    ## 9462       0.0      0.0      0.0     -0.1      0.0
    ## 9469       0.0      0.0     -0.1      0.0      0.0
    ## 9600       0.0      0.0      0.0      0.0      0.0
    ## 9610       0.9      0.7      0.0      0.0      0.0
    ## 9618       0.0      0.0      0.0      0.0      0.0
    ## 9630       0.4      0.3      0.0      0.3      0.1
    ## 9800       0.0      0.0      0.0      0.0      0.0
    ## 9815       0.0      0.1     -1.7      0.0      0.0
    ## 9831       0.0      0.0      0.0      0.0      0.0
    ## 9950       0.0      0.0     -0.5      0.0      0.0
    ## 10149      0.0      0.0      0.0      0.0      0.0
    ## 10238      0.0      0.0      0.0      0.1      0.0
    ## 10277      0.0      0.0     -0.2      0.0      0.0
    ## 10429      0.0      0.0      0.0      0.0      0.0
    ## 10573      0.0      0.0     -0.5      0.0      0.0
    ## 11372      0.1      0.9      0.0      0.0      0.0
    ## 11387      0.2      0.2      0.0      0.0      0.0
    ## 11745      0.0      0.0      0.0      0.0      0.0
    ## 12149      0.0      0.0      0.0      0.0      0.0
    ## 12795      0.0      0.1      0.0      0.0      0.0
    ## 12800      0.0      0.0      0.0      0.0      0.0
    ## 12808      0.0      0.0      0.0      0.0      0.0
    ## 12909      0.0      0.0      0.0      0.0      0.0
    ## 12912      0.0      0.0      0.0      0.0      0.0
    ## 12915      0.0      0.0      0.0      0.0      0.0
    ## 12917      0.0      0.0      0.0      0.0      0.0
    ## 12963      0.0      0.0      0.0      0.0      0.0
    ## 12966      0.0      0.0      0.0      0.0      0.0
    ## 13013      0.0      0.0      0.9      0.0      0.0
    ## 13016      0.0      0.0      0.1      0.0      0.0
    ## 13191      0.0      0.0      0.0      0.0      0.0
    ## 13407      0.0      0.0      0.0      0.0      0.0
    ## 13437      0.0      0.3      0.0      0.0      0.1
    ## 13861      0.0      0.0      0.0      0.0      0.0
    ## 13944      0.0      0.0      0.0     -0.1      0.2
    ## 15099      0.0      0.0      0.0      0.0      0.0
    ## 16620      0.0      0.0      0.0      0.9      0.7
    ## 16649      0.0      0.0      0.7      0.6      0.9
    ## 16658      0.0      0.0      0.0      0.0      0.0
    ## 18022      0.0      0.0      0.1      0.0      0.0

The `betaTable()` result can also be viewed with `ionHeatmap()`. Ions that are related to the extracted ions are simply denoted with: *1* if their relationship is directly proportional, *-1* if their relationship is inversely proportional and *0* if there is no relationship. The extracted ions from `ionx()` can be viewed with `ionxPlot()` and the ion-ion relationship from `ionr()` can be viewed with `ionrPlot()`.

Therefore,

``` r
ion_heatmap<-ionHeatmap(extract_ions,relate_ions)
```

![](mineRva_Tools_Workflow/unnamed-chunk-3-1.png)

``` r
extracted_ion_plots<-ionxPlot(xm1=x_method1,inx=extract_ions,sgrp=sample_group)
```

![](mineRva_Tools_Workflow/unnamed-chunk-3-2.png)

    ## Press [enter] to continue

![](mineRva_Tools_Workflow/unnamed-chunk-3-3.png)

    ## Press [enter] to continue

![](mineRva_Tools_Workflow/unnamed-chunk-3-4.png)

    ## Press [enter] to continue

![](mineRva_Tools_Workflow/unnamed-chunk-3-5.png)

    ## Press [enter] to continue

![](mineRva_Tools_Workflow/unnamed-chunk-3-6.png)

    ## Press [enter] to continue

![](mineRva_Tools_Workflow/unnamed-chunk-3-7.png)

    ## Press [enter] to continue

![](mineRva_Tools_Workflow/unnamed-chunk-3-8.png)

    ## Press [enter] to continue

![](mineRva_Tools_Workflow/unnamed-chunk-3-9.png)

    ## Press [enter] to continue

![](mineRva_Tools_Workflow/unnamed-chunk-3-10.png)

    ## Press [enter] to continue

``` r
related_ion_plots<-ionrPlot(xm2=x_method2,inr=relate_ions,sgrp=sample_group,frac=0.2) 
```

![](mineRva_Tools_Workflow/unnamed-chunk-3-11.png)

    ## Press [enter] to continue

![](mineRva_Tools_Workflow/unnamed-chunk-3-12.png)

    ## Press [enter] to continue

![](mineRva_Tools_Workflow/unnamed-chunk-3-13.png)

    ## Press [enter] to continue

![](mineRva_Tools_Workflow/unnamed-chunk-3-14.png)

    ## Press [enter] to continue

![](mineRva_Tools_Workflow/unnamed-chunk-3-15.png)

    ## Press [enter] to continue

![](mineRva_Tools_Workflow/unnamed-chunk-3-16.png)

    ## Press [enter] to continue

![](mineRva_Tools_Workflow/unnamed-chunk-3-17.png)

    ## Press [enter] to continue

![](mineRva_Tools_Workflow/unnamed-chunk-3-18.png)

    ## Press [enter] to continue

![](mineRva_Tools_Workflow/unnamed-chunk-3-19.png)

    ## Press [enter] to continue

![](mineRva_Tools_Workflow/unnamed-chunk-3-20.png)

    ## Press [enter] to continue

![](mineRva_Tools_Workflow/unnamed-chunk-3-21.png)

    ## Press [enter] to continue

![](mineRva_Tools_Workflow/unnamed-chunk-3-22.png)

    ## Press [enter] to continue

![](mineRva_Tools_Workflow/unnamed-chunk-3-23.png)

    ## Press [enter] to continue

![](mineRva_Tools_Workflow/unnamed-chunk-3-24.png)

    ## Press [enter] to continue

![](mineRva_Tools_Workflow/unnamed-chunk-3-25.png)

    ## Press [enter] to continue

![](mineRva_Tools_Workflow/unnamed-chunk-3-26.png)

    ## Press [enter] to continue

![](mineRva_Tools_Workflow/unnamed-chunk-3-27.png)

    ## Press [enter] to continue

![](mineRva_Tools_Workflow/unnamed-chunk-3-28.png)

    ## Press [enter] to continue

![](mineRva_Tools_Workflow/unnamed-chunk-3-29.png)

    ## Press [enter] to continue

![](mineRva_Tools_Workflow/unnamed-chunk-3-30.png)

    ## Press [enter] to continue

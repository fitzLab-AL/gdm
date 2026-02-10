# Summarize a Fitted Generalized Dissimilarity Model

This function summarizes the gdm model object returned from
[`gdm`](https://mfitzpatrick.al.umces.edu/gdm/reference/gdm.md).

## Usage

``` r
# S3 method for class 'gdm'
summary(object, ...)
```

## Arguments

- object:

  A gdm model object resulting from a call to
  [`gdm`](https://mfitzpatrick.al.umces.edu/gdm/reference/gdm.md).

- ...:

  Ignored.

## Value

summary prints its output to the R Console window and returns no value.

## See also

[`gdm`](https://mfitzpatrick.al.umces.edu/gdm/reference/gdm.md)

## Examples

``` r
##set up site-pair table using the southwest data set
sppData <- southwest[, c(1,2,14,13)]
envTab <- southwest[, c(2:ncol(southwest))]
sitePairTab <- formatsitepair(sppData, 2, XColumn="Long", YColumn="Lat", sppColumn="species",
                              siteColumn="site", predData=envTab)
#> Warning: No abundance column was specified, so the biological data are assumed to be presences.
#> Aggregation function missing: defaulting to length

##create GDM
gdmMod <- gdm(sitePairTab, geo=TRUE)

##summary of GDM
summary(gdmMod)
#> [1] 
#> [1] 
#> [1] GDM Modelling Summary
#> [1] Creation Date:  Tue Feb 10 15:51:25 2026
#> [1] 
#> [1] Name:  gdmMod
#> [1] 
#> [1] Data:  sitePairTab
#> [1] 
#> [1] Samples:  4371
#> [1] 
#> [1] Geographical distance used in model fitting?  TRUE
#> [1] 
#> [1] NULL Deviance:  651.914
#> [1] GDM Deviance:  129.025
#> [1] Percent Deviance Explained:  80.208
#> [1] 
#> [1] Intercept:  0.277
#> [1] 
#> [1] PREDICTOR ORDER BY SUM OF I-SPLINE COEFFICIENTS:
#> [1] 
#> [1] Predictor 1: bio19
#> [1] Splines: 3
#> [1] Min Knot: 114.394
#> [1] 50% Knot: 172.416
#> [1] Max Knot: 554.771
#> [1] Coefficient[1]: 0.941
#> [1] Coefficient[2]: 0.868
#> [1] Coefficient[3]: 0
#> [1] Sum of coefficients for bio19: 1.809
#> [1] 
#> [1] Predictor 2: phTotal
#> [1] Splines: 3
#> [1] Min Knot: 277.978
#> [1] 50% Knot: 584.609
#> [1] Max Knot: 1860.37
#> [1] Coefficient[1]: 1.127
#> [1] Coefficient[2]: 0.23
#> [1] Coefficient[3]: 0
#> [1] Sum of coefficients for phTotal: 1.357
#> [1] 
#> [1] Predictor 3: bio5
#> [1] Splines: 3
#> [1] Min Knot: 25.571
#> [1] 50% Knot: 32.16
#> [1] Max Knot: 36.188
#> [1] Coefficient[1]: 0.127
#> [1] Coefficient[2]: 0.453
#> [1] Coefficient[3]: 0.114
#> [1] Sum of coefficients for bio5: 0.694
#> [1] 
#> [1] Predictor 4: solumDepth
#> [1] Splines: 3
#> [1] Min Knot: 705.02
#> [1] 50% Knot: 1017.628
#> [1] Max Knot: 1247.705
#> [1] Coefficient[1]: 0.682
#> [1] Coefficient[2]: 0
#> [1] Coefficient[3]: 0
#> [1] Sum of coefficients for solumDepth: 0.682
#> [1] 
#> [1] Predictor 5: awcA
#> [1] Splines: 3
#> [1] Min Knot: 12.975
#> [1] 50% Knot: 22.186
#> [1] Max Knot: 50.7
#> [1] Coefficient[1]: 0
#> [1] Coefficient[2]: 0
#> [1] Coefficient[3]: 0.523
#> [1] Sum of coefficients for awcA: 0.523
#> [1] 
#> [1] Predictor 6: Geographic
#> [1] Splines: 3
#> [1] Min Knot: 0.452
#> [1] 50% Knot: 2.46
#> [1] Max Knot: 6.532
#> [1] Coefficient[1]: 0.014
#> [1] Coefficient[2]: 0.372
#> [1] Coefficient[3]: 0
#> [1] Sum of coefficients for Geographic: 0.386
#> [1] 
#> [1] Predictor 7: sandA
#> [1] Splines: 3
#> [1] Min Knot: 56.697
#> [1] 50% Knot: 72.951
#> [1] Max Knot: 83.993
#> [1] Coefficient[1]: 0.092
#> [1] Coefficient[2]: 0
#> [1] Coefficient[3]: 0.139
#> [1] Sum of coefficients for sandA: 0.231
#> [1] 
#> [1] Predictor 8: shcA
#> [1] Splines: 3
#> [1] Min Knot: 78.762
#> [1] 50% Knot: 179.351
#> [1] Max Knot: 521.985
#> [1] Coefficient[1]: 0
#> [1] Coefficient[2]: 0.156
#> [1] Coefficient[3]: 0
#> [1] Sum of coefficients for shcA: 0.156
#> [1] 
#> [1] Predictor 9: bio6
#> [1] Splines: 3
#> [1] Min Knot: 4.373
#> [1] 50% Knot: 5.509
#> [1] Max Knot: 9.224
#> [1] Coefficient[1]: 0.121
#> [1] Coefficient[2]: 0
#> [1] Coefficient[3]: 0
#> [1] Sum of coefficients for bio6: 0.121
#> [1] 
#> [1] Predictor 10: bio15
#> [1] Splines: 3
#> [1] Min Knot: 29.167
#> [1] 50% Knot: 55.008
#> [1] Max Knot: 87.143
#> [1] Coefficient[1]: 0.027
#> [1] Coefficient[2]: 0
#> [1] Coefficient[3]: 0
#> [1] Sum of coefficients for bio15: 0.027
#> [1] 
#> [1] Predictor 11: bio18
#> [1] Splines: 3
#> [1] Min Knot: 0
#> [1] 50% Knot: 0
#> [1] Max Knot: 52
#> [1] Coefficient[1]: 0
#> [1] Coefficient[2]: 0
#> [1] Coefficient[3]: 0
#> [1] Sum of coefficients for bio18: 0
```

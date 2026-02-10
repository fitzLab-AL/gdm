# Fit a Generalized Dissimilarity Model to Tabular Site-Pair Data

The gdm function is used to fit a generalized dissimilarity model to
tabular site-pair data formatted as follows using the
[`formatsitepair`](https://mfitzpatrick.al.umces.edu/gdm/reference/formatsitepair.md)
function: distance, weights, s1.xCoord, s1.yCoord, s2.xCoord, s2.yCoord,
s1.Pred1, s1.Pred2, ...,s1.PredN, s2.Pred1, s2.Pred2, ..., s2.PredN. The
distance column contains the response variable must be any ratio-based
dissimilarity (distance) measure between Site 1 and Site 2. The weights
column defines any weighting to be applied during fitting of the model.
If equal weighting is required, then all entries in this column should
be set to 1.0 (default). The third and fourth columns, s1.xCoord and
s1.yCoord, represent the spatial coordinates of the first site in the
site pair (s1). The fifth and sixth columns, s2.xCoord and s2.yCoord,
represent the coordinates of the second site (s2). Note that the first
six columns are REQUIRED, even if you do not intend to use geographic
distance as a predictor (in which case these columns can be loaded with
dummy data if the actual coordinates are unknown - though that would be
weird, no?). The next N\*2 columns contain values for N predictors for
Site 1, followed by values for the same N predictors for Site 2.  

The following is an example of a GDM input table header with three
environmental predictors (Temp, Rain, Bedrock):  

distance, weights, s1.xCoord, s1.yCoord, s2.xCoord, s2.yCoord, s1.Temp,
s1.Rain, s1.Bedrock, s2.Temp, s2.Rain, s2.Bedrock

## Usage

``` r
gdm(data, geo=FALSE, splines=NULL, knots=NULL)
```

## Arguments

- data:

  A data frame containing the site pairs to be used to fit the GDM
  (obtained using the
  [`formatsitepair`](https://mfitzpatrick.al.umces.edu/gdm/reference/formatsitepair.md)
  function). The observed response data must be located in the first
  column. The weights to be applied to each site pair must be located in
  the second column. If geo is TRUE, then the s1.xCoord, s1.yCoord and
  s2.xCoord, s2.yCoord columns will be used to calculate the geographic
  distance between site pairs for inclusion as the geographic predictor
  term in the model. Site coordinates ideally should be in a projected
  coordinate system (i.e., not longitude-latitude) to ensure proper
  calculation of geographic distances. If geo is FALSE (default), then
  the s1.xCoord, s1.yCoord, s2.xCoord and s2.yCoord data columns must
  still be included, but are ignored in fitting the model. Columns
  containing the predictor data for Site 1, and the predictor data for
  Site 2, follow.

- geo:

  Set to TRUE if geographic distance between sites is to be included as
  a model term. Set to FALSE if geographic distance is to be omitted
  from the model. Default is FALSE.

- splines:

  An optional vector of the number of I-spline basis functions to be
  used for each predictor in fitting the model. If supplied, it must
  have the same length as the number of predictors (including geographic
  distance if geo is TRUE). If this vector is not provided
  (splines=NULL), then a default of 3 basis functions is used for all
  predictors.

- knots:

  An optional vector of knots in *units of the predictor variables* to
  be used in the fitting process. If knots are supplied and
  splines=NULL, then the knots argument must have the same length as the
  number of predictors \* n, where n is the number of knots (default=3).
  If both knots and the number of splines are supplied, then the length
  of the knots argument must be the same as the sum of the values in the
  splines vector. Note that the default values for knots when the
  default three I-spline basis functions are 0 (minimum), 50 (median),
  and 100 (maximum) quantiles.

## Value

gdm returns a gdm model object. The function
[`summary.gdm`](https://mfitzpatrick.al.umces.edu/gdm/reference/summary.gdm.md)
can be used to obtain or print a synopsis of the results. A gdm model
object is a list containing at least the following components:

- dataname:

  The name of the table used as the data argument to the model.

- geo:

  Whether geographic distance was used as a predictor in the model.

- gdmdeviance:

  The deviance of the fitted GDM model.

- nulldeviance:

  The deviance of the null model.

- explained:

  The percentage of null deviance explained by the fitted GDM model.

- intercept:

  The fitted value for the intercept term in the model.

- predictors:

  A list of the names of the predictors that were used to fit the model,
  in order of the amount of turnover associated with each predictor
  (based on the sum of the I-spline coefficients).

- coefficients:

  A list of the coefficients for each spline for each of the predictors
  considered in model fitting.

- knots:

  A vector of the knots derived from the x data (or user defined), for
  each predictor.

- splines:

  A vector of the number of I-spline basis functions used for each
  predictor.

- creationdate:

  The date and time of model creation.

- observed:

  The observed response for each site pair (from data column 1).

- predicted:

  The predicted response for each site pair, from the fitted model
  (after applying the link function).

- ecological:

  The linear predictor (ecological distance) for each site pair, from
  the fitted model (before applying the link function).

## References

Ferrier S, Manion G, Elith J, Richardson, K (2007) Using generalized
dissimilarity modelling to analyse and predict patterns of beta
diversity in regional biodiversity assessment. *Diversity &
Distributions* 13, 252-264.

## See also

[`formatsitepair`](https://mfitzpatrick.al.umces.edu/gdm/reference/formatsitepair.md)`, `[`summary.gdm`](https://mfitzpatrick.al.umces.edu/gdm/reference/summary.gdm.md)`, `[`plot.gdm`](https://mfitzpatrick.al.umces.edu/gdm/reference/plot.gdm.md)`, `[`predict.gdm`](https://mfitzpatrick.al.umces.edu/gdm/reference/predict.gdm.md)`, `[`gdm.transform`](https://mfitzpatrick.al.umces.edu/gdm/reference/gdm.transform.md)

## Examples

``` r
 ##fit table environmental data
 # format site-pair table using the southwest data table
 head(southwest)
#>   species site    awcA  phTotal   sandA     shcA solumDepth     bio5     bio6
#> 1    spp1 1066 14.4725 546.1800 71.3250 178.8650   875.1725 31.43824 5.058823
#> 2    spp1 1026 16.2575 470.9950 68.8975 105.8400   928.4925 33.14412 4.852941
#> 3    spp1 1025 23.1375 459.7425 71.4700  88.3550   892.2275 32.84000 4.817143
#> 4    spp1 1026 16.2575 470.9950 68.8975 105.8400   928.4925 33.14412 4.852941
#> 5    spp1 1027 17.0175 489.3950 74.6775 147.2125   951.9050 33.17813 4.590625
#> 6    spp1 1047 17.3625 515.0825 75.7525 164.1875   981.4750 32.61579 4.676316
#>      bio15 bio18    bio19       Lat     Long
#> 1 40.38235     0 132.6471 -32.99425 118.7573
#> 2 48.20588     0 140.2941 -32.04285 118.3495
#> 3 53.88571    43 145.0571 -31.99067 117.8260
#> 4 48.20588     0 140.2941 -32.04285 118.3495
#> 5 44.00000     0 135.6875 -32.09326 118.8736
#> 6 42.00000     0 134.0263 -32.54354 118.8157
 sppData <- southwest[, c(1,2,13,14)]
 envTab <- southwest[, c(2:ncol(southwest))]

 sitePairTab <- formatsitepair(sppData, 2, XColumn="Long", YColumn="Lat", sppColumn="species",
                               siteColumn="site", predData=envTab)
#> Warning: No abundance column was specified, so the biological data are assumed to be presences.
#> Aggregation function missing: defaulting to length

 ##fit table GDM
 gdmTabMod <- gdm(sitePairTab, geo=TRUE)
 summary(gdmTabMod)
#> [1] 
#> [1] 
#> [1] GDM Modelling Summary
#> [1] Creation Date:  Tue Feb 10 15:51:13 2026
#> [1] 
#> [1] Name:  gdmTabMod
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

 ##fit raster environmental data
 ##sets up site-pair table
 rastFile <- system.file("./extdata/swBioclims.grd", package="gdm")
 envRast <- terra::rast(rastFile)

 ##environmental raster data
 sitePairRast <- formatsitepair(sppData, 2, XColumn="Long",
                                YColumn="Lat", sppColumn="species",
                                siteColumn="site", predData=envRast)
#> Warning: No abundance column was specified, so the biological data are assumed to be presences.
#> Aggregation function missing: defaulting to length
#> Warning: When using rasters for environmental covariates (predictors), each site is assigned to the
#>               raster cell in which the site is located. If more than one site occurs within the same raster cell,
#>               the biological data of those sites are aggregated (more likely as raster resolution decreases).
 ##sometimes raster data returns NA in the site-pair table, these rows will
 ##have to be removed before fitting gdm
 sitePairRast <- na.omit(sitePairRast)

 ##fit raster GDM
 gdmRastMod <- gdm(sitePairRast, geo=TRUE)
 summary(gdmRastMod)
#> [1] 
#> [1] 
#> [1] GDM Modelling Summary
#> [1] Creation Date:  Tue Feb 10 15:51:14 2026
#> [1] 
#> [1] Name:  gdmRastMod
#> [1] 
#> [1] Data:  sitePairRast
#> [1] 
#> [1] Samples:  4278
#> [1] 
#> [1] Geographical distance used in model fitting?  TRUE
#> [1] 
#> [1] NULL Deviance:  634.418
#> [1] GDM Deviance:  190.801
#> [1] Percent Deviance Explained:  69.925
#> [1] 
#> [1] Intercept:  0.395
#> [1] 
#> [1] PREDICTOR ORDER BY SUM OF I-SPLINE COEFFICIENTS:
#> [1] 
#> [1] Predictor 1: bio19
#> [1] Splines: 3
#> [1] Min Knot: 112
#> [1] 50% Knot: 166
#> [1] Max Knot: 621
#> [1] Coefficient[1]: 1.416
#> [1] Coefficient[2]: 0.869
#> [1] Coefficient[3]: 0
#> [1] Sum of coefficients for bio19: 2.285
#> [1] 
#> [1] Predictor 2: bio5
#> [1] Splines: 3
#> [1] Min Knot: 257
#> [1] 50% Knot: 324
#> [1] Max Knot: 362
#> [1] Coefficient[1]: 0.025
#> [1] Coefficient[2]: 0.514
#> [1] Coefficient[3]: 0.39
#> [1] Sum of coefficients for bio5: 0.929
#> [1] 
#> [1] Predictor 3: Geographic
#> [1] Splines: 3
#> [1] Min Knot: 0.419
#> [1] 50% Knot: 2.451
#> [1] Max Knot: 6.511
#> [1] Coefficient[1]: 0.134
#> [1] Coefficient[2]: 0.523
#> [1] Coefficient[3]: 0
#> [1] Sum of coefficients for Geographic: 0.657
#> [1] 
#> [1] Predictor 4: bio15
#> [1] Splines: 3
#> [1] Min Knot: 28
#> [1] 50% Knot: 54
#> [1] Max Knot: 88
#> [1] Coefficient[1]: 0.191
#> [1] Coefficient[2]: 0
#> [1] Coefficient[3]: 0
#> [1] Sum of coefficients for bio15: 0.191
#> [1] 
#> [1] Predictor 5: bio6
#> [1] Splines: 3
#> [1] Min Knot: 41
#> [1] 50% Knot: 55
#> [1] Max Knot: 95
#> [1] Coefficient[1]: 0.07
#> [1] Coefficient[2]: 0
#> [1] Coefficient[3]: 0.013
#> [1] Sum of coefficients for bio6: 0.083
#> [1] 
#> [1] Predictor 6: bio18
#> [1] Splines: 3
#> [1] Min Knot: 29
#> [1] 50% Knot: 48
#> [1] Max Knot: 89
#> [1] Coefficient[1]: 0
#> [1] Coefficient[2]: 0
#> [1] Coefficient[3]: 0
#> [1] Sum of coefficients for bio18: 0
```

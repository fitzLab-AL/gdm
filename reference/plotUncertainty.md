# Plot I-splines With Error Bands Using Bootstrapping.

This function estimates uncertainty in the fitted I-splines by fitting
many GDMs using a subsample of the data. The function can run in
parallel on multicore machines to reduce computation time (recommended
for large number of iterations). I-spline plots with error bands (+/-
one standard deviation) are produced showing (1) the variance of
I-spline coefficients and (2) a rug plot indicating how sites used in
model fitting are distributed along each gradient. Function result
optionally can be saved to disk as a csv for custom plotting, etc. The
result output table will have 6 columns per predictor, three each for
the x and y values containing the lower bound, full model, and upper
bound.

## Usage

``` r
plotUncertainty(spTable, sampleSites, bsIters, geo=FALSE,
splines=NULL, knots=NULL, splineCol="blue", errCol="grey80",
plot.linewidth=2.0, plot.layout=c(2,2), parallel=FALSE, cores=2, save=FALSE,
fileName="gdm.plotUncertainy.csv")
```

## Arguments

- spTable:

  A site-pair table, same as used to fit a
  [`gdm`](https://mfitzpatrick.al.umces.edu/gdm/reference/gdm.md).

- sampleSites:

  The fraction (0-1) of sites to retain from the full site-pair table
  when subsampling.

- bsIters:

  The number of bootstrap iterations to perform.

- geo:

  Same as the
  [`gdm`](https://mfitzpatrick.al.umces.edu/gdm/reference/gdm.md) geo
  argument.

- splines:

  Same as the
  [`gdm`](https://mfitzpatrick.al.umces.edu/gdm/reference/gdm.md)
  splines argument.

- knots:

  Same as the
  [`gdm`](https://mfitzpatrick.al.umces.edu/gdm/reference/gdm.md) knots
  argument.

- splineCol:

  The color of the plotted mean spline. The default is "blue".

- errCol:

  The color of shading for the error bands (+/- one standard deviation
  around the mean line). The default is "grey80".

- plot.linewidth:

  The line width of the plotted mean spline line. The default is 2.

- plot.layout:

  Same as the
  [`plot.gdm`](https://mfitzpatrick.al.umces.edu/gdm/reference/plot.gdm.md)
  plot.layout argument.

- parallel:

  Perform the uncertainty assessment using multiple cores? Default =
  FALSE.

- cores:

  When the parallel argument is set to TRUE, the number of cores to be
  registered for the foreach loop. Must be \<= the number of cores in
  the machine running the function.

- save:

  Save the function result (e.g., for custom plotting)? Default=FALSE.

- fileName:

  Name of the csv file to save the data frame that contains the function
  result. Default = gdm.plotUncertainy.csv. Ignored if save=FALSE.

## Value

plotUncertainty returns NULL. Saves a csv to disk if save=TRUE.

## References

Shryock, D. F., C. A. Havrilla, L. A. DeFalco, T. C. Esque, N. A.
Custer, and T. E. Wood. 2015. Landscape genomics of *Sphaeralcea
ambigua* in the Mojave Desert: a multivariate, spatially-explicit
approach to guide ecological restoration. *Conservation Genetics*
16:1303-1317.

## See also

[`plot.gdm`](https://mfitzpatrick.al.umces.edu/gdm/reference/plot.gdm.md)`, `[`formatsitepair`](https://mfitzpatrick.al.umces.edu/gdm/reference/formatsitepair.md)`, `[`subsample.sitepair`](https://mfitzpatrick.al.umces.edu/gdm/reference/subsample.sitepair.md)

## Examples

``` r
##set up site-pair table using the southwest data set
sppData <- southwest[, c(1,2,13,14)]
envTab <- southwest[, c(2:ncol(southwest))]
sitePairTab <- formatsitepair(sppData, 2, XColumn="Long", YColumn="Lat",
                              sppColumn="species", siteColumn="site", predData=envTab)
#> Warning: No abundance column was specified, so the biological data are assumed to be presences.
#> Aggregation function missing: defaulting to length

##plot GDM uncertainty using one core
#not run
#plotUncertainty(sitePairTab, sampleSites=0.70, bsIters=5, geo=TRUE, plot.layout=c(3,3))

##plot GDM uncertainty in parallel
#not run
#plotUncertainty(sitePairTab, sampleSites=0.70, bsIters=50, geo=TRUE, plot.layout=c(3,3),
                 #parallel=T, cores=10)
```

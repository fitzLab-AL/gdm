# Overview of the functions in the gdm package

Generalized Dissimilarity Modeling is a statistical technique for
modelling variation in biodiversity between pairs of geographical
locations or through time. The gdm package provides functions to fit,
evaluate, summarize, and plot Generalized Dissimilarity Models and to
make predictions (across space and/or through time) and map biological
patterns by transforming environmental predictor variables.

## Details

The functions in the gdm package provide the tools necessary for fitting
GDMs, including functions to prepare biodiversity and environmental
data. Major functionality includes:

- Formatting various types of biodiversity and environmental data to
  gdm's site-pair format used in model fitting

- Fitting GDMs using geographic and environmental distances between
  sites

- Plotting fitted functions & extracting I-spline values

- Estimating predictor importance using matrix permutation and predictor
  contributions using deviance paritioning

- Using cross-validation to evaluate models

- Predicting pairwise dissimiliarites between sites or times and
  transforming environmental predictors to biological importance and
  mapping these patterns.

To see the preferable citation of the package, type `citation("gdm")`.

## I. Formatting input data

GDM fits biological distances to pairwise site geographical and
environmental distances. Most users will need to first format their data
to gdm's site-pair table format:

|                                                                                       |                                                                    |
|---------------------------------------------------------------------------------------|--------------------------------------------------------------------|
| [`formatsitepair`](https://mfitzpatrick.al.umces.edu/gdm/reference/formatsitepair.md) | To convert biodiversity and environmental data to site-pair format |

## II. Model fitting, evaluation, and summary

|                                                                                                       |                                                        |
|-------------------------------------------------------------------------------------------------------|--------------------------------------------------------|
| [`gdm`](https://mfitzpatrick.al.umces.edu/gdm/reference/gdm.md)                                       | To fit a GDM                                           |
| [`gdm.crossvalidation`](https://mfitzpatrick.al.umces.edu/gdm/reference/gdm.crossvalidation.md)       | To evaluate a GDM                                      |
| [`gdm.partition.deviance`](https://mfitzpatrick.al.umces.edu/gdm/reference/gdm.partition.deviance.md) | To asses predictor contributions to deviance explained |
| [`gdm.varImp`](https://mfitzpatrick.al.umces.edu/gdm/reference/gdm.varImp.md)                         | To asses model significance and predictor importance   |
| [`summary`](https://rdrr.io/r/base/summary.html)                                                      | To summarize a GDM                                     |

## III. Model prediction and transformation of environmental data

|                                                                                     |                                                                                      |
|-------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------|
| [`predict`](https://rdrr.io/r/stats/predict.html)                                   | To predict biological dissimilarities between sites in space or between time periods |
| [`gdm.transform`](https://mfitzpatrick.al.umces.edu/gdm/reference/gdm.transform.md) | To transform each environmental predictor to biological importance                   |

## IV. Plotting model output and fitted functions

|                                                                                         |                                                            |
|-----------------------------------------------------------------------------------------|------------------------------------------------------------|
| [`plot`](https://rdrr.io/r/graphics/plot.default.html)                                  | To plot model fit and I-splines                            |
| [`isplineExtract`](https://mfitzpatrick.al.umces.edu/gdm/reference/isplineExtract.md)   | To extract I-spline values to allow for custom plotting    |
| [`plotUncertainty`](https://mfitzpatrick.al.umces.edu/gdm/reference/plotUncertainty.md) | To estimate and plot model sensitivity using bootstrapping |

## Author

The gdm development team is Matt Fitzpatrick and Karel Mokany. The R
package is based on code originally developed by Glenn Manion under the
direction of Simon Ferrier. Where others have contributed to individual
functions, credits are provided in function help pages.

The maintainer of the R version of gdm is Matt Fitzpatrick
\<mfitzpatrick@umces.edu\>.

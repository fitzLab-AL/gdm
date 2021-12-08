<!-- See http://style.tidyverse.org/news.html for advice on writing news -->

# gdm 1.5.0-2
* Added *verbose* argument to the `formatsitepair` function to control whether site-pair attributes are printed.

### bug fixes
* Fix bug that prevented printing of I-spline coefficients when using the `summary` function.
* Fix bug that prevented the *sampleSites* argument in the `sitepairformat` function from actually subsampling the site-pair table.

# gdm 1.5.0-1
### bug fixes
* Fix bug associated with the formatting of the `gdm` model object that broke the `predict.gdm` function.

# gdm 1.5.0

* Substantial update to the package structure to use `devtools`, `roxygen2`, etc. development tools.

* Fixed a number of bugs associated with using non-default numebrs of splines. 

* Add two new functions `gdm.partition.deviance` and `gdm.crossvalidation`. The `gdm.partition.deviance` functions partitions deviance explained by different sets of predictors into unique/shared components. The `gdm.crossvalidation` function produces a number of model evaluation metrics using bootstrapping.

* `summary` now provides the sum of the I-spline coefficients (indicator of predictor importance) and lists variables from highest sum to lowest.

* `plot` now plots I-splines in order from highest sum of coefficients to least.

<!-- See http://style.tidyverse.org/news.html for advice on writing news -->

# gdm 1.6.0-2
### bug fixes
* Fix negative deviance explained values returned by `gdm.crossvalidation`.  

### Updates
* Address new CRAN C++11 requirements.

# gdm 1.6.0
* Transition to `terra` functionality for much faster processing of rasters.  
* Fix bug in the `formatsitepair` function for bio data provided as distance matrix.

# gdm 1.5.0-9
### bug fixes
* Fix bug in the `varImp` function related to NULL models

# gdm 1.5.0-8
### bug fixes
* Fix bug in the `formatsitepair` function when using custom weights. 

# gdm 1.5.0-7
### bug fixes
* Fix bug in the `varImp` function when predSelect = F.

# gdm 1.5.0-6
### bug fixes
* Fix bug in the `varImp` function related to testing models with few (<4) variables.

# gdm 1.5.0-5
### bug fixes
* Fix bug in the `varImp` function if sum of I-spline coeffs = 0 for geographic distance variable.

# gdm 1.5.0-4
### bug fixes
* Fix minor bug in the `varImp` function.

# gdm 1.5.0-3
* New, much faster version of the `varImp` function.

### bug fixes
* Fix bug that broke `varImp` if using non-default number of I-splines.

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

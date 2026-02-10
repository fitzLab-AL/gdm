# Calculate GDM Deviance for Observed & Predicted Dissimilarities

Calculate GDM deviance for observed & predicted dissimilarities. Can be
used for assessing cross-validation data. Translated from the c++
function CalcGDMDevianceDouble() in the file NNLS_Double.cpp from the
GDM R package.

## Usage

``` r
calculate.gdm.deviance(predDiss, obsDiss)
```

## Arguments

- predDiss:

  (float) A vector of predicted dissimilarity values, of same length as
  obsDiss.

- obsDiss:

  (float) A vector of observed dissimilarity values, of same length as
  predDiss.

## Value

A single value (float) being the deviance.

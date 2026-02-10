# Create Site-Pair Table

Creates a site-pair table from the lower half of a site-by-site distance
(dissimilarity) matrix. This function is called from the
[`formatsitepair`](https://mfitzpatrick.al.umces.edu/gdm/reference/formatsitepair.md)
function and not needed by the user.

## Usage

``` r
createsitepair(dist, spdata, envInfo, dXCol, dYCol, siteCol, weightsType,
custWeights)
```

## Arguments

- dist:

  The lower half of a site-by-site distance (dissimilarity) matrix,
  provided by the
  [`formatsitepair`](https://mfitzpatrick.al.umces.edu/gdm/reference/formatsitepair.md)
  function.

- spdata:

  Input species data, the same as the bioData input to the
  [`formatsitepair`](https://mfitzpatrick.al.umces.edu/gdm/reference/formatsitepair.md)
  function.

- envInfo:

  Input environmental data. Only accepts data tables as input. If the
  environmental data for
  [`formatsitepair`](https://mfitzpatrick.al.umces.edu/gdm/reference/formatsitepair.md)
  are rasters, the data would have been extracted into table format
  within
  [`formatsitepair`](https://mfitzpatrick.al.umces.edu/gdm/reference/formatsitepair.md).

- dXCol:

  Input x coordinate, the same as the XColumn input to the
  [`formatsitepair`](https://mfitzpatrick.al.umces.edu/gdm/reference/formatsitepair.md)
  function.

- dYCol:

  Input y coordinate, the same as the YColumn input to the
  [`formatsitepair`](https://mfitzpatrick.al.umces.edu/gdm/reference/formatsitepair.md)
  function.

- siteCol:

  Site column, taken from either the species or environmental tables.

- weightsType:

  The method of determining the site-pair weights used in model fitting.

- custWeights:

  Custom weights, as a vector, if given by the user.

## Value

A site-pair table with appropriate distance (dissimilarity) and weight
columns used for fitting GDM.

## Note

This function is called from the
[`formatsitepair`](https://mfitzpatrick.al.umces.edu/gdm/reference/formatsitepair.md)
function and not needed by the user.

## See also

[`formatsitepair`](https://mfitzpatrick.al.umces.edu/gdm/reference/formatsitepair.md)

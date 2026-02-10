# Permutate Site-Pair Table Rows, Internal Function

A function which randomizes the rows of a given site-pair table. This
function is called from the
[`gdm.varImp`](https://mfitzpatrick.al.umces.edu/gdm/reference/gdm.varImp.md)
function and not needed by the user. be called directly by the user.

## Usage

``` r
permutateSitePair(spTab, siteVarTab, indexTab, vNames)
```

## Arguments

- spTab:

  A site-pair table.

- siteVarTab:

  A site x variable table.

- indexTab:

  A table of index values for the site-pair table.

- vNames:

  Vector of variable names in both the site-pair and site x variable
  tables.

## Value

A new site-pair table with rows randomized.

## Note

This function is called from the
[`gdm.varImp`](https://mfitzpatrick.al.umces.edu/gdm/reference/gdm.varImp.md)
function and not needed by the user.

## See also

[`gdm.varImp`](https://mfitzpatrick.al.umces.edu/gdm/reference/gdm.varImp.md)

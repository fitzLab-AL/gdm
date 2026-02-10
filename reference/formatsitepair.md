# Combines Biological and Environmental Data to Produce a GDM-formatted Site-Pair Table

This function takes input biological data and environmental, geographic,
and other predictor data and builds a site-pair table required for
fitting a Generalized Dissimilarity Model using the
[`gdm`](https://mfitzpatrick.al.umces.edu/gdm/reference/gdm.md)
function. NOTE: x-y coordinates of sites MUST be present in either the
biological or the environmental data. Site coordinates ideally should be
in a projected coordinate system (i.e., not longitude-latitude) to
ensure proper calculation of geographic distances.  

The input biological data can be in one of the following four formats.
Note that the general term "species" is used, but any classification of
biological entities (e.g. functional types, haplotypes, etc) can be used
as long as an appropriate distance metric is also supplied (see "dist"
argument):  

1.  site-by-species matrix  

2.  x, y, species list  

3.  site-by-site biological distance (dissimilarity) matrix  

4.  an existing site-pair table (see Details)  

Predictor data can be provided in three formats:  

1.  a site-by-predictor matrix with a column for each predictor variable
    and a row for each site  

2.  a terra object SpatRaster, with one raster for each predictor
    variable  

3.  one or more site-by-site distance matrices using the "distPreds"
    argument (see below).

## Usage

``` r
formatsitepair(bioData, bioFormat, dist="bray", abundance=FALSE, siteColumn=NULL,
XColumn, YColumn, sppColumn=NULL, abundColumn=NULL, sppFilter=0, predData,
distPreds=NULL, weightType="equal", custWeights=NULL, sampleSites=1, verbose=FALSE)
```

## Arguments

- bioData:

  The input biological (the response variable) data table, in one of the
  four formats defined above (see Details).

- bioFormat:

  An integer code specifying the format of bioData. Acceptable values
  are 1, 2, 3, or 4 (see Details).

- dist:

  Default = "bray". A character code indicating the metric to quantify
  pairwise site distances / dissimilarities. Calls the
  [`vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.html)
  function from the vegan package to calculate dissimilarity and
  therefore accepts any method available from
  [`vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.html).

- abundance:

  Default = FALSE. Indicates whether the biological data are abundance
  data (TRUE) or presence-absence (0, 1) data (FALSE).

- siteColumn:

  The name of the column in either the biological or environmental data
  table containing a unique site identifier. If a site column is
  provided in both the biological and environmental data, the site
  column name must be the same in both tables.

- XColumn:

  The name of the column containing x-coordinates of sites.
  X-coordinates can be provided in either the biological or
  environmental data tables, but MUST be in at least one of them. If an
  x-coordinate column is provided in both the biological and
  environmental data, the column name must be identical. Site
  coordinates ideally should be in a projected coordinate system (i.e.,
  not longitude-latitude) to ensure proper calculation of geographic
  distances. Note that if you are using rasters, they must be in the
  same coordinate system as the site coordinates.

- YColumn:

  The name of the column containing y-coordinates of sample sites.
  Y-coordinates can be provided in either the biological or
  environmental data tables, but MUST be in at least one of them. If a
  y-coordinate column is provided in both the biological and
  environmental data, the column name must be identical. Site
  coordinates ideally should be in a projected coordinate system (i.e.,
  not longitude-latitude) to ensure proper calculation of geographic
  distances. Note that if you are using rasters, they must be in the
  same coordinate system as the site coordinates.

- sppColumn:

  Only used if bioFormat = 2 (x, y, species list). The name of the
  column containing unique name / identifier for each species.

- abundColumn:

  If abundance = TRUE, this parameter identifies the column containing
  the measure of abundance at each site. Only used if bioFormat = 2
  (i.e., x, y, species list), though in the case of abundance data, the
  table would have four columns: x, y, species, abundance.

- sppFilter:

  Default = 0. To account for limited sampling effort at some sites,
  sppFilter removes all sites at which the number of recorded species
  (i.e., observed species richness) is less than the specified value.
  For example, if sppFilter = 5, all sites with fewer than 5 recorded
  species will be removed.

- predData:

  The environmental predictor data. Accepts either a site-by-predictor
  table or a terra object SpatRaster.

- distPreds:

  An optional list of distance matrices to be used as predictors in
  combination with predData. For example, a site-by-site dissimilarity
  matrix for one biological group (e.g., trees) can be used as a
  predictor for another group (e.g., ferns). Each distance matrix must
  have as the first column the names of the sites (therefore the matrix
  will not be square). The name of the column containing the site names
  should have the same name as that provided for the siteColumn
  argument. Site IDs are required here to ensure correct ordering of
  sites in the construction of the site-pair table. Note that the
  formatsitepair function will not accept distance matrices in the as
  the only predictors (i.e., at least one additional, non-distPreds
  predictor variable is required). If you wish to fit GDM using only
  distance matrices provided using distPreds, provide one fake predictor
  (e.g., with all sites having the same value), plus site and coordinate
  columns. The s1 and s2 columns for this fake variable can then be
  removed by hand before fitting the GDM.

- weightType:

  Default = "equal". Defines the weighting for sites. Can be either: (1)
  "equal" (weights for all sites set = 1), (2) "richness" (each site
  weighted according to number of species recorded), or (3) "custom"
  (user defined). If weightType="custom", the user must provide a vector
  of site weights equal to the number of rows in the full site-pair
  table (i.e., before species filtering (sppFilter argument) or
  sub-sampling is taken into account (sampleSites argument)).

- custWeights:

  A two column matrix or data frame of user-defined site weights. The
  first column should be the site name and should be named the same as
  that provided for the siteColumn argument. The second column should be
  numeric weight values and should be named "weights". The weight values
  represent the importance of each site in model fitting, and the values
  in the output site-pair table is an average of the two sites in each
  site-pair. Required when weightType = "custom". Ignored otherwise.

- sampleSites:

  Default = 1. A number between 0-1 indicating the fraction of sites to
  be used to construct the site-pair table. This argument can be used to
  reduce the number of sites to overcome possible memory limitations
  when fitting models with very large numbers of sites.

- verbose:

  Default = FALSE. If TRUE, summary of information regarding dimensions
  of the site-pair table will be printed that can be useful for
  diagnostics.

## Value

A formatted site-pair table containing the response (biological distance
or dissimilarity), predictors, and weights as required for fitting
Generalized Dissimilarity Models.

## Details

bioData and bioFormat: The function accepts biological data in the
following formats:

bioData = site-by-species matrix; bioFormat = 1: assumes that the
response data are provided with a site ID column (specified by siteCol)
and, optionally, two columns for the x & y coordinates of the sites. All
remaining columns contain the biological data, with a column for each
biological entity (most commonly species). In the case that a raster
stack (a terra object SpatRaster) is provided for the environmental data
(predData), x-y coordinates MUST be provided in bioData to allow
extraction of the environmental data at site locations. The x-y
coordinates will be intersected with the raster stack and, if the number
of unique cells intersected by the points is less than the number of
unique site IDs (i.e. multiple sites fall within a single cell), the
function will use the raster cell as the site ID and aggregate sites
accordingly. Therefore, model fitting will be sensitive to raster cell
size. If the environmental data are in tabular format, they should have
the same number of sites (i.e., same number of rows) as bioData. The x-y
coordinate and site ID columns must have the same names in bioData and
predData.

bioData = x, y, species list (optionally a fourth column with abundance
can be provided); bioFormat = 2: assumes a table of 3 or 4 columns, the
first two being the x & y coordinates of species records, the third
(sppCol) being the name / identifier of the species observed at that
location, and optionally a fourth column indicating a measure of
abundance. If an abundance column is not provided, presence-only data
are assumed. In the case that a raster stack (a terra object SpatRaster)
is provided for the environmental data (predData), the x-y coordinates
will be intersected with the raster stack and, if the number of unique
cells intersected by the points is less than the number of unique site
IDs (i.e. multiple sites fall within a single cell), the function will
use the raster cell as the site ID and aggregate sites accordingly.
Therefore, model fitting will be sensitive to raster cell size.

bioData = site-site distance (dissimilarity) matrix; bioFormat = 3. This
option allows the use of an existing site-site distance (dissimilarity)
matrix, such as genetic distance matrix calculated outside of the gdm
package. Only the lower triangle of the matrix is required to create the
site-pair table, but the function will automatically removes the upper
triangle if present. The code checks and aligns the order of sites in
the distance matrix and the predictor data to ensure they match. To do
so, (1) a site column is required in both the distance matrix and the
predictor data and (2) site IDs are required to be a number. This is the
only bioFormat in which the environmental data MAY NOT be a raster
stack.

bioData = site-pair table; bioFormat = 4: with an already created
site-pair table, this option allows the user to add one or more distance
matrices (see distPreds above) to the existing site-pair table and/or
sub-sample the site-pair table (see sample above). If the site-pair
table was not created using the formatsitepair function, the user will
need to ensure the order of the sites matches that in other tables being
provided to the function.

NOTES: (1) The function assumes that the x-y coordinates and the raster
stack (if used) are in the same coordinate system. No checking is
performed to confirm this is the case. (2) The function assumes that the
association between the provided site and x-y coordinate columns are
singular and unique. Therefore, the function will fail should a given
site has more than one sets of coordinates associated with it, as well
as multiple sites being given the exact same coordinates.

## Examples

``` r
## tabular data
# start with the southwest data table
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

#########table type 1
## site-species table without coordinates
testData1a <- reshape2::dcast(sppData, site~species)
#> Using Long as value column: use value.var to override.
#> Aggregation function missing: defaulting to length
##site-species table with coordinates
coords <- unique(sppData[, 2:ncol(sppData)])
testData1b <- merge(testData1a, coords, by="site")
## site-species, table-table
exFormat1a <- formatsitepair(testData1a, 1, siteColumn="site", XColumn="Long",
YColumn="Lat", predData=envTab)

#' # next, let's try environmental raster data
## not run
# rastFile <- system.file("./extdata/swBioclims.grd", package="gdm")
# envRast <- terra::rast(rastFile)

## site-species, table-raster
## not run
# exFormat1b <- formatsitepair(testData1b, 1, siteColumn="site", XColumn="Long",
# YColumn="Lat", predData=envRast)

#########table type 2
## site xy spp list, table-table
exFormat2a <- formatsitepair(sppData, 2, XColumn="Long", YColumn="Lat",
sppColumn="species", siteColumn="site", predData=envTab)
#> Warning: No abundance column was specified, so the biological data are assumed to be presences.
#> Aggregation function missing: defaulting to length

## site xy spp list, table-raster
## not run
# exFormat2b <- formatsitepair(sppData, 2, XColumn="Long", YColumn="Lat",
# sppColumn="species", siteColumn="site", predData=envRast)

#########table type 3
## It is possible to format a site-pair table by starting
# with a pre-calculated matrix of biological distances
dim(gdmDissim) # pairwise distance matrix + 1 column for site IDs
#> [1] 94 95
gdmDissim[1:5, 1:5]
#>    site         1         2         3         4
#> V1  881 0.0000000 0.4485981 0.7575758 0.8939394
#> V2  882 0.4485981 0.0000000 0.5837563 0.8170732
#> V3  883 0.7575758 0.5837563 0.0000000 0.4782609
#> V4  884 0.8939394 0.8170732 0.4782609 0.0000000
#> V5  885 0.9178082 0.8202247 0.5813953 0.4375000
# now we can format the table:
exFormat3 <- formatsitepair(gdmDissim, 3, XColumn="Long", YColumn="Lat",
                            predData=envTab, siteColumn="site")

#########table type 4
## adds a predictor matrix to an existing site-pair table, in this case,
## predData needs to be provided, but is not actually used
exFormat4 <- formatsitepair(exFormat2a, 4, predData=envTab, siteColumn="site",
                            distPreds=list(as.matrix(gdmDissim)))
```

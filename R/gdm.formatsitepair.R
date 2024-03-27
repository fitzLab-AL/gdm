#' @title Combines Biological and Environmental Data to Produce a GDM-formatted
#' Site-Pair Table
#'
#' @description This function takes input biological data and environmental,
#' geographic, and other predictor data and builds a site-pair table required
#' for fitting a Generalized Dissimilarity Model using the \code{\link[gdm]{gdm}}
#' function. NOTE: x-y coordinates of sites MUST be present in either the
#' biological or the environmental data. Site coordinates ideally should be in a
#' projected coordinate system (i.e., not longitude-latitude) to ensure proper
#' calculation of geographic distances. \cr
#'
#' The input biological data can be in one of the following four formats. Note
#' that the general term "species" is used, but any classification of biological
#' entities (e.g. functional types, haplotypes, etc) can be used as long as an
#' appropriate distance metric is also supplied (see "dist" argument): \cr
#'
#' \enumerate{
#'   \item site-by-species matrix \cr
#'
#'   \item x, y, species list \cr
#'
#'   \item site-by-site biological distance (dissimilarity) matrix \cr
#'
#'   \item an existing site-pair table (see Details) \cr
#'   }
#'
#' Predictor data can be provided in three formats: \cr
#'
#' \enumerate{
#'   \item a site-by-predictor matrix with a column for each predictor variable
#'   and a row for each site \cr
#'
#'   \item a terra object SpatRaster, with one raster for each predictor variable \cr
#'
#'   \item one or more site-by-site distance matrices using the "distPreds"
#'   argument (see below).
#'   }
#'
#' @usage
#' formatsitepair(bioData, bioFormat, dist="bray", abundance=FALSE, siteColumn=NULL,
#' XColumn, YColumn, sppColumn=NULL, abundColumn=NULL, sppFilter=0, predData,
#' distPreds=NULL, weightType="equal", custWeights=NULL, sampleSites=1, verbose=FALSE)
#'
#' @param bioData The input biological (the response variable) data table, in
#' one of the four formats defined above (see Details).
#'
#' @param bioFormat An integer code specifying the format of bioData. Acceptable
#' values are 1, 2, 3, or 4 (see Details).
#'
#' @param dist Default = "bray". A character code indicating the metric to
#' quantify pairwise site distances / dissimilarities. Calls the
#' \code{\link[vegan]{vegdist}} function from the \pkg{vegan} package to
#' calculate dissimilarity and therefore accepts any method available from
#' \code{\link[vegan]{vegdist}}.
#'
#' @param abundance Default = FALSE. Indicates whether the biological data are
#' abundance data (TRUE) or presence-absence (0, 1) data (FALSE).
#'
#' @param siteColumn The name of the column in either the biological or
#' environmental data table containing a unique site identifier. If a site column
#' is provided in both the biological and environmental data, the site column name
#' must be the same in both tables.
#'
#' @param XColumn The name of the column containing x-coordinates of sites.
#' X-coordinates can be provided in either the biological or environmental data
#' tables, but MUST be in at least one of them. If an x-coordinate column is
#' provided in both the biological and environmental data, the column name must
#' be identical. Site coordinates ideally should be in a projected coordinate
#' system (i.e., not longitude-latitude) to ensure proper calculation of
#' geographic distances. Note that if you are using rasters, they must be in the
#' same coordinate system as the site coordinates.
#'
#' @param YColumn The name of the column containing y-coordinates of sample sites.
#' Y-coordinates can be provided in either the biological or environmental data
#' tables, but MUST be in at least one of them. If a y-coordinate column is
#' provided in both the biological and environmental data, the column name must
#' be identical. Site coordinates ideally should be in a projected coordinate
#' system (i.e., not longitude-latitude) to ensure proper calculation of
#' geographic distances. Note that if you are using rasters, they must be in the
#' same coordinate system as the site coordinates.
#'
#' @param sppColumn Only used if bioFormat = 2 (x, y, species list). The name of
#' the column containing unique name / identifier for each species.
#'
#' @param abundColumn If abundance = TRUE, this parameter identifies the column
#' containing the measure of abundance at each site. Only used if bioFormat = 2
#' (i.e., x, y, species list), though in the case of abundance data, the table
#' would have four columns: x, y, species, abundance.
#'
#' @param sppFilter Default = 0. To account for limited sampling effort at some
#' sites, sppFilter removes all sites at which the number of recorded species
#' (i.e., observed species richness) is less than the specified value. For
#' example, if sppFilter = 5, all sites with fewer than 5 recorded species
#' will be removed.
#'
#' @param predData The environmental predictor data. Accepts either a
#' site-by-predictor table or a terra object SpatRaster.
#'
#' @param distPreds An optional list of distance matrices to be used as predictors
#' in combination with predData. For example, a site-by-site dissimilarity matrix
#' for one biological group (e.g., trees) can be used as a predictor for another
#' group (e.g., ferns). Each distance matrix must have as the first column the
#' names of the sites (therefore the matrix will not be square). The name of the
#' column containing the site names should have the same name as that provided
#' for the siteColumn argument. Site IDs are required here to ensure correct ordering
#' of sites in the construction of the site-pair table. Note that the formatsitepair
#' function will not accept distance matrices in the as the only predictors
#' (i.e., at least one additional, non-distPreds predictor variable is required). If you wish to fit GDM
#' using only distance matrices provided using distPreds, provide one fake predictor (e.g., with all sites
#' having the same value), plus site and coordinate columns. The s1 and
#' s2 columns for this fake variable can then be removed by hand before fitting the GDM.
#'
#' @param weightType Default = "equal". Defines the weighting for sites. Can be
#' either: (1) "equal" (weights for all sites set = 1), (2) "richness" (each
#' site weighted according to number of species recorded), or (3) "custom"
#' (user defined). If weightType="custom", the user must provide a vector of
#' site weights equal to the number of rows in the full site-pair table (i.e.,
#' before species filtering (sppFilter argument) or sub-sampling is taken into
#' account (sampleSites argument)).
#'
#' @param custWeights A two column matrix or data frame of user-defined site
#' weights. The first column should be the site name and should be named the same
#' as that provided for the siteColumn argument. The second column should be numeric
#' weight values and should be named "weights". The weight values represent the
#' importance of each site in model fitting, and the values in the output
#' site-pair table is an average of the two sites in each site-pair. Required
#' when weightType = "custom". Ignored otherwise.
#'
#' @param sampleSites Default = 1. A number between 0-1 indicating the fraction
#' of sites to be used to construct the site-pair table. This argument can be
#' used to reduce the number of sites to overcome possible memory limitations
#' when fitting models with very large numbers of sites.
#'
#' @param verbose Default = FALSE. If TRUE, summary of information regarding
#' dimensions of the site-pair table will be printed that can be useful for diagnostics.
#'
#' @details
#' bioData and bioFormat:
#' The function accepts biological data in the following formats:
#'
#' bioData = site-by-species matrix; bioFormat = 1: assumes that the response
#' data are provided with a site ID column (specified by siteCol) and, optionally,
#' two columns for the x & y coordinates of the sites. All remaining columns
#' contain the biological data, with a column for each biological entity (most
#' commonly species). In the case that a raster stack (a terra object SpatRaster) is provided for the
#' environmental data (predData), x-y coordinates MUST be provided in bioData
#' to allow extraction of the environmental data at site locations. The x-y
#' coordinates will be intersected with the raster stack and, if the number of
#' unique cells intersected by the points is less than the number of unique site
#' IDs (i.e. multiple sites fall within a single cell), the function will use
#' the raster cell as the site ID and aggregate sites accordingly. Therefore,
#' model fitting will be sensitive to raster cell size. If the environmental
#' data are in tabular format, they should have the same number of sites
#' (i.e., same number of rows) as bioData. The x-y coordinate and site ID
#' columns must have the same names in bioData and predData.
#'
#' bioData = x, y, species list (optionally a fourth column with abundance can
#' be provided); bioFormat = 2: assumes a table of 3 or 4 columns, the first two
#' being the x & y coordinates of species records, the third (sppCol) being the
#' name / identifier of the species observed at that location, and optionally a
#' fourth column indicating a measure of abundance.  If an abundance column is
#' not provided, presence-only data are assumed. In the case that a raster stack
#' (a terra object SpatRaster) is provided for the environmental data (predData),
#' the x-y coordinates will be intersected with the raster stack and, if the
#' number of unique cells intersected by the points is less than the number of
#' unique site IDs (i.e. multiple sites fall within a single cell), the function
#' will use the raster cell as the site ID and aggregate sites accordingly.
#' Therefore, model fitting will be sensitive to raster cell size.
#'
#' bioData = site-site distance (dissimilarity) matrix; bioFormat = 3. This option
#' allows the use of an existing site-site distance (dissimilarity) matrix, such as
#' genetic distance matrix calculated outside of the gdm package. Only the lower
#' triangle of the matrix is required to create the site-pair table, but the
#' function will automatically removes the upper triangle if present. The code
#' checks and aligns the order of sites in the distance matrix and the predictor
#' data to ensure they match. To do so, (1) a site column is required in both
#' the distance matrix and the predictor data and (2) site IDs are required to
#' be a number. This is the only bioFormat in which the environmental data MAY
#' NOT be a raster stack.
#'
#' bioData = site-pair table; bioFormat = 4: with an already created site-pair
#' table, this option allows the user to add one or more distance matrices (see
#' distPreds above) to the existing site-pair table and/or sub-sample the
#' site-pair table (see sample above). If the site-pair table was not created
#' using the formatsitepair function, the user will need to ensure the order of
#' the sites matches that in other tables being provided to the function.
#'
#' NOTES: (1) The function assumes that the x-y coordinates and the raster stack
#' (if used) are in the same coordinate system. No checking is performed to
#' confirm this is the case. (2) The function assumes that the association between
#'  the provided site and x-y coordinate columns are singular and unique.
#'  Therefore, the function will fail should a given site has more than one sets of
#'  coordinates associated with it, as well as multiple sites being given the
#'  exact same coordinates.
#'
#' @return A formatted site-pair table containing the response (biological
#' distance or dissimilarity), predictors, and weights as required for fitting
#' Generalized Dissimilarity Models.
#'
#' @examples
#' ## tabular data
#' # start with the southwest data table
#'  head(southwest)
#'  sppData <- southwest[c(1,2,13,14)]
#'  envTab <- southwest[c(2:ncol(southwest))]
#'
#' #########table type 1
#' ## site-species table without coordinates
#' testData1a <- reshape2::dcast(sppData, site~species)
#' ##site-species table with coordinates
#' coords <- unique(sppData[, 2:ncol(sppData)])
#' testData1b <- merge(testData1a, coords, by="site")
#' ## site-species, table-table
#' exFormat1a <- formatsitepair(testData1a, 1, siteColumn="site", XColumn="Long",
#' YColumn="Lat", predData=envTab)
#'
#' #' # next, let's try environmental raster data
#' ## not run
#' # rastFile <- system.file("./extdata/swBioclims.grd", package="gdm")
#' # envRast <- terra::rast(rastFile)
#'
#' ## site-species, table-raster
#' ## not run
#' # exFormat1b <- formatsitepair(testData1b, 1, siteColumn="site", XColumn="Long",
#' # YColumn="Lat", predData=envRast)
#'
#' #########table type 2
#' ## site xy spp list, table-table
#' exFormat2a <- formatsitepair(sppData, 2, XColumn="Long", YColumn="Lat",
#' sppColumn="species", siteColumn="site", predData=envTab)
#'
#' ## site xy spp list, table-raster
#' ## not run
#' # exFormat2b <- formatsitepair(sppData, 2, XColumn="Long", YColumn="Lat",
#' # sppColumn="species", siteColumn="site", predData=envRast)
#'
#' #########table type 3
#' ## It is possible to format a site-pair table by starting
#' # with a pre-calculated matrix of biological distances
#' dim(gdmDissim) # pairwise distance matrix + 1 column for site IDs
#' gdmDissim[1:5, 1:5]
#' # now we can format the table:
#' exFormat3 <- formatsitepair(gdmDissim, 3, XColumn="Long", YColumn="Lat",
#'                             predData=envTab, siteColumn="site")
#'
#' #########table type 4
#' ## adds a predictor matrix to an existing site-pair table, in this case,
#' ## predData needs to be provided, but is not actually used
#' exFormat4 <- formatsitepair(exFormat2a, 4, predData=envTab, siteColumn="site",
#'                             distPreds=list(as.matrix(gdmDissim)))
#'
#' @keywords gdm
#'
#' @importFrom methods is
#' @importFrom stats as.dist
#' @importFrom vegan vegdist
#' @importFrom reshape2 dcast
#'
#' @export
formatsitepair <- function(bioData, bioFormat, dist="bray", abundance=FALSE,
                           siteColumn=NULL, XColumn, YColumn, sppColumn=NULL,
                           abundColumn=NULL, sppFilter=0, predData, distPreds=NULL,
                           weightType="equal", custWeights=NULL, sampleSites=1,
                           verbose=FALSE){
  ###########################
  ##lines used to quickly test function
  #bioData <- southwest[, c(1,2,13,14)]
  #bioFormat <- 2
  #dist <- "bray"
  #abundance <- F
  #siteColumn <- "site"
  #XColumn <- "Long"
  #YColumn <- "Lat"
  #sppColumn <- "species"
  #sppFilter <- 0
  #abundColumn <- NULL
  #predData <- southwest[, c(2:ncol(southwest))]
  #distPreds <- NULL
  #weightType <- "equal"
  #custWeights <- NULL
  #sampleSites <- 1
  #################
  #bioData <- bio
  #bioFormat <- 2
  #dist <- "bray"
  #abundance <- T
  #siteColumn <- "site"
  #XColumn <- "X"
  #YColumn <- "Y"
  #sppColumn <- "locus"
  #sppFilter <- 0
  #abundColumn <- "abund"
  #predData <- preds
  #distPreds <- NULL
  #weightType <- "equal"
  #custWeights <- NULL
  #sampleSites <- 1
  ###########################
  ##input error checking
  ##makes sure bioData has the required format
  if(!(is(bioData, "data.frame") | is(bioData, "matrix") | is(bioData, "gdmData"))){
    "bioData object needs to be either of class data.frame or matrix in one of
    the acceptable formats listed in the help document."
  }
  ##transforms bioData into a data frame in order to proceed without error, just to make sure,
  ##basically remove the affects of factors
  if(is(bioData, "data.frame") | is(bioData, "matrix")){
    bioData <- as.data.frame(bioData, stringsAsFactors=F)
  }

  ##makes sure predData is in an acceptable format
  if(!(is(predData, "data.frame") | is(predData, "matrix") | .is_raster(predData))){
    "predData object needs to either of class data.frame, matrix, or raster."
  }
  if(is(predData, "data.frame") | is(predData, "matrix")){
    predData <- as.data.frame(predData, stringsAsFactors=F)
  }

  # make sure predData is a terra object or convert it to terra
  if (.is_raster(predData)) {
    # check terra package is available
    .check_pkgs("terra")
    predData <- .check_rast(predData, "predData")
  }

  ##if bioFormat is not an acceptable number, exit function
  if(bioFormat %in% c(1:4)){} else{
    stop("Acceptable values for the bioFormat argument are: 1, 2, 3, or 4")
  }
  ##checks that geo has either TRUE or FALSE
  if(!(abundance==TRUE | abundance==FALSE)){
    stop("abundance argument must be either TRUE or FALSE")
  }

  ##if sampleSites is not a number, then exit function
  #if(is.null(sampleSites)==TRUE){
  #  stop("sampleSites argument must be a number between 0-1")
  #}
  if (is.numeric(sampleSites) == FALSE | sampleSites <= 0 |
      sampleSites > 1) {
    stop("sampleSites argument must be a number 0 < x <= 1")
  }

  ##makes sure that sppFilter is a number, if not exit function
  if((is.numeric(sppFilter)==FALSE & is.null(sppFilter)==FALSE) | sppFilter<0){
    stop("sppFilter argument must be a positive integer")
  }
  ##makes sure a proper weightType is used
  if(weightType %in% c("equal", "richness", "custom")){} else{
    stop("Acceptable values for the weightType argument are: equal, richness, or custom")
  }
  ##if weightType == custom, makes sure a custWeights is attached
  if(weightType=="custom" & is.null(custWeights)==T){
    stop("weightType argument = 'custom', but no custWeights vector provided.")
  }
  ##with bioFormat 2, makes sure spp data has been included
  if(bioFormat==2 & is.null(sppColumn)==TRUE){
    stop("Need to define sppColumn argument when bioFormat==2")
  }
  ##makes sure that the sppColumn name can be found in the bioData with bioFormat 2
  #if(bioFormat==2 & (sppColumn %in% names(bioData))){} else{
  #  stop("Cannot find sppColumn in bioData - check name?")
  #}
  ##makes sure that a site column is provided when using table type 2 and raster environmental data
  if(bioFormat==2 & is.null(siteColumn)==TRUE){
    if(!.is_raster(predData)){
      stop("A siteColumn needs to be provided in either the bioData or predData objects.")
    }
  }
  ##when a site column is provided
  if(is.null(siteColumn)==FALSE){
    ##makes sure the site column name is of type character
    if(!is(siteColumn, "character")){
      stop("siteColumn argument needs to be of class = 'character'.")
      ##checks to see if siteColumn exists in the bioData for bioFormats 1 and 2
    }else if(!(siteColumn %in% colnames(bioData)) & (bioFormat==1 | bioFormat==2)){
      stop("Cannot find a match for the name of the siteColumn in the columns
           of the bioData object.")
    }
    ##if the siteColumn is provided with input type 3, remove it
    #if(bioFormat==3 & siteColumn %in% colnames(bioData)){
    #  wSite <- which(colnames(bioData)==siteColumn)
    #  bioData <- bioData[-wSite]
    #}
  }
  ##checks to make sure that the coordinate columns are characters and can be
  # found in either the biological or environmental data
  if(bioFormat!=4){
    if(!is(XColumn, "character")){
      stop("XColumn argument needs to be of class 'character'.")
    }else if(!is(YColumn, "character")){
      stop("YColumn argument needs to be of class 'character'.")
    }else if(!(XColumn %in% colnames(bioData) | XColumn %in% colnames(predData))){
      stop("XColumn not found in either the bioData or predData arguments")
    }else if(!(YColumn %in% colnames(bioData) | YColumn %in% colnames(predData))){
      stop("YColumn not found in either the bioData or predData arguments")
    }
  }
  ##checks table type 3 specific requirements
  if(bioFormat==3){
    if(weightType=="richness"){
      stop("Cannot weight by site richness when supplying the biological data
           as a distance matrix.")
    }else if(nrow(bioData)!=(ncol(bioData)-1)){
      stop("Check dimensions of the bioData object. Does the biological dissimilarity matrix have a column for site IDs as required?")
    }
  }

  ##warns if distPreds are not matrices
  for(mat in distPreds){
    if(!is(mat, "matrix") & !is(mat, "data.frame")){
      warning("One or more of the provided predictor distance matrices are not
              of class 'matrix'.")
    }
  }
  ##if a custom weight vector is provided, makes sure it is a vector
  if(is.null(custWeights)==FALSE & (!is(custWeights, "data.frame") & !is(custWeights, "matrix"))){
    stop("The argument custWeights needs to be of class 'data.frame' or 'matrix'.")
  }

  ##sets up variables to be used later
  toRemove <- NULL       ##removed from sppfilter
  removeRand <- NULL     ##remove from random sample out sites
  distData <- NULL
  ##checks input data format
  ##species data as site-species matrix
  if(bioFormat==1 | bioFormat==2){
    ##first, if bioFormat=2, then transforms it into a site-spp matrix (bioFormat 1)
    if(bioFormat==2){
      ##makes sure that the sppColumn name can be found in the bioData with bioFormat 2
      if((sppColumn %in% colnames(bioData))){} else{
        stop("Cannot find sppColumn in bioData - check name?")
      }

      ##insert a site column if one was not given
      if(is.null(siteColumn)){
        colnames(bioData)[which(colnames(bioData)==XColumn)] <- "bioData_FSP_xColumn"
        colnames(bioData)[which(colnames(bioData)==YColumn)] <- "bioData_FSP_yColumn"
        bioData <- transform(bioData,
                             griddedSiteID=as.numeric(interaction(bioData$bioData_FSP_xColumn,
                                                                         bioData$bioData_FSP_yColumn,
                                                                         drop=TRUE)))
        siteColumn <- "griddedSiteID"
        colnames(bioData)[which(colnames(bioData)=="bioData_FSP_xColumn")] <- XColumn
        colnames(bioData)[which(colnames(bioData)=="bioData_FSP_yColumn")] <- YColumn
      }

      ##insert presence if abundance was not given
      if(is.null(abundColumn)){
        warning("No abundance column was specified, so the biological data are assumed to be presences.")
        bioData[, "bioData_FSP_abundColumn"] <- 1
        abundColumn <- "bioData_FSP_abundColumn"
      }

      ##rename the siteColumn and sppColumn in order to cast the data into a siteXspp matrix
      preCastBio <- bioData
      colnames(preCastBio)[which(colnames(preCastBio)==siteColumn)] <- "griddedSiteID"
      colnames(preCastBio)[which(colnames(preCastBio)==sppColumn)] <- "bioData_FSP_sppColumn"
      castData <- dcast(preCastBio, fill=0, griddedSiteID~bioData_FSP_sppColumn, value.var=abundColumn)
      ##adds coordinates to the cast data
      uniqueCoords <- unique(preCastBio[which(colnames(preCastBio) %in% c("griddedSiteID", XColumn, YColumn))])
      bioData <- merge(castData, uniqueCoords, by="griddedSiteID")
      colnames(bioData)[which(colnames(bioData)=="griddedSiteID")] <- siteColumn
    }

    ##checks to see if the coordinates can be found in bioData, if not, checks to see if
    ##they can be found in envData
    if((XColumn %in% colnames(bioData))==FALSE | (YColumn %in% colnames(bioData)==FALSE)){
      xCol <- which(colnames(predData)==XColumn)
      yCol <- which(colnames(predData)==YColumn)
      locs <- predData[c(xCol,yCol)]
    }else{
      xCol <- which(colnames(bioData)==XColumn)
      yCol <- which(colnames(bioData)==YColumn)
      locs <- bioData[c(xCol,yCol)]
    }

    # checks unique sites against rasters
    if(.is_raster(predData)){

      # when using rasters, uses the cell as the site
      warning("When using rasters for environmental covariates (predictors), each site is assigned to the
              raster cell in which the site is located. If more than one site occurs within the same raster cell,
              the biological data of those sites are aggregated (more likely as raster resolution decreases).")
      # gets the cell location of the given coordinates
      cellID <- data.frame(cellName = terra::cellFromXY(predData, locs))

      # if none of the points intersected with the predictor raster
      if (all(is.na(cellID$cellName))) {
        stop("None of the data points provided intersect with the rasters. Double check spatial data.")
      }

      cellLocs <- as.data.frame(terra::xyFromCell(predData, cellID$cellName))
      ##temporarily keeps old site in to identify what to remove from other objects
      rastBioData <- cbind(cellID, cellLocs, bioData[-c(which(colnames(bioData) %in% c(XColumn, YColumn)))])

      ##if custom weights selected, gives weight table new cell site names
      if(weightType=="custom" & !is.null(custWeights)){
        nameTab <- unique(rastBioData[c("cellName", siteColumn)])
        tempWeightTab <- merge(x=nameTab, y=custWeights, by=siteColumn)
        siteNum <- which(colnames(tempWeightTab)=="cellName")
        custWeights <- tempWeightTab[-siteNum]
        colnames(custWeights)[1] <- siteColumn
      }
      ##removes original site column from the rastBioData table
      siteNum <- which(colnames(rastBioData)==siteColumn)
      rastBioData <- rastBioData[-siteNum]

      ##aggregates species data by cell
      cellNum <- which(colnames(rastBioData)=="cellName")
      bioData <- stats::aggregate(rastBioData, rastBioData[cellNum], FUN=mean)
      bioData <- bioData[-cellNum]

      # extracts raster data into environmental predictor data table
      rastEx <- terra::extract(predData, bioData$cellName)

      ##renames bioData columns which have been updated from rasters
      colnames(bioData)[which(colnames(bioData)=="cellName")] <- siteColumn
      colnames(bioData)[which(colnames(bioData)=="x")] <- XColumn
      colnames(bioData)[which(colnames(bioData)=="y")] <- YColumn

      ##updates locs object based on raster coordinates
      xCol <- which(colnames(bioData)==XColumn)
      yCol <- which(colnames(bioData)==YColumn)
      locs <- bioData[c(xCol,yCol)]

      ##recreates the predData
      siteCol <- which(colnames(bioData)==siteColumn)
      predData <- cbind(bioData[siteCol], locs, rastEx)
    }

    ##filters out sites with low species count
    ##first isolates the species data
    siteCol <- which(colnames(bioData)==siteColumn)
    xCol <- which(colnames(bioData)==XColumn)
    yCol <- which(colnames(bioData)==YColumn)
    sppDat <- bioData[, -c(siteCol, xCol, yCol)]
    ##totals the number of species per site
    sppDat[sppDat>=1] <- 1
    sppDat[sppDat==0] <- 0
    sppDat[is.na(sppDat)] <- 0
    sppTotals <- cbind.data.frame(bioData[,siteCol], apply(sppDat, 1, function(m){sum(as.numeric(m))}))
    colnames(sppTotals) <- c(siteColumn, "totals")
    ##filters out sites with less species than filter
    filterBioDat <- subset(sppTotals, sppTotals[colnames(sppTotals)[2]] >= sppFilter)
    toRemove <- bioData[,siteCol][which(!(bioData[,siteCol] %in% filterBioDat[,siteColumn]))]
    ##reassembles bioData after filtering
    spSiteCol <- filterBioDat[1]
    bioData <- unique(merge(spSiteCol, bioData, by=siteColumn))

    ##subsample sites using random sampling
    if (sampleSites < 1) {
      fullSites <- bioData[, siteCol]
      randRows <- sort(sample(1:nrow(bioData),
                              round(nrow(bioData) * sampleSites, 0)))
      ##actual selection of the random rows to keep
      bioData <- bioData[c(randRows),]
      #removeRand <- fullLength[-(randRows)]
      ##records the sites that have been removed, for distPreds later in function
      removeRand <- fullSites[which(! (fullSites %in% bioData[,siteCol]))]
    }

    ##identifies and removes filtered out sites and sampled sites from predData
    ##renames siteColumn in order to access objects correctly
    colnames(bioData)[colnames(bioData)==siteColumn] <- "theSiteColumn"
    colnames(predData)[colnames(predData)==siteColumn] <- "theSiteColumn"
    predData <- unique(predData)
    #predData3 <- predData[which(predData$theSiteColumn %in%
    #                             as.character(as.numeric(bioData$theSiteColumn))),]
    predData <- predData[which(predData$theSiteColumn %in%
                                 as.character(bioData$theSiteColumn)),]

    ##remove custom weights from any sites removed by species filtering and sampling
    if(weightType=="custom" & !is.null(custWeights)){
      colnames(custWeights)[colnames(custWeights)==siteColumn] <- "theSiteColumn"
      custWeights <- custWeights[which(predData$theSiteColumn %in% custWeights[,"theSiteColumn"]),]
      hwap <- custWeights[,"theSiteColumn"]#[,1]
      hwap <- order(hwap)
      custWeights <- custWeights[hwap,]
      colnames(custWeights)[colnames(custWeights)=="theSiteColumn"] <- siteColumn
    }

    ##rename site columns for results
    colnames(bioData)[colnames(bioData)=="theSiteColumn"] <- siteColumn
    colnames(predData)[colnames(predData)=="theSiteColumn"] <- siteColumn

    ##as a final check, makes sure bioData and predData sites are in same order
    predSite <- which(names(predData) == siteColumn)
    bioSite <- which(names(bioData)==siteColumn)
    hwap <- predData[,predSite]
    hwap <- order(hwap)
    predData <- predData[hwap,]
    rosetta <- bioData[,bioSite]
    rosetta <- order(rosetta)
    bioData <- bioData[rosetta,]

    ##sets up species data for calculating dissimilarity
    bx <- which(names(bioData)==XColumn)
    by <- which(names(bioData)==YColumn)
    sppData <- bioData[-c(bioSite, bx, by)]

    ##creates distance matrix
    if(abundance==F){
      sppData[sppData>=1] <- 1
      sppData[sppData==0] <- 0
      sppData[is.na(sppData)] <- 0
      distData <- vegdist(sppData, dist, binary=T)
    }else{
      sppData[is.na(sppData)] <- 0
      distData <- vegdist(sppData, dist, binary=F)
    }

    ########################################################################
    ##species data as site-by-site distance matrix
  }else if(bioFormat==3){
    holdSite <- bioData[, which(siteColumn %in% colnames(bioData))]
    if(!is.numeric(holdSite)){
      stop("Site IDs must be numeric when providing a site-by-site distance matrix (bioFormat=3).")
    }
    bioData <- bioData[, -which(siteColumn %in% colnames(bioData))]
    orderedData <- as.matrix(as.dist(bioData[order(holdSite),
                                             order(holdSite)]))
    distData <- lower.tri(as.matrix(orderedData), diag = FALSE)
    distData <- as.vector(orderedData[distData])
    predData <- unique(predData)
    hwap <- predData[siteColumn][, 1]
    hwap <- order(hwap)
    predData <- predData[order(holdSite), ]
    ########################################################################
    ##site pair table, already prepped
  }else if(bioFormat==4){
    ##site-pair distance value
    outTable <- bioData
    ########################################################################

  }else{
    ##return error, bioFormat argument out of bounds
    stop(paste("bioFormat argument of '", as.character(bioFormat), "' is not an accepted input value", sep=""))
  }
  ########################################################################

  ##With the dissim distance calculated, creates and fills the table in gdm format
  if(bioFormat!=4){
    ##creates base site-pair table
    outTable <- as.data.frame(createsitepair(dist=distData,
                                             spdata=bioData,
                                             envInfo=predData,
                                             dXCol=XColumn,
                                             dYCol=YColumn,
                                             siteCol=siteColumn,
                                             weightsType=weightType,
                                             custWeights=custWeights))
  }else{
    outTable <- bioData
  }

  ##first checks that the size of the dissimilarity matrices, if any were provided
  ##then pastes any dissimilarity matrices onto the created site-pair table
  if(length(distPreds)>0){

    baseMat <- distPreds[[1]]
    ##checks to size of each dissimilarity matrix, to make sure they are all the same
    lapply(distPreds, function(mat, mat1){
      if((dim(mat1)[1]!=dim(mat)[1]) & (dim(mat1)[2]!=dim(mat)[2])){
        stop("The dimensions of your predictor matrices are not the same.")
      }
    }, mat1=baseMat)
    #print(dim(baseMat))

    ##hold site columns
    holdSiteCols <- lapply(distPreds, function(dP){dP[,which(siteColumn %in% colnames(dP))]})
    #remove site column from matrices
    distPreds <- lapply(distPreds, function(dP){dP[,-which(siteColumn %in% colnames(dP))]})
    #print(dim(distPreds[[1]]))
    ##orders the distance matrices of distPreds
    distPreds <- mapply(function(dP, hSC){as.matrix(as.dist(dP[order(hSC),order(hSC)]))},
                        dP=distPreds, hSC=holdSiteCols, SIMPLIFY=FALSE)
    #print(dim(distPreds[[1]]))
    ##orders the site columns to match the distance matrices
    orderSiteCols <- lapply(holdSiteCols, function(hSC){hSC[order(hSC)]})

    ##compares sites with sites removed, for one reason or another
    rmSites <- c(toRemove, removeRand)
    ##removes the sites removed above when creating the site-pair table
    if(length(rmSites)>0){
      rmIndex <- lapply(orderSiteCols, function(hSC, tR){which((hSC %in% tR))}, tR=rmSites)
      distPreds <- mapply(function(mat, tR){mat[-c(tR), -c(tR)]}, mat=distPreds, tR=rmIndex, SIMPLIFY=FALSE)
    }
    ##set new baseMat
    baseMat <- distPreds[[1]]
    #print(dim(baseMat))

    ##checks the size of the dissimilarity matrices against the size of distData
    baseMatDat <- lower.tri(as.matrix(baseMat),diag=FALSE)
    baseMatDat <- as.vector(baseMat[baseMatDat])
    #print(nrow(outTable))
    #print(length(baseMatDat))
    if(nrow(outTable)!=length(baseMatDat)){
      stop("The dimensions of the distance predictor matrices do not match the biological data.")
    }

    ##adds any associated distance predictors to the sitepair table
    for(num in 1:length(distPreds)){
      #num <- 2
      ##isolate matrix
      matrixDat <- lower.tri(as.matrix(distPreds[[num]], diag=FALSE))
      sweetFreakinPredMatrix <- as.vector(distPreds[[num]][matrixDat])
      ##add matrix to table
      if(ncol(outTable)>6){
        ##break table up to insert matrix data columns
        baseSitePair <- outTable[,1:6]
        otherSitePair <- outTable[,7:ncol(outTable)]
        otherNames <- colnames(otherSitePair)
        s1SitePair <- as.data.frame(otherSitePair[,1:(ncol(otherSitePair)/2)])
        colnames(s1SitePair) <- otherNames[1:(ncol(otherSitePair)/2)]
        s2SitePair <- as.data.frame(otherSitePair[,(ncol(otherSitePair)/2+1):ncol(otherSitePair)])
        colnames(s2SitePair) <- otherNames[(ncol(otherSitePair)/2+1):ncol(otherSitePair)]
        ##formats data from dissimilarity matrices
        s1Zeros <- as.data.frame(rep(0,length(sweetFreakinPredMatrix)))
        colnames(s1Zeros) <- paste("s1.matrix_", num, sep="")
        s2Mat <- as.data.frame(sweetFreakinPredMatrix)
        colnames(s2Mat) <- paste("s2.matrix_", num, sep="")
        ##restructures the output table
        outTable <- cbind(baseSitePair, s1SitePair, s1Zeros, s2SitePair, s2Mat)
      }else{
        ##formats data from dissimilarity matrices
        s1Zeros <- as.data.frame(rep(0,length(sweetFreakinPredMatrix)))
        colnames(s1Zeros) <- paste("s1.matrix_", num, sep="")
        s2Mat <- as.data.frame(sweetFreakinPredMatrix)
        colnames(s2Mat) <- paste("s2.matrix_", num, sep="")
        ##restructures the output table
        outTable <- cbind(outTable, s1Zeros, s2Mat)
      }
    }
  }

  if(verbose){
    if(weightType[1]=="equal"){
      print("Site weighting type: Equal")
    }else if(weightType[1]=="custom"){
      print("Site weighting type: Custom")
    }else{
      print("Site weighting type: Richness")
    }
    print(paste0("Site-pair table created with ", nrow(outTable), " rows ",
                 "(", nrow(unique(outTable[,3:4]))+1, " unique sites)", " and ",
                 ncol(outTable) , " columns (", (ncol(outTable)-6)/2,
                 " environmental variables)."))
  }

  ##return output table
  class(outTable) <- c("gdmData", "data.frame")
  return(outTable)
}

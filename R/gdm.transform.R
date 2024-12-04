#' @title Transform Environmental Data Using a Fitted Generalized Dissimilarity Model
#'
#' @description This function transforms geographic and environmental predictors using (1) the
#' fitted functions from a model object returned from \code{\link[gdm]{gdm}} and (2) a
#' data frame or raster object containing predictor data for a set of sites.
#'
#' @usage gdm.transform(model, data, filename = "", ...)
#'
#' @param model A gdm model object resulting from a call to \code{\link[gdm]{gdm}}.
#'
#' @param data Either (i) a data frame containing values for each predictor variable in the model, formatted as follows: X, Y, var1, var2, var3, ..., varN or
#'  (ii) a terra object SpatRaster with one layer per predictor variable used in the model,
#'  excluding X and Y (rasters for x- and y-coordinates are built automatically from the input
#'  rasters if the model was fit with geo=TRUE). The order of the columns (data frame) or raster layers (SpatRaster) MUST be the same as the order of the predictors in
#'  the site-pair table used in model fitting. There is currently no checking to ensure that the order
#'  of the variables to be transformed are the same as those in the site-pair table used in model fitting.
#'  If geographic distance was not used as a predictor in model fitting, the x- and y-columns
#'  need to be removed from the data to be transformed.
#'  Output is provided in the same format as the input data.
#'
#' @param filename character. Output filename for rasters. When provided the raster layers are
#' written to file directly.
#'
#' @param ... additional arguments to pass to terra \code{\link[terra]{predict}} function.
#'
#' @return
#' gdm.transform returns either a data frame with the same number of rows as the input data frame or a SpatRaster,
#' depending on the format of the input data. If the model uses geographic distance as a predictor the output object
#' will contain columns or layers for the transformed X and Y values for each site.
#' The transformed environmental data will be in the remaining columns or layers.
#'
#' @references Ferrier S, Manion G, Elith J, Richardson, K (2007) Using generalized dissimilarity modelling to analyse and predict patterns of beta diversity in regional biodiversity assessment. \emph{Diversity & Distributions} 13, 252-264.
#'
#' Fitzpatrick MC, Keller SR (2015) Ecological genomics meets community-level modeling of biodiversity: Mapping the genomic landscape of current and future environmental adaptation. \emph{Ecology Letters} 18: 1-16
#'
#' @examples
#' # start with the southwest data set
#' # grab the columns with xy, site ID, and species data
#' sppTab <- southwest[, c("species", "site", "Lat", "Long")]
#'
#' ##fit gdm using rasters
#' rastFile <- system.file("./extdata/swBioclims.grd", package="gdm")
#' envRast <- terra::rast(rastFile)
#'
#' sitePairRast <- formatsitepair(sppTab, 2, XColumn="Long", YColumn="Lat",
#'                                sppColumn="species", siteColumn="site",
#'                                predData=envRast)
#' ##remove NA values
#' sitePairRast <- na.omit(sitePairRast)
#'
#' ##fit raster GDM
#' gdmRastMod <- gdm(sitePairRast, geo=TRUE)
#'
#' ##raster input, raster output
#' transRasts <- gdm.transform(gdmRastMod, envRast)
#'
#' # map biological patterns; increase maxcell if using large rasters
#' pcaSamp <- terra::prcomp(transRasts, maxcell = 1e4)
#'
#' # note the use of the 'index' argument
#' pcaRast <- terra::predict(transRasts, pcaSamp, index=1:3)
#'
#' # stretch the PCA rasters to make full use of the colour spectrum
#' pcaRast <- terra::stretch(pcaRast)
#'
#' terra::plotRGB(pcaRast, r=1, g=2, b=3)
#'
#' @keywords gdm
#'
#' @export
gdm.transform <- function(model, data, filename = "", ...){
  options(warn.FPU = FALSE)

  # error checking of inputs
  # checks to make sure a gdm model is given
  if(!is(model, "gdm")){
    stop("model object must be of class 'gdm'.")
  }

  if(!(.is_raster(data) | is(data, "data.frame"))){
    stop("Data to be transformed must either be a raster object or data frame.")
  }

  if(!is(data, "data.frame")) {
    # check whether terra package is available
    .check_pkgs("terra")
    # if data is raster (and not dataframe), convert to terra if it's not
    data <- .check_rast(data, "data")
  }

  # checks rather geo was T or F in the model object
  geo <- model$geo

  # produce outputs
  if (.is_raster(data)) {
    # transform the env data using gdm model and terra package
    # if geo, use interpolate to get xy; otherwise use predict without xy
    if (geo) {
      # get the extent to pass it to the function for correcting the transformed XY coords
      rast_ext <- terra::ext(data[[1]])[1:4]

      output <- terra::interpolate(
        object = data,
        model = model,
        fun = gdm_trans,
        xy_range = rast_ext,
        xyNames = c("xCoord", "yCoord"),
        na.rm = FALSE,
        filename = filename,
        ...
      )
    } else {
      output <- terra::predict(
        object = data,
        model = model,
        fun = gdm_trans,
        na.rm = FALSE,
        filename = filename,
        ...
      )
    }

    # get the predictors with non-zero sum of coefficients
    splineindex <- 1
    predInd <- c()
    for(i in 1:length(model$predictors)){
      # only if the sum of the coefficients associated with this predictor is > 0.....
      numsplines <- model$splines[i]
      coeff_sum <- sum(model$coefficients[splineindex:(splineindex + numsplines - 1)])

      if (coeff_sum > 0) {
        predInd <- c(predInd, i)
      }

      splineindex <- splineindex + numsplines
    }

    if(geo){
      predInd <- c(1,2,predInd[-1]+1)
    }

    # subset non-zero raster layers
    output <- terra::subset(output, predInd)

    return(output)

  } else {

    # transform the env data using gdm model using a data.frame
    output <- gdm_trans(
      mod = model,
      dat = data
    )

    colnames(output) <- colnames(data)

    return(output)
  }

}


# create a general predict function to benefit from terra::predict function
gdm_trans <- function(mod, dat, xy_range = NULL, ...) {
  # check for NAs
  has_na <- anyNA(dat)
  nm <- colnames(dat)

  if (has_na) {
    # output matrix
    out <- matrix(NA, nrow(dat), ncol(dat))
    colnames(out) <- nm
    # get the complete case and filter the data
    cc <- which(stats::complete.cases(dat))
    dat <- dat[cc, ]

    if (nrow(dat) < 1) {
      return(out)
    }
  }

  # Use the XY range to correctly calculate transformed XY coords in chunks
  # only for model with geo = true
  if (!is.null(xy_range)) {
    # get the first row twice in case there's only one row
    endrows <- dat[c(1, 1), ]
    endrows$xCoord <- xy_range[1:2]
    endrows$yCoord <- xy_range[3:4]
    # add the two row at the end of the data.frame
    dat <- rbind(dat, endrows)
  }

  nr <- nrow(dat)
  nc <- ncol(dat)

  z <- .C("GDM_TransformFromTable", as.integer(nr), as.integer(nc),
          as.integer(mod$geo), as.integer(length(mod$predictors)),
          as.integer(mod$splines), as.double(mod$knots), as.double(mod$coefficients),
          as.matrix(dat), trandata = as.double(matrix(0.0, nr, nc)), PACKAGE = "gdm")

  # prepare the output matrix
  transformed <- matrix(z$trandata, nrow = nr, byrow = FALSE)
  colnames(transformed) <- nm

  # remove the xy range if it's a model with geo
  if (!is.null(xy_range)) {
    transformed <- transformed[1:(nr - 2), ]
  }

  if (has_na) {
    out[cc, ] <- transformed
    return(out)
  } else {
    return(transformed)
  }
}


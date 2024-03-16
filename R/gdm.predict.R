#' @title Predict Biological Dissimilarities Between Sites or Times Using a
#' Fitted Generalized Dissimilarity Model
#'
#' @description This function predicts biological distances between sites or times using a
#'  model object returned from \code{\link[gdm]{gdm}}. Predictions between site
#'  pairs require a data frame containing the values of predictors for pairs
#'  of locations, formatted as follows: distance, weights, s1.X, s1.Y, s2.X,
#'  s2.Y, s1.Pred1, s1.Pred2, ..., s1.PredN, s2.Pred1, s2.Pred2, ..., s2.PredN, ...,
#'  Predictions of biological change through time require two raster stacks or
#'  bricks for environmental conditions at two time periods, each with a
#'  layer for each environmental predictor in the fitted model.
#'
#' @usage \method{predict}{gdm}(object, data, time=FALSE, predRasts=NULL, ...)
#'
#' @param object A gdm model object resulting from a call to \code{\link[gdm]{gdm}}.
#'
#' @param data Either a data frame containing the values of predictors for pairs
#' of sites, in the same format and structure as used to fit the model using
#' \code{\link[gdm]{gdm}} or a raster stack if a prediction of biological change
#' through time is needed.
#'
#' For a data frame, the first two columns - distance and weights - are required
#' by the function but are not used in the prediction and can therefore be filled
#' with dummy data (e.g. all zeros). If geo is TRUE, then the s1.X, s1.Y and s2.X,
#' s2.Y columns will be used for calculating the geographical distance between
#' each site for inclusion of the geographic predictor term into the GDM model.
#' If geo is FALSE, then the s1.X, s1.Y, s2.X and s2.Y data columns are ignored.
#' However these columns are still REQUIRED and can be filled with dummy data
#' (e.g. all zeroes). The remaining columns are for N predictors for Site 1 and
#' followed by N predictors for Site 2. The order of the columns must match those
#' in the site-pair table used to fit the model.
#'
#' A raster stack should be provided only when time=T and should contain one
#' layer for each environmental predictor in the same order as the columns in
#' the site-pair table used to fit the model.
#'
#' @param time TRUE/FALSE: Is the model prediction for biological change through time?
#'
#' @param predRasts A raster stack characterizing environmental conditions for a
#' different time in the past or future, with the same extent, resolution, and
#' layer order as the data object. Required only if time=T.
#'
#' @param ... additional arguments that can be accepted by terra::predict.
#'
#' @return predict returns either a response vector with the same length as the
#'  number of rows in the input data frame or a raster depicting change through time across the study region.
#'
#' @seealso \code{\link[gdm]{gdm.transform}}
#'
#' @examples
#' ##set up site-pair table using the southwest data set
#' sppData <- southwest[, c(1,2,14,13)]
#' envTab <- southwest[, c(2:ncol(southwest))]
#'
#' # remove soils (no rasters for these)
#' envTab <- envTab[,-c(2:6)]
#' sitePairTab <- formatsitepair(sppData, 2, XColumn="Long", YColumn="Lat", sppColumn="species",
#'                              siteColumn="site", predData=envTab)
#'
##create GDM
#' gdmMod <- gdm(sitePairTab, geo=TRUE)
#'
#' ##predict GDM
#' predDiss <- predict(gdmMod, sitePairTab)
#'
#' ##time example
#' rastFile <- system.file("./extdata/swBioclims.grd", package="gdm")
#' envRast <- terra::rast(rastFile)
#'
#' ##make some fake climate change data
#' futRasts <- envRast
#' ##reduce winter precipitation by 25%
#' futRasts[[3]] <- futRasts[[3]]*0.75
#'
#' timePred <- predict(gdmMod, envRast, time=TRUE, predRasts=futRasts)
#' terra::plot(timePred)
#'
#' @keywords gdm
#'
#' @importFrom methods is
#'
#' @export
predict.gdm <- function(object, data, time=FALSE, predRasts=NULL, filename="", ...){
  #################
  ##lines used to quickly test function
  ##object = gdm model
  ##data = a sitepair table
  #object <- gdm.rastF
  #data <- envRast
  #time <- T
  #predRasts <- futRasts
  #################
  options(warn.FPU = FALSE)

  ##if making a time prediction, check if all data are in the correct format,
  ##and then transforms the raster data into data tables in order to utilize
  ##the C predict utility
  if (time) {
    .check_pkgs("terra")
    # checks to make sure the inputs are correct
    if (is.null(predRasts)) {
      stop("Prediction rasters required when time = TRUE")
    }
    if (!.is_raster(data)) {
      stop("Prediction data need to be a raster object when time = TRUE")
    }
    if (!.is_raster(predRasts)) {
      stop("predRasts need to be a raster object when time = TRUE")
    }
    data <- .check_rast(data, "data")
    predRasts <- .check_rast(predRasts, "predRasts")

    # tests to see if raster data is stackable
    tryCatch(
      {
        terra::compareGeom(data, predRasts, lyrs=TRUE, crs=TRUE, ext=TRUE, rowcol=TRUE)
      },
      error = function(cond) {
        message(
          "Current and prediction raster objects must have matching layer counts and identical extent, resolution, and origin to be stackable."
        )
      }
    )

    for(i in 1:terra::nlyr(data)){
      if(names(data)[i]!=names(predRasts)[i]){
        stop("Layer names do not match the variables used to fit the model.")
      }
    }
    if (object$geo) {
      if(terra::nlyr(data)!=length(object$predictors)-1 | terra::nlyr(predRasts)!=length(object$predictors)-1){
        stop("Number of variables supplied for prediction does not equal the number used to fit the model.")
      }
    } else {
      if(terra::nlyr(data)!=length(object$predictors) | terra::nlyr(predRasts)!=length(object$predictors)){
        stop("Number of variables supplied for prediction does not equal the number used to fit the model.")
      }
    }

    # create XY rasters; data and predRasts must have the same XY
    x <- terra::init(data, fun = "x")
    y <- terra::init(data, fun = "y")
    dummData <- terra::init(data[[1]], fun = 0L)

    # sets the correct names to the data
    names(data) <- paste0("s1.", names(data))
    names(predRasts) <- paste0("s2.", names(predRasts))

    # stack all the raster layers to for prediction
    data <- c(
      setNames(dummData, "distance"),
      setNames(dummData, "weights"),
      setNames(x, "s1.xCoord"),
      setNames(y, "s1.yCoord"),
      setNames(x, "s2.xCoord"),
      setNames(y, "s2.yCoord"),
      data,
      predRasts
    )

  }

  # makes the prediction based on the data object
  gdm_predict <- function(mod, dat, ...) {

    nr <- nrow(dat)
    predicted <- rep(0, times = nr)

    z <- .C( "GDM_PredictFromTable",
             as.matrix(dat),
             as.integer(mod$geo),
             as.integer(length(mod$predictors)),
             as.integer(nr),
             as.double(mod$knots),
             as.integer(mod$splines),
             as.double(c(mod$intercept, mod$coefficients)),
             preddata = as.double(predicted),
             PACKAGE = "gdm")

    return(
      z$preddata
    )
  }

  # if a time prediction, maps the predicted values to a raster and returns
  # the layer, otherwise returns a dataframe of the predicted values
  if(time){
    # predict using gdm model and terra package
    output <- terra::predict(
      object = data,
      model = object,
      fun = gdm_predict,
      na.rm = TRUE,
      filename = filename,
      ...
    )

    return(output)

  } else{
    # predict using a data.frame
    output <- gdm_predict(
      mod = object,
      dat = data
    )

    return(output)
  }
}

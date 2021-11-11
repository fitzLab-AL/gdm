#' Predict Biological Dissimilarities Between Sites or Times Using a Generalized Dissimilarity Model
#'
#' This function predicts biological distances between sites or times using a
#'  model object returned from \code{\link{gdm}}. Predictions between site
#'  pairs require a data frame containing the values of predictors for pairs
#'  of locations, formatted as follows: distance, weights, s1.X, s1.Y, s2.X,
#'  s2.Y, s1.Pred1, s1.Pred2, ..., s1.PredN, s2.Pred1, s2.Pred2, ..., s2.PredN, ...,
#'  Predictions of biological change through time require two raster stacks or
#'  bricks for environmental conditions at two time periods, each with a
#'  layer for each environmental predictor in the fitted model.
#'
#' @usage \method{predict}{gdm}(object, data, time=FALSE, predRasts=NULL, ...)
#'
#' @param object A gdm model object resulting from a call to \code{\link{gdm}}.
#'
#' @param data Either a data frame containing the values of predictors for pairs of sites, in the same format and structure as used to fit the model using \code{\link{gdm}} or a raster stack if a prediction of biological change through time is needed.
#'
#'For a data frame, the first two columns - distance and weights - are required by the function but are not used in the prediction and can therefore be filled with dummy data (e.g. all zeros). If geo is TRUE, then the s1.X, s1.Y and s2.X, s2.Y columns will be used for calculating the geographical distance between each site for inclusion of the geographic predictor term into the GDM model. If geo is FALSE, then the s1.X, s1.Y, s2.X and s2.Y data columns are ignored. However these columns are still REQUIRED and can be filled with dummy data (e.g. all zeroes). The remaining columns are for N predictors for Site 1 and followed by N predictors for Site 2. The order of the columns must match those in the site-pair table used to fit the model.
#'
#'A raster stack should be provided only when time=T and should contain one layer for each environmental predictor in the same order as the columns in the site-pair table used to fit the model.
#'
#' @param time TRUE/FALSE: Is the model prediction for biological change through time?
#'
#' @param predRasts A raster stack characterizing environmental conditions for a different time in the past or future, with the same extent, resolution, and layer order as the data object. Required only if time=T.
#'
#' @param ... Ignored.
#'
#' @return
#' predict returns either a response vector with the same length as the number of rows in the input data frame or a raster depicting change through time across the study region.
#'
#' @seealso \code{\link{gdm.transform}}
#'
#' @examples
#' ##sets up site-pair table
#' load(system.file("./data/gdm.RData", package="gdm"))
#' sppData <- gdmExpData[, c(1,2,14,13)]
#' envTab <- gdmExpData[, c(2:ncol(gdmExpData))]
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
#' rastFile <- system.file("./extdata/stackedVars.grd", package="gdm")
#' envRast <- stack(rastFile)
#'
#' ##make some fake climate change data
#' futRasts <- envRast
#' ##reduce winter precipitation by 25%
#' futRasts[[3]] <- futRasts[[3]]*0.75
#'
#' timePred <- predict(gdmMod, envRast, time=TRUE, predRasts=futRasts)
#' plot(timePred)
#'
#' @keywords gdm
#'
#' @export
##function to either predict the biological dissimilarities between sites,
##or to predict the dissimilarity of the same sites between two time periods,
##based on a gdm
predict.gdm <- function (object, data, time=FALSE, predRasts=NULL, ...){
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

  ##if making a time prediction, makes sure all data is in the correct format,
  ##and then transforms the raster data into data tables in order to utalize
  ##the C predict utility
  if(time==TRUE){
    ##checks to make sure the inputs are correct
    if(is.null(predRasts)==TRUE){
      stop("Prediction rasters required when time = TRUE")
    }
    if(!is(data, "RasterStack") & !is(data, "RasterLayer") & !is(data, "RasterBrick")){
      stop("Prediction data need to be a raster object when time = TRUE")
    }
    if(!is(predRasts, "RasterStack") & !is(predRasts, "RasterLayer") & !is(predRasts, "RasterBrick")){
      stop("predRasts need to be a raster object when time = TRUE")
    }
    if(raster::nlayers(data)!=raster::nlayers(predRasts)){
      stop("Current and future raster objects must have the same number of layers")
    }
    if(object$geo==TRUE){
      if(raster::nlayers(data)!=length(object$predictors)-1 | raster::nlayers(predRasts)!=length(object$predictors)-1){
        stop("Number of variables supplied for prediction does not equal the number used to fit the model.")
      }
    }else{
      if(raster::nlayers(data)!=length(object$predictors) | raster::nlayers(predRasts)!=length(object$predictors)){
        stop("Number of variables supplied for prediction does not equal the number used to fit the model.")
      }
    }

    for(i in 1:raster::nlayers(data)){
      if(names(data)[i]!=names(predRasts)[i]){
        stop("Layer names do not match the variables used to fit the model.")
      }
    }
    ##tests to see if raster data is stackable
    tryRasts <- try(raster::stack(data[[1]], predRasts[[1]]), silent=TRUE)
    if(is(tryRasts, "try-error")){
      stop("Current and prediction rasters differ in extent, resolution, and / or origin and therefore are not stackable.")
    }

    ##sets up sitepair table with current and future data
    predLayer <- data[[1]]
    currXY <- as.data.frame(na.omit(rasterToPoints(data, progress='text')))
    predXY <- as.data.frame(na.omit(rasterToPoints(predRasts, progress='text')))
    cells <- raster::cellFromXY(predLayer, cbind(currXY$x, currXY$y))
    dummData <- rep.int(0, nrow(currXY))
    data <- cbind(dummData, dummData, currXY[,1:2], currXY, predXY[,-c(1,2)])
    ##adds s1 or s2 to the variables name of the data
    t1var <- paste("s1.", colnames(currXY)[-c(1,2)], sep="")
    t2var <- paste("s2.", colnames(predXY)[-c(1,2)], sep="")
    ##sets the correct names to the data
    colnames(data) <- c("distance", "weights", "s1.xCoord", "s1.yCoord",
                        "s2.xCoord", "s2.yCoord", t1var, t2var)
  }

  ##makes the prediction based on the data object
  predicted <- rep(0,times=nrow(data))
  z <- .C( "GDM_PredictFromTable",
           as.matrix(data),
           as.integer(object$geo),
           as.integer(length(object$predictors)),
           as.integer(nrow(data)),
           as.double(object$knots),
           as.integer(object$splines),
           as.double(c(object$intercept,object$coefficients)),
           preddata = as.double(predicted),
           PACKAGE = "gdm")

  ##if a time prediction, maps the predicted values to a raster and returns
  ##the layer, otherwise returns a dataframe of the predicted values
  if(time==FALSE){
    return(z$preddata)
  }else{
    predLayer[cells] <- z$preddata
    return(predLayer)
  }
}

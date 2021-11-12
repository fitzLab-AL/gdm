#' Fit a Generalized Dissimilarity Model to Tabular Site-Pair Data
#'
#' @description The gdm function is used to fit a generalized dissimilarity model to tabular
#' site-pair data formatted as follows using the \code{\link[gdm]{formatsitepair}}
#' function: distance, weights, s1.xCoord, s1.yCoord, s2.xCoord, s2.yCoord,
#' s1.Pred1, s1.Pred2, ...,s1.PredN, s2.Pred1, s2.Pred2, ..., s2.PredN. The
#' distance column contains the response variable must be any ratio-based
#' dissimilarity (distance) measure between Site 1 and Site 2. The weights column
#' defines any weighting to be applied during fitting of the model. If equal
#' weighting is required, then all entries in this column should be set to 1.0
#' (default). The third and fourth columns, s1.xCoord and s1.yCoord, represent
#' the spatial coordinates of the first site in the site pair (s1). The fifth
#' and sixth columns, s2.xCoord and s2.yCoord, represent the coordinates of the
#' second site (s2). Note that the first six columns are REQUIRED, even if you
#' do not intend to use geographic distance as a predictor (in which case these
#' columns can be loaded with dummy data if the actual coordinates are
#' unknown - though that would be weird, no?). The next N*2 columns contain values
#' for N predictors for Site 1, followed by values for the same N predictors for
#' Site 2. \cr
#'
#' The following is an example of a GDM input table header with three
#' environmental predictors (Temp, Rain, Bedrock): \cr
#'
#' distance, weights, s1.xCoord, s1.yCoord, s2.xCoord, s2.yCoord, s1.Temp,
#' s1.Rain, s1.Bedrock, s2.Temp, s2.Rain, s2.Bedrock
#'
#' @usage gdm(data, geo=FALSE, splines=NULL, knots=NULL)
#'
#' @param data A data frame containing the site pairs to be used to fit the GDM
#'   (obtained using the \code{\link[gdm]{formatsitepair}} function). The
#'   observed response data must be located in the first column. The weights to
#'   be applied to each site pair must be located in the second column. If geo
#'   is TRUE, then the s1.xCoord, s1.yCoord and s2.xCoord, s2.yCoord columns
#'   will be used to calculate the geographic distance between site pairs for
#'   inclusion as the geographic predictor term in the model. Site coordinates
#'   ideally should be in a projected coordinate system (i.e., not longitude-latitude)
#'   to ensure proper calculation of geographic distances. If geo is FALSE
#'   (default), then the s1.xCoord, s1.yCoord, s2.xCoord and s2.yCoord data
#'   columns must still be included, but are ignored in fitting the model.
#'   Columns containing the predictor data for Site 1, and the predictor data
#'   for Site 2, follow.
#'
#' @param geo Set to TRUE if geographic distance between sites is to be included
#'   as a model term. Set to FALSE if geographic distance is to be omitted from
#'   the model. Default is FALSE.
#'
#' @param splines An optional vector of the number of I-spline basis functions
#'   to be used for each predictor in fitting the model. If supplied, it must
#'   have the same length as the number of predictors (including geographic
#'   distance if geo is TRUE). If this vector is not provided (splines=NULL),
#'   then a default of 3 basis functions is used for all predictors.
#'
#' @param knots An optional vector of knots in  \emph{units of the predictor
#'   variables} to be used in the fitting process. If knots are supplied and
#'   splines=NULL, then the knots argument must have the same length as the
#'   number of predictors * n, where n is the number of knots (default=3). If both
#'   knots and the number of splines are supplied, then the length of the knots
#'   argument must be the same as the sum  of the values in the splines vector.
#'   Note that the default values for knots when the default three I-spline basis
#'   functions are 0 (minimum), 50 (median), and 100 (maximum) quantiles.
#'
#' @return gdm returns a gdm model object. The function
#'   \code{\link[gdm]{summary.gdm}} can be used to obtain or print a synopsis of the
#'   results. A gdm model object is a list containing at least the following
#'   components:
#'
#'   \describe{ \item{dataname}{The name of the table used as the data argument
#'   to the model.} \item{geo}{Whether geographic distance was used as a
#'   predictor in the model.} \item{gdmdeviance}{The deviance of the fitted GDM
#'   model.} \item{nulldeviance}{ The deviance of the null model.}
#'   \item{explained}{The percentage of null deviance explained by the fitted
#'   GDM model.} \item{intercept}{The fitted value for the intercept term in the
#'   model.} \item{predictors}{A list of the names of the predictors that were
#'   used to fit the model.} \item{coefficients}{A list of the coefficients for
#'   each spline for each of the predictors considered in model fitting.}
#'   \item{knots}{A vector of the knots derived from the x data (or user
#'   defined), for each predictor.} \item{splines}{ A vector of the number of
#'   I-spline basis functions used for each predictor.} \item{creationdate}{ The
#'   date and time of model creation.} \item{observed}{The observed response for
#'   each site pair (from data column 1).} \item{predicted}{ The predicted
#'   response for each site pair, from the fitted model (after applying the link
#'   function).} \item{ecological}{ The linear predictor (ecological distance)
#'   for each site pair, from the fitted model (before applying the link
#'   function).} }
#'
#' @references Ferrier S, Manion G, Elith J, Richardson, K (2007) Using
#'   generalized dissimilarity modelling to analyse and predict patterns of beta
#'   diversity in regional biodiversity assessment. \emph{Diversity &
#'   Distributions} 13, 252-264.
#'
#' @seealso \code{\link[gdm]{formatsitepair}, \link[gdm]{summary.gdm},
#'   \link[gdm]{plot.gdm}, \link[gdm]{predict.gdm}, \link[gdm]{gdm.transform}}
#'
#' @examples
#'  ##fit table environmental data
#'  ##sets up site-pair table, environmental tabular data
#'  load(system.file("./data/southwest.RData", package="gdm"))
#'  sppData <- southwest[c(1,2,13,14)]
#'  envTab <- southwest[c(2:ncol(southwest))]
#'
#'  sitePairTab <- formatsitepair(sppData, 2, XColumn="Long", YColumn="Lat", sppColumn="species",
#'                                siteColumn="site", predData=envTab)
#'
#'  ##fit table GDM
#'  gdmTabMod <- gdm(sitePairTab, geo=TRUE)
#'  summary(gdmTabMod)
#'
#'  ##fit raster environmental data
#'  ##sets up site-pair table
#'  rastFile <- system.file("./extdata/stackedVars.grd", package="gdm")
#'  envRast <- raster::stack(rastFile)
#'
#'  ##environmental raster data
#'  sitePairRast <- formatsitepair(sppData, 2, XColumn="Long",
#'                                 YColumn="Lat", sppColumn="species",
#'                                 siteColumn="site", predData=envRast)
#'  ##sometimes raster data returns NA in the site-pair table, these rows will
#'  ##have to be removed before fitting gdm sitePairRast <-
#'  na.omit(sitePairRast)
#'
#'  ##fit raster GDM
#'  gdmRastMod <- gdm(sitePairRast, geo=TRUE)
#'  summary(gdmRastMod)
#'
#' @keywords gdm
#'
#' @importFrom stats median quantile
#'
#' @export
# fit a gdm object from a sitepair table
gdm <- function (data, geo=FALSE, splines=NULL, knots=NULL){
  #################
  ##lines used to quickly test function
  #data <- sitepairs
  #geo <- F
  #splines <- NULL
  # knots <- NULL
  #################
  options(warn.FPU = FALSE)

  ##adds error checking to gdm function
  ##checks to see if in site-pair format from formatsitepair function
  if(!is(data, "gdmData")){
    warning("Site-pair table provided for model fitting is not of class 'gdmData'. Did you format your data using formatsitepair?")
  }
  ##checks to makes sure data is a matrix or data frame
  if(!(is(data, "gdmData") | is(data, "matrix") | is(data, "data.frame"))){
    stop("Site-pair table provided for model fitting needs to be either of class 'gdmData', 'matrix', or 'data.frame'.")
  }

  ##sanity check on the data table
  if(ncol(data) < 6){
    stop("Site-pair table provided for model fitting requires at least 6 columns: Observed, weights, s1.xCoord, s1.yCoord, s2.xCoord, s2.yCoord")
  }
  if(nrow(data) < 1){
    stop("Site-pair table provided for model fitting has < 1 row.")
  }

  ##checks that geo has either TRUE or FALSE
  if(!(geo==TRUE | geo==FALSE)){
    stop("geo argument must be either TRUE or FALSE")
  }
  ##makes sure splines is a numeric vector
  if(is.null(splines)==FALSE & !is(splines, "numeric")){
    stop("splines argument is not of class = 'numeric'.")
  }
  ##checks knots inputs
  if(is.null(knots)==FALSE & !is(knots, "numeric")){
    stop("knots argument is not of class = 'numeric'.")
  }

  ##check that the response data is [0..1]
  rtmp <- data[,1]
  if(length(rtmp[rtmp<0]) > 0){
    stop("Response data have negative values. Must be between 0 - 1.")
  }
  if (length(rtmp[rtmp>1]) > 0){
    stop("Response data have values greater than 1. Must be between 0 - 1.")
  }

  ##current data format is response,weights,X0,Y0,X1,Y1 before any predictor data (thus 6 leading columns)
  LEADING_COLUMNS <- 6
  if(geo){
    nPreds <- ( ncol(data) - LEADING_COLUMNS ) / 2 + 1
  }else{
    nPreds <- ( ncol(data) - LEADING_COLUMNS ) / 2
  }

  ##checks to make sure at least one predictor is available
  if(nPreds < 1){
    stop("Data has no predictor varaibles.")
  }

  ##setup the predictor name list, and removes the "s1." and "s2." to make resulting names more intuitive
  if(geo==TRUE){
    if(nPreds > 1){
      predlist <- c("Geographic", sapply(strsplit(names(data)[(LEADING_COLUMNS+1):(LEADING_COLUMNS+nPreds-1)], "s1.", fixed=T), "[[", 2))
    }else{
      predlist <- c("Geographic")
    }
  }else{
    predlist <- sapply(strsplit(names(data)[(LEADING_COLUMNS+1):(LEADING_COLUMNS+nPreds)], "s1.", fixed=T), "[[", 2)
  }

  ##deal with the splines and knots
  if(is.null(knots)){
    ##generate knots internally from the data
    if(is.null(splines)){
      nSplines <- 3
      quantvec <- rep(0, nPreds * nSplines)
      splinvec <- rep(nSplines, nPreds)

      if(geo==TRUE){
        ##get knots for the geographic distance
        v <- sqrt((data[,3]-data[,5])^2 + (data[,4]-data[,6])^2)
        quantvec[1] <- min(v)
        quantvec[2] <- stats::median(v)
        quantvec[3] <- max(v)

        if(nPreds > 1){
          ##get knots for the environmental predictors
          for(i in seq(from = 1, to = nPreds-1, by = 1)){
            v <- c(data[,i+LEADING_COLUMNS], data[,i+LEADING_COLUMNS+nPreds-1])
            index = i * nSplines
            quantvec[index+1] <- min(v)
            quantvec[index+2] <- stats::median(v)
            quantvec[index+3] <- max(v)
          }
        }
      }else{
        ## get knots for the environmental predictors after skipping geographic preds
        for(i in seq(from = 1, to = nPreds, by = 1)){
          v <- c(data[,i+LEADING_COLUMNS], data[,i+LEADING_COLUMNS+nPreds])
          index = (i-1) * nSplines
          quantvec[index+1] <- min(v)
          quantvec[index+2] <- stats::median(v)
          quantvec[index+3] <- max(v)
        }
      }
    }else{
      ##otherwise check that the supplied splines vector has enough data and minumum spline values of 3
      if(length(splines) != nPreds){
        stop(paste("Number of splines does not equal the number of predictors.
              Splines argument has", length(splines), "items but needs", nPreds, "items."))
      }

      ##count the total number of user defined splines to dimension the knots vector
      quantvec <- rep(0, sum(splines))
      splinvec <- splines

      if(geo==T){
        if(splines[1] < 3){
          stop("Must have at least 3 splines per predictor.")
        }

        ## get knots for the geographic distance
        v <- sqrt((data[,3]-data[,5])^2 + (data[,4]-data[,6])^2)
        quantvec[1] <- min(v)		## 0% knot
        quantvec[splines[1]] <- max(v)	## 100% knot

        quant_increment <- 1.0 / (splines[1]-1)
        this_increment <- 1
        for (i in seq(from = 2, to = (splines[1]-1), by = 1)){
          ## mid % knots
          quantvec[i] <- stats::quantile(v,quant_increment*this_increment)
          this_increment <- this_increment + 1
        }

        if(nPreds > 1){
          ##get knots for the environmental predictors
          current_quant_index <- splines[1]
          for(i in seq(from = 1, to = nPreds-1, by = 1)){
            num_splines <- splines[i+1]
            if(num_splines < 3){
              stop("Must have at least 3 splines per predictor.")
            }

            v <- c(data[,i+LEADING_COLUMNS], data[,i+LEADING_COLUMNS+nPreds-1])
            quantvec[current_quant_index+1] <- min(v)	            ## 0% knot
            quantvec[current_quant_index+num_splines] <- max(v)	    ## 100% knot

            quant_increment <- 1.0 / (num_splines-1)
            this_increment <- 1
            for(i in seq(from = 2, to = (num_splines-1), by = 1)){
              ##mid % knots
              quantvec[current_quant_index+i] <- stats::quantile(v,quant_increment*this_increment)
              this_increment <- this_increment + 1
            }

            current_quant_index <- current_quant_index + num_splines
          }
        }
      }else{
        ##get knots for the environmental predictors
        current_quant_index <- 0
        for(i in seq(from = 1, to = nPreds, by = 1)){
          num_splines <- splines[i]
          if(num_splines < 3){
            stop("Must have at least 3 splines per predictor.")
          }

          v <- c(data[,i+LEADING_COLUMNS], data[,i+LEADING_COLUMNS+nPreds])
          quantvec[current_quant_index+1] <- min(v)	        ## 0% knot
          quantvec[current_quant_index+num_splines] <- max(v)	## 100% knot

          quant_increment <- 1.0 / (num_splines-1)
          this_increment <- 1
          for(i in seq(from = 2, to = (num_splines-1), by = 1)){
            ##mid % knots
            quantvec[current_quant_index+i] <- stats::quantile(v,quant_increment*this_increment)
            this_increment <- this_increment + 1
          }
          current_quant_index <- current_quant_index + num_splines
        }
      }
    }
  }else{
    ##user defined knots supplied as an argument
    if(is.null(splines)){
      ##check that there are nPreds * 3 knots in the user defined vector
      if(length(knots) != (nPreds * 3)){
        stop(paste("When knots are supplied by the user, there should be", (nPreds * 3), "items in the knots argument, not", length(knots), "items."))
      }

      ## now check that each of the three knots for each predictor are in ascending order
      for(i in seq(from = 1, to = nPreds, by = 1)){
        index = i * 3
        if((knots[index-1] < knots[index-2]) ||
           (knots[index] < knots[index-2]) ||
           (knots[index] < knots[index-1])){
          stop(paste("Knots for ", predlist[i], "are not in ascending order."))
        }
      }

      nSplines <- 3
      quantvec <- knots
      splinvec <- rep(nSplines, nPreds)
    }else{
      ##check that there are sum(splines) knots in the user defined vector
      if(length(knots) != sum(splines)){
        stop(paste("When knots are supplied by the user, there should be", sum(splines), "items in the knots argument, not", length(knots), "items."))
      }

      ##now check that each of the knots for each predictor are in ascending order
      index = 0
      for(i in seq(from = 1, to = nPreds, by = 1)){
        for(j in seq(from = 2, to = splines[i], by = 1)){
          if(knots[index+j] < knots[index+j-1]){
            stop(paste("Knots for ", predlist[i], "are not in ascending order."))
          }
        }
        index <- index + splines[i]
      }

      quantvec <- knots
      splinvec <- splines
    }
  }

  p1 <- 0
  p2 <- 0
  p3 <- 0
  p4 <- 0
  p5 <- rep(0,times=length(quantvec))
  p6 <- rep(0,times=nrow(data))
  p7 <- rep(0,times=nrow(data))
  p8 <- rep(0,times=nrow(data))

  ##Call the dll function, fitting the gdm model
  z <- .C( "GDM_FitFromTable",
           tempfile("gdmtmp_", tmpdir=getwd()), #paste(getwd()),
           as.matrix(data),
           as.integer(geo),
           as.integer(nPreds),
           as.integer(nrow(data)),
           as.integer(ncol(data)),
           as.integer(splinvec),
           as.double(quantvec),
           gdmdev = as.double(p1),
           nulldev = as.double(p2),
           expdev = as.double(p3),
           intercept = as.double(p4),
           coeffs = as.double(p5),
           response = as.double(p6),
           preddata = as.double(p7),
           ecodist = as.double(p8),
           PACKAGE = "gdm")

  m <- match.call(expand.dots = F)

  ##creates the gdm object, and fills its parts
  gdmModOb <- structure(list(dataname = m[[2]],
                             geo = geo,
                             sample = nrow(data),
                             gdmdeviance = z$gdmdev,
                             nulldeviance = z$nulldev,
                             explained = z$expdev,
                             intercept = z$intercept,
                             predictors = predlist,
                             coefficients = z$coeffs,
                             knots = quantvec,
                             splines = splinvec,
                             creationdate = date(),
                             observed = z$response,
                             predicted = z$preddata,
                             ecological = z$ecodist))
  ##sets gdm object class
  class(gdmModOb) <- c("gdm", "list")

  ##reports a warning should the model "fit", yet the sum of coefficients = 0
  if(sum(gdmModOb$coefficients)==0){
    warning("The algorithm was unable to fit a model to your data. The sum of the spline coefficients = 0 and deviance explained = NULL. Returning NULL object.")
    ##sets the deviance explained to NULL, to reflect that the model didn't fit correctly
    gdmModOb <- NULL
  }

  ##returns gdm object
  return(gdmModOb)
}

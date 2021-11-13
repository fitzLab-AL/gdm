#' @title Extract I-spline Values From a gdm Object.
#'
#' @description Extracts the I-spline values from a gdm object. There is one
#' I-spline for each predictor that has at least one non-zero coefficient in
#' the fitted model.
#'
#' @usage isplineExtract(model)
#'
#' @param model A gdm object from \code{\link[gdm]{gdm}}.
#'
#' @return A list with two items. The first item contains the x-values (actual
#' values of the predictors) of the I-splines and the second item contains the
#' y-values (partial ecological distances) of the fitted I-splines.
#'
#' @references Ferrier S, Manion G, Elith J, Richardson, K (2007) Using
#' generalized dissimilarity modelling to analyse and predict patterns of beta
#' diversity in regional biodiversity assessment. \emph{Diversity & Distributions}
#' 13, 252-264.
#'
#' Fitzpatrick MC, Sanders NJ, Normand S, Svenning J-C, Ferrier S, Gove AD,
#' Dunn RR (2013). Environmental and historical imprints on beta diversity: insights
#' from variation in rates of species turnover along gradients. Proceedings of the
#' Royal Society: Series B 280, art. 1768
#'
#' @examples
#' ##set up site-pair table using the southwest data set
#' sppData <- southwest[, c(1,2,14,13)]
#' envTab <- southwest[, c(2:ncol(southwest))]
#' sitePairTab <- formatsitepair(sppData, 2, XColumn="Long", YColumn="Lat", sppColumn="species",
#'                               siteColumn="site", predData=envTab)
#'
#' ##create GDM
#' gdmMod <- gdm(sitePairTab, geo=TRUE)
#'
#' ##extracts splines
#' exSplines <- isplineExtract(gdmMod)
#'
#' ##plot spline(s)
#' #spline for winter precip (bio19)
#' plot(exSplines[[1]][,"bio19"], exSplines[[2]][,"bio19"], type="l",
#'      lwd=3, xlab="Winter precipitation (mm)", ylab="Partial Ecological Distance")
#'
#' @keywords gdm
#'
#' @export
isplineExtract <- function (model){
  ###########################
  #model = gdmOb
  ###########################
  ##error checking
  ##checks to make sure a gdm model is given
  if(!is(model, "gdm")){
    stop("model argument must be of class = 'gdm'.")
  }

  ##Collects or sets simple data
  options(warn.FPU = FALSE)
  PSAMPLE <- 200
  preddata <- rep(0, times = PSAMPLE)
  pn <- model$predictors
  nPreds <- length(pn)
  yDat <- xDat <- matrix(0,PSAMPLE,nPreds)
  colnames(yDat) <- colnames(xDat) <- pn
  pmin <- 1
  pmax <- PSAMPLE

  ##cycles through each predictor and fills the spline matrices
  splineindex <- 1
  for (i in 1:nPreds){
    #i<-1
    numsplines <- model$splines[i]
    z <- .C("GetPredictorPlotData",
            pdata = as.double(preddata),
            as.integer(PSAMPLE),
            as.double(model$coefficients[splineindex:(splineindex + numsplines - 1)]),
            as.double(model$knots[splineindex:(splineindex + numsplines - 1)]),
            as.integer(numsplines),
            PACKAGE = "gdm")
    yDat[,i] <- z$pdata
    pmin <- pmin + PSAMPLE
    pmax <- pmax + PSAMPLE
    xDat[,i] <-  seq(from=model$knots[[(i*numsplines)-(numsplines-1)]],to=model$knots[[(i*numsplines)]], length=PSAMPLE)
    splineindex <- splineindex + numsplines
  }

  ##lists and returns matrices
  outData <- list(x=xDat,y=yDat)
  return(outData)
}

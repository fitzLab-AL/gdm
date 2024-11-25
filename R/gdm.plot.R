#' @title Plot Model Fit and I-splines from a Fitted Generalized Dissimilarity
#'  Model.
#'
#' @description plot is used to plot the I-splines and fit of a generalized
#' dissimilarity model created using the \code{\link[gdm]{gdm}} function.
#'
#' @usage \method{plot}{gdm}(x, plot.layout = c(2, 2), plot.color = "blue",
#'   plot.linewidth = 2, include.rug = FALSE, rug.sitepair = NULL, ...)
#'
#' @param x A gdm model object returned from \code{\link[gdm]{gdm}}.
#'
#' @param plot.layout This argument specifies the row and column layout for the
#' plots, including: (1) a single page plot of observed response data against
#' the raw linear predictor (ecological distance) from the model, and (2) a
#' single page plot of the observed response against the predicted response
#' from the model, i.e. after applying the link function, 1.0 - exp(-y), to the
#' linear predictor, and (3) the I-splines fitted to the individual predictors.
#' Default is 2 rows by 2 columns. To produce one predictor plot per page set
#' plot.layout to c(1,1). The first two model plots are always produced on a
#' single page each and therefore the layout parameter affects only the layout
#' of the I-spline plots for those predictors that featured in the model fitting
#' process (i.e., predictors with all-zero I-spline coefficients are not plotted).
#'
#' @param plot.color Color of the data points that are plotted for the overall plots.
#'
#' @param plot.linewidth The line width for the regression line over-plotted in
#' the two overall plots to optimize the display of the line over the data points.
#'
#' @param include.rug Whether or not to include a rug plot of the predictor
#' values used to fit the gdm in the I-spline plots. When set to TRUE, a s
#' ite-pair table must be supplied for the rug.sitepair argument. Default is FALSE.
#'
#' @param rug.sitepair A site-pair table used to add a rug plot of the predictor
#' values used to fit the gdm in the I-spline plots. This should be the same
#' site-pair table used to fit the gdm model being plotted. The function does
#' not check whether the supplied site-pair table matches that used in model fitting.
#'
#' @param ... Ignored.
#'
#' @return plot returns NULL. Use \code{\link[gdm]{summary.gdm}} to obtain a
#' synopsis of the model object.
#'
#' @references Ferrier S, Manion G, Elith J, Richardson, K (2007) Using
#' generalized dissimilarity modelling to analyse and predict patterns of
#' beta diversity in regional biodiversity assessment.
#' \emph{Diversity & Distributions} 13:252-264.
#'
#' @seealso \code{\link[gdm]{isplineExtract}}
#'
#' @examples
#' ##set up site-pair table using the southwest data set
#' sppData <- southwest[, c(1,2,13,14)]
#' envTab <- southwest[, c(2:ncol(southwest))]
#' sitePairTab <- formatsitepair(sppData, 2, XColumn="Long", YColumn="Lat",
#'                               sppColumn="species", siteColumn="site",
#'                               predData=envTab)
#'
#' ##create GDM
#' gdmMod <- gdm(sitePairTab, geo=TRUE)
#'
#' ##plot GDM
#' plot(gdmMod, plot.layout=c(3,3))
#'
#' @keywords gdm
#'
#' @importFrom stats na.omit
#' @importFrom grDevices dev.new
#' @importFrom grDevices dev.next
#' @importFrom graphics lines
#' @importFrom graphics points
#' @importFrom graphics rug
#' @importFrom graphics par
#'
#' @export
##function to plot the splines of a gdm object
plot.gdm <- function (x, plot.layout = c(2,2), plot.color = "blue",
                      plot.linewidth=2.0, include.rug=FALSE, rug.sitepair=NULL, ...){
  #################
  ##lines used to quickly test function
  #x <- model.2
  #plot.layout <- c(2,2)
  #plot.color <- "green"
  #plot.linewidth <- 2.0
  #include.rug <- T
  #rug.sitepair <- GDM_input_table
  #################
  ##checks to make sure that a site-pair table has been included

  options(warn.FPU = FALSE)
  PSAMPLE <- 200
  preddata <- rep(0,times=PSAMPLE)

  ##establish what plot layout to use
  thisplot <- 0
  one_page_per_plot <- FALSE
  if ((plot.layout[1]==1) && (plot.layout[2]==1)){
    one_page_per_plot <- TRUE
  }else{
    par(mfrow=plot.layout)
  }

  ##plots the compositional dissimilarity spline plot
  plot(x$ecological, x$observed, xlab="Predicted Ecological Distance", ylab="Observed Compositional Dissimilarity", type="n", ylim=c(0,1))
  points(x$ecological, x$observed, pch=20, cex=0.25, col=plot.color)
  overlayX <- seq( from=min(x$ecological), to=max(x$ecological), length=PSAMPLE )
  overlayY <- 1 - exp( - overlayX )
  lines( overlayX, overlayY, lwd=plot.linewidth )
  thisplot <- thisplot + 1

  ##determines rather or not to put multiple plots on one page or not
  if(one_page_per_plot){
    dev.new()
    dev.next()
  }
  ##plots the second compositional dissimilarity spline plot
  plot(x$predicted, x$observed, xlab="Predicted Compositional Dissimilarity", ylab="Observed Compositional Dissimilarity", type="n", ylim=c(0,1))
  points( x$predicted, x$observed, pch=20, cex=0.25, col=plot.color )
  overlayX <- overlayY <- seq( from=min(x$predicted), to=max(x$predicted), length=PSAMPLE )
  lines( overlayX, overlayY, lwd=plot.linewidth )
  thisplot <- thisplot + 1

  #########################
  # reorder the predictors, splines, coeffs in order of
  # importance based on sum(coeffs)
  thiscoeff <- 1
  thisquant <- 1
  sumCoeff <- NULL
  for(i in 1:length(x$predictors)){
    numsplines <- x$splines[[i]]
    holdCoeff <- NULL
    for(j in 1:numsplines){
      holdCoeff[j] <- x$coefficients[[thiscoeff]]
      thiscoeff <- thiscoeff + 1
    }
    sumCoeff[i] <- sum(holdCoeff)
  }

  lateralus <- NULL
  schism <- NULL
  orderPreds <- order(sumCoeff, decreasing = T)
  for(op in orderPreds){
    parabol <- 1+(cumsum(x$splines)[op]-x$splines[op])#1+(op*x$splines[op]-(x$splines[op]))
    parabola <- cumsum(x$splines)[op]#op*x$splines[op]
    lateralus <- c(lateralus, x$coefficients[parabol:parabola])
    schism <- c(schism, x$knots[parabol:parabola])
  }

  x$predictors <- x$predictors[orderPreds]
  x$splines <- x$splines[orderPreds]
  x$coefficients <- lateralus
  x$knots <- schism
  x$sumCoeff <- sumCoeff[orderPreds]
  #########################

  ##determine the max of all the predictor data, to be used in the plotting below
  preds <- length(x$predictors)
  predmax <- 0
  splineindex <- 1
  for(i in 1:preds){
    ##only if the sum of the coefficients associated with this predictor is > 0.....
    numsplines <- x$splines[i]
    if(sum(x$coefficients[splineindex:(splineindex+numsplines-1)]) > 0){
      ## get predictor plot Y-data
      z <- .C( "GetPredictorPlotData",
               pdata = as.double(preddata),
               as.integer(PSAMPLE),
               as.double(x$coefficients[splineindex:(splineindex+numsplines-1)]),
               as.double(x$knots[splineindex:(splineindex+numsplines-1)]),
               as.integer(numsplines),
               PACKAGE = "gdm" )

      v <- max(z$pdata)
      if(v > predmax){
        predmax <- v
      }
    }
    ##update the spline index
    splineindex <- splineindex + numsplines
  }

  ##plot the predictors with non-zero sum of coefficients
  splineindex <- 1
  for(i in 1:preds){
    #i <- 3
    ##only if the sum of the coefficients associated with this predictor is > 0.....
    numsplines <- x$splines[i]
    if(sum(x$coefficients[splineindex:(splineindex+numsplines-1)]) > 0){
      ##plots one graph per page, unless specified otherwise
      if (one_page_per_plot){
        dev.new()
        dev.next()
      }else{
        thisplot <- thisplot + 1
        if(thisplot > (plot.layout[1] * plot.layout[2])){
          par(ask=T)
          #dev.new()
          #x11()
          #dev.next()
          thisplot <- 1
          par(mfrow=plot.layout)
        }
      }

      ##get predictor plot Y-data
      z <- .C( "GetPredictorPlotData",
               pdata = as.double(preddata),
               as.integer(PSAMPLE),
               as.double(x$coefficients[splineindex:(splineindex+numsplines-1)]),
               as.double(x$knots[splineindex:(splineindex+numsplines-1)]),
               as.integer(numsplines),
               PACKAGE = "gdm")

      if(x$geo & x$predictors[i]=="Geographic"){
        varNam <- "Geographic Distance"
        ##calculates rug plot data
        if(include.rug==TRUE){
          rugData <- unique(sqrt(((rug.sitepair$s1.xCoord-rug.sitepair$s2.xCoord)^2)+((rug.sitepair$s1.yCoord-rug.sitepair$s2.yCoord)^2)))
        }
      }else{
        varNam <- x$predictors[i]
        ##gets rug plot data
        if(include.rug==TRUE){
          varDat <- grep(varNam, colnames(rug.sitepair))
          rugData <- unique(c(rug.sitepair[,c(varDat[1])], rug.sitepair[,c(varDat[2])]))
        }
      }

      parabol <- 1+(cumsum(x$splines)[i]-x$splines[i])
      parabola <- cumsum(x$splines)[i]

      plot(seq(from=x$knots[parabol],
               to=x$knots[parabola], length=PSAMPLE),
           z$pdata,
           xlab=varNam,
           ylab=paste("f(", varNam, ")", sep="" ),
           ylim=c(0,predmax),
           type="l")

      #plot(seq(from=x$knots[[(i*3)-2]], to=x$knots[[(i*3)]], length=PSAMPLE), z$pdata,
      #     xlab=varNam, ylab=paste("f(", varNam, ")", sep="" ), ylim=c(0,predmax), type="l")
      if(include.rug==TRUE){
        rug(rugData)
      }
    }
    ##update the spline index
    splineindex <- splineindex + numsplines
  }
}


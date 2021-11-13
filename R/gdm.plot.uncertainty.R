#' @title Plot I-splines with error bands using bootstrapping.
#'
#' @description This function estimates uncertainty in the fitted I-splines
#' by fitting many GDMs using a subsample of the data. The function can run in parallel
#' on multicore machines to reduce computation time (recommended for large number
#' of iterations). I-spline plots with error bands (+/- one standard deviation)
#' are produced showing (1) the variance of I-spline coefficients and (2) a rug plot
#' indicating how sites used in model fitting are distributed along each gradient.
#' Function result optionally can be saved to disk as a csv for custom plotting, etc.
#' The result output table will have 6 columns per predictor, three each for the
#' x and y values containing the lower bound, full model, and upper bound.
#'
#' @usage plotUncertainty(spTable, sampleSites, bsIters, geo=FALSE,
#' splines=NULL, knots=NULL, splineCol="blue", errCol="grey80",
#' plot.linewidth=2.0, plot.layout=c(2,2), parallel=FALSE, cores=2, save=FALSE,
#' fileName="gdm.plotUncertainy.csv")
#'
#' @param spTable A site-pair table, same as used to fit a \code{\link[gdm]{gdm}}.
#'
#' @param sampleSites The fraction (0-1) of sites to retain from the full
#' site-pair table when subsampling.
#'
#' @param bsIters The number of bootstrap iterations to perform.
#'
#' @param geo Same as the \code{\link[gdm]{gdm}} geo argument.
#'
#' @param splines Same as the \code{\link[gdm]{gdm}} splines argument.
#'
#' @param knots Same as the \code{\link[gdm]{gdm}} knots argument.
#'
#' @param splineCol The color of the plotted mean spline. The default is "blue".
#'
#' @param errCol The color of shading for the error bands (+/- one standard
#' deviation around the mean line). The default is "grey80".
#'
#' @param plot.linewidth The line width of the plotted mean spline line. The
#' default is 2.
#'
#' @param plot.layout Same as the \code{\link[gdm]{plot.gdm}} plot.layout argument.
#'
#' @param parallel Perform the uncertainty assessment using multiple
#' cores? Default = FALSE.
#'
#' @param cores When the parallel argument is set to TRUE, the number of
#' cores to be registered for the foreach loop. Must be <= the number of cores
#' in the machine running the function.
#'
#' @param save Save the function result (e.g., for custom plotting)? Default=FALSE.
#'
#' @param fileName Name of the csv file to save the data frame that contains the function
#' result. Default = gdm.plotUncertainy.csv. Ignored if save=FALSE.
#'
#' @return plotUncertainty returns NULL. Saves a csv to disk if save=TRUE.
#'
#' @references Shryock, D. F., C. A. Havrilla, L. A. DeFalco, T. C. Esque,
#' N. A. Custer, and T. E. Wood. 2015. Landscape genomics of \emph{Sphaeralcea ambigua}
#' in the Mojave Desert: a multivariate, spatially-explicit approach to guide
#' ecological restoration. \emph{Conservation Genetics} 16:1303-1317.
#'
#' @seealso \code{\link[gdm]{plot.gdm}, \link[gdm]{formatsitepair}, \link[gdm]{subsample.sitepair}}
#'
#' @examples
#' ##set up site-pair table using the southwest data set
#' sppData <- southwest[c(1,2,13,14)]
#' envTab <- southwest[c(2:ncol(southwest))]
#' sitePairTab <- formatsitepair(sppData, 2, XColumn="Long", YColumn="Lat",
#'                               sppColumn="species", siteColumn="site", predData=envTab)
#'
#' ##plot GDM uncertainty using one core
#' #not run
#' #plotUncertainty(sitePairTab, sampleSites=0.70, bsIters=5, geo=TRUE, plot.layout=c(3,3))
#'
#' ##plot GDM uncertainty in parallel
#' #not run
#' #plotUncertainty(sitePairTab, sampleSites=0.70, bsIters=50, geo=TRUE, plot.layout=c(3,3),
#'                  #parallel=T, cores=10)
#'
#' @keywords gdm
#'
#' @importFrom utils write.csv
#' @importFrom graphics par
#' @importFrom graphics points
#' @importFrom graphics polygon
#' @importFrom graphics rug
#' @importFrom stats sd
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @importFrom doParallel registerDoParallel
#'
#' @export
plotUncertainty <- function(spTable, sampleSites, bsIters, geo=FALSE, splines=NULL,
                            knots=NULL, splineCol="blue", errCol="grey80",
                            plot.linewidth=2.0, plot.layout=c(2,2), parallel=FALSE,
                            cores=2, save=FALSE, fileName="gdm.plotUncertainy.csv"){
  #################
  #spTable <- sitePairTab          ##the input site-pair table to subsample from
  #sampleSites <- 0.9     ##fraction of sites that should be retained from site pair table
  #bsIters <- 50       ##the number of times the site-pair table should be sampled
  #geo <- F              ##rather or not the gdm model takes geography into account, see gdm
  #splines <- NULL       ##splines gdm setting, see gdm
  #knots <- NULL         ##knots gdm setting, see gdm
  #splineCol <- "blue"    ##color of the center line
  #errCol <- "grey80"        ##color of the uncertainty polygon
  #plot.linewidth <- 2.0    ##line width of the center line
  #plot.layout <- c(3,3)    ##number of plots per page
  #parallel <- F            ##rather or not the sampling should happen in parallel processing, to speed it up
  #cores <- 6               ##number of cores to if parallel processing
  #################
  ##function breaks and warnings
  ##makes sure that table is a properly formatted site-pair table
  if(!is(spTable, "gdmData")){
    warning("The spTable object is not of class 'gdmData'. See the formatsitepair function for help.")
  }
  ##checks to makes sure table is a matrix or table frame
  if(!(is(spTable, "gdmData") | is(spTable, "matrix") | is(spTable, "data.frame"))){
    stop("The spTable object must be of class 'gdmData', 'matrix', or 'data.frame'.")
  }

  ##sanity check on the data table
  if(ncol(spTable) < 6){
    stop("spTable object requires at least 6 columns: Observed, weights, s1.xCoord, s1.yCoord, s2.xCoord, s2.yCoord")
  }
  if(nrow(spTable) < 1){
    stop("spTable object has < 1 rows.")
  }

  ##checks that geo has either TRUE or FALSE
  if(!(geo==TRUE | geo==FALSE)){
    stop("geo argument must be either TRUE or FALSE")
  }
  ##makes sure splines is a numeric vector
  if(is.null(splines)==FALSE & !is(splines, "numeric")){
    stop("splines object must of of class = 'numeric'.")
  }
  ##checks knots inputs
  if(is.null(knots)==FALSE & !is(knots, "numeric")){
    stop("knots object must of of class = 'numeric'.")
  }

  ##checks that parallel has either TRUE or FALSE
  if(!(parallel==TRUE | parallel==FALSE)){
    stop("parallel argument must be either TRUE or FALSE")
  }
  ##makes sure that cores has a value when parallel is true
  if(parallel==TRUE & is.null(cores)==TRUE){
    stop("If parallel==TRUE, the number of cores must be specified")
  }
  ##makes sure that cores is a positive integer
  if((is.null(cores)==FALSE & is.numeric(cores)==FALSE) | cores<1){
    stop("argument cores needs to be a positive integer")
  }

  ##makes sure that bsIters is a positive integer
  if((is.null(bsIters)==FALSE & is.numeric(bsIters)==FALSE) | bsIters<1){
    stop("argument bsIters needs to be a positive integer")
  }
  ##makes sure sampleSites is a number
  if(is.numeric(sampleSites)==FALSE){
    stop("sampleSites must be a number between 0 and 1")
  }
  ##makes sure that sampleSites is between 0 and 1
  if(sampleSites < 0){
    stop("sampleSites must be a number between 0 and 1")
  }
  if(sampleSites > 1){
    stop("sampleSites must be a number between 0 and 1")
  }
  if(sampleSites==0){
    stop("a sampleSites value of 0 will remove all sites from the analysis (bad).")
  }

  if(save & is.null(fileName)){
    stop("Save is TRUE, but no fileName provided.")
  }

  ##check again that these values are integers, seems to truncate if not
  cores <- as.integer(cores)
  bsIters <- as.integer(bsIters)

  ##assign k to prevent issues to cran checking
  k <- NULL

  ##makes copies of the site-pair table in order to randomly subsample each one differently
  lstSP <- lapply(1:bsIters, function(i){spTable})

  ##runs parallel if selected
  if(parallel==TRUE){

    ##sets cores
    cl <- makeCluster(cores, outfile="")
    registerDoParallel(cl)
    ##first removes a number of sites according to input
    subSamps <- foreach(k=1:length(lstSP), .verbose=F, .packages=c("gdm")) %dopar%
      subsample.sitepair(lstSP[[k]], sampleSites=sampleSites)
    ##models the subsamples
    gdmMods <- foreach(k=1:length(subSamps), .verbose=F, .packages=c("gdm")) %dopar%
      #gdmMods <- try(foreach(k=1, .verbose=F, .packages=c("gdm")) %dopar%
      gdm(subSamps[[k]], geo=geo, splines=splines, knots=knots)
    stopCluster(cl)
  }else{
    ##first removes a number of sites according to input
    subSamps <- lapply(lstSP, subsample.sitepair, sampleSites=sampleSites)
    ##models the subsamples
    gdmMods <- lapply(subSamps, gdm, geo=geo, splines=splines, knots=knots)
  }

  ##models the full gdm
  fullGDMmodel <- gdm(spTable, geo=geo, splines=splines, knots=knots)

  #get deviance explained by each model
  devExps <- lapply(gdmMods, function(x){
    x$explained
  })
  devExps <- unlist(devExps)

  ##Extracts the splines for each model
  exUncertSplines <- lapply(gdmMods, isplineExtract)
  fullGDMsplines <- isplineExtract(fullGDMmodel)

  ##get the names of the predictor variables
  predVars <- colnames(exUncertSplines[[1]][[1]])

  ##establish what plot layout to use
  thisplot <- 0
  one_page_per_plot <- FALSE
  if ((plot.layout[1]==1) && (plot.layout[2]==1)){
    one_page_per_plot <- TRUE
  }else{
    par(mfrow=plot.layout)
  }

  ##sets the plotting minimum and maximum y-values
  totalYmin <- Inf
  totalYmax <- -Inf

  ##determines the bounds of the plots
  for(p in 1:length(predVars)){
    predV <- predVars[p]

    ##gets the minimum and maximum, to set the ploting extent
    for(nm in 1:length(exUncertSplines)){
      #nm=1
      selPlot <- exUncertSplines[[nm]]

      spYmax <- max(selPlot[[2]][,predV])
      spYmin <- min(selPlot[[2]][,predV])

      totalYmax <- max(c(totalYmax, spYmax))
      totalYmin <- min(c(totalYmin, spYmin))
    }
  }

  ##plots by variable
  outData <- data.frame(matrix(nrow=200, ncol=length(predVars)*6))
  for(p in 1:length(predVars)){
    predV <- predVars[p]

    ##sets the plotting minimum and maximum x-values
    totalXmin <- Inf
    totalXmax <- -Inf
    ##gets the minimum and maximum, to set the ploting extent
    for(nm in 1:length(exUncertSplines)){
      #nm=1
      selPlot <- exUncertSplines[[nm]]

      spXmax <- max(selPlot[[1]][,predV])
      spXmin <- min(selPlot[[1]][,predV])

      if(spXmax > totalXmax){totalXmax = spXmax}
      if(spXmin < totalXmin){totalXmin = spXmin}
    }

    ##checks to make sure that there is some significance, if there is the the data is plotted
    if(totalYmax!=0){
      ##add mean of subsets to plot
      plotX <- NULL
      plotY <- NULL
      byVarMatX <- NULL
      byVarMatY <- NULL
      ##create matrices based on the variable and its x and y location for each model iteration
      for(nn in 1:length(exUncertSplines)){
        plotX[[nn]] <- exUncertSplines[[nn]][[1]]
        plotY[[nn]] <- exUncertSplines[[nn]][[2]]
        byVarMatY <- cbind(byVarMatY, plotY[[nn]][,predV])
        byVarMatX <- cbind(byVarMatX, plotX[[nn]][,predV])
      }
      ##gets spline data from the full gdm model
      fullPlotX <- fullGDMsplines[[1]]
      fullPlotX <- fullPlotX[,predV]
      fullPlotY <- fullGDMsplines[[2]]
      fullPlotY <- fullPlotY[,predV]

      ##calculates the confidence intervals for plotting
      sdX <- apply(as.matrix(byVarMatX), 1, sd)
      sdY <- apply(as.matrix(byVarMatY), 1, sd)
      highBoundX <- fullPlotX + sdX
      lowBoundX <- fullPlotX - sdX
      highBoundY <- fullPlotY + sdY
      lowBoundY <- fullPlotY - sdY

      if(p==1){
        start <- p
        end <- p*6
      } else {
        start <- end+1
        end <- p*6
      }
      outData[,start:end] <- cbind(lowBoundX,
                                   fullPlotX,
                                   highBoundX,
                                   lowBoundY,
                                   fullPlotY,
                                   highBoundY)
      colnames(outData)[start:end] <- paste(predV,
                                          c("minusSD_X",
                                            "fullModel_X",
                                            "plusSD_X",
                                            "minusSD_Y",
                                            "fullModel_Y",
                                            "plusSD_Y"), sep="_")

      ##collects the data to be used in the rug plot
      if(predV=="Geographic"){
        ##calculates unique euclidean distance between sites
        rugData <- unique(sqrt(((spTable$s1.xCoord-spTable$s2.xCoord)^2)+((spTable$s1.yCoord-spTable$s2.yCoord)^2)))
      }else{
        ##gets unique values of variable data
        varDat <- grep(predV, colnames(spTable))
        rugData <- unique(c(spTable[,c(varDat[1])], spTable[,c(varDat[2])]))
      }

      ##plots one graph per page, unless specified otherwise
      if (one_page_per_plot){
        dev.new()
        dev.next()
      }else{
        thisplot <- thisplot + 1
        if(thisplot > (plot.layout[1] * plot.layout[2])){
          thisplot <- 1
          par(mfrow=plot.layout)
        }
      }
      #settings <- par(pars)
      ##plots mean data line and polygon of uncertainty
      plot(NULL, xlim=c(totalXmin, totalXmax), ylim=c(totalYmin, totalYmax), xlab=predV, ylab="Partial Ecological Distance")
      polygon(c(lowBoundX, rev(highBoundX)), c(lowBoundY, rev(highBoundY)), col=errCol, border=NA)
      lines(fullPlotX, fullPlotY, col=splineCol, lwd=plot.linewidth)
      rug(rugData)
    }
  }

  if(save){
    write.csv(outData, fileName, row.names = F)
  }
}

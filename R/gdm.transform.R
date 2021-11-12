#' @title Transform Environmental Data Using a Generalized Dissimilarity Model
#'
#' @description This function transforms geographic and environmental predictors using (1) the
#' fitted functions from a model object returned from \code{\link[gdm]{gdm}} and (2) a
#' data frame or raster stack containing predictor data for a set of sites.
#'
#' @usage gdm.transform(model, data)
#'
#' @param model A gdm model object resulting from a call to \code{\link[gdm]{gdm}}.
#'
#' @param data Either (i) a data frame containing values for each predictor variable in the model, formatted as follows: X, Y, var1, var2, var3, ..., varN or (ii) a raster stack with one layer per predictor variable used in the model, excluding X and Y (rasters for x- and y-coordinates are built automatically from the input rasters if the model was fit with geo=T). The order of the columns (data frame) or raster layers (raster stack) MUST be the same as the order of the predictors in the site-pair table used in model fitting. There is currently no checking to ensure that the order of the variables to be transformed are the same as those in the site-pair table used in model fitting. If geographic distance was not used as a predictor in model fitting, the x- and y-columns need to be removed from the data to be transformed. Output is provided in the same format as the input data.
#'
#' @return
#' gdm.transform returns either a data frame with the same number of rows as the input data frame or a raster stack, depending on the format of the input data. If the model uses geographic distance as a predictor the output object will contain columns or layers for the transformed X and Y values for each site. The transformed environmental data will be in the remaining columns or layers.
#'
#' @references Ferrier S, Manion G, Elith J, Richardson, K (2007) Using generalized dissimilarity modelling to analyse and predict patterns of beta diversity in regional biodiversity assessment. \emph{Diversity & Distributions} 13, 252-264.
#'
#' Fitzpatrick MC, Keller SR (2015) Ecological genomics meets community-level modeling of biodiversity: Mapping the genomic landscape of current and future environmental adaptation. \emph{Ecology Letters} 18: 1-16
#'
#' @examples
#' load(system.file("./data/southwest.RData", package="gdm"))
#' # grab the columns with xy, site ID, and species data
#' sppTab <- southwest[, c("species", "site", "Lat", "Long")]
#'
#' ##fit gdm using rasters
#' rastFile <- system.file("./extdata/stackedVars.grd", package="gdm")
#' envRast <- raster::stack(rastFile)
#' sitePairRast <- formatsitepair(sppTab, 2, XColumn="Long", YColumn="Lat", sppColumn="species",
#'                                siteColumn="site", predData=envRast)
#' ##remove NA values
#' sitePairRast <- na.omit(sitePairRast)
#'
#' ##fit raster GDM
#' gdmRastMod <- gdm(sitePairRast, geo=TRUE)
#'
#' ##raster input, raster output
#' transRasts <- gdm.transform(gdmRastMod, envRast)
#'
#' # map biological patterns
#' rastDat <- sampleRandom(transRasts, 10000)
#' pcaSamp <- prcomp(rastDat)
#'
#' # note the use of the 'index' argument
#' pcaRast <- predict(transRasts, pcaSamp, index=1:3)
#'
#' # scale rasters
#' pcaRast[[1]] <- (pcaRast[[1]]-pcaRast[[1]]@data@min) /
#'   (pcaRast[[1]]@data@max-pcaRast[[1]]@data@min)*255
#' pcaRast[[2]] <- (pcaRast[[2]]-pcaRast[[2]]@data@min) /
#'   (pcaRast[[2]]@data@max-pcaRast[[2]]@data@min)*255
#' pcaRast[[3]] <- (pcaRast[[3]]-pcaRast[[3]]@data@min) /
#'   (pcaRast[[3]]@data@max-pcaRast[[3]]@data@min)*255
#'
#' plotRGB(pcaRast, r=1, g=2, b=3)
#'
#' @keywords gdm
#'
#' @export
gdm.transform <- function(model, data){
  #################
  ##lines used to quickly test function
  #model <- gdmModel
  #data <- climCurrExt
  #data <- cropRasts[[3:nlayers(cropRasts)]]
  #################
  options(warn.FPU = FALSE)
  rastDat <- NULL
  dataCheck <- class(data)

  ##error checking of inputs
  ##checks to make sure a gdm model is given
  if(!is(model, "gdm")){
    stop("model object must be of class 'gdm'.")
  }
  ##checks to make sure data is a correct format
  if(!(dataCheck=="RasterStack" | dataCheck=="RasterLayer" | dataCheck=="RasterBrick" | dataCheck=="data.frame")){
    stop("Data to be transformed must either be a raster object or data frame.")
  }

  ##checks rather geo was T or F in the model object
  geo <- model$geo

  ##turns raster data into dataframe
  if(dataCheck=="RasterStack" | dataCheck=="RasterLayer" | dataCheck=="RasterBrick"){
    ##converts the raster object into a dataframe, for the gdm transformation
    rastDat <- data
    data <- raster::rasterToPoints(rastDat)
    ##determines the cell number of the xy coordinates
    rastCells <- raster::cellFromXY(rastDat, xy=data[,1:2])

    ##checks for NA in the
    checkNAs <- as.data.frame(which(is.na(data), arr.ind=T))
    if(nrow(checkNAs)>0){
      warning("Extracted data from rasters contained NAs. These were automatically removed from the data object to be transformed.")
      data <- na.omit(data)
      rastCells <- rastCells[-c(checkNAs$row)]
    }

    ##if geo was not T in the model, removes the coordinates from the data frame
    if(geo==FALSE){
      data <- data[,3:ncol(data)]
    }
  }

  sizeVal <- 10000000
  ##sets up the data to be transformed into pieces to be transformed
  holdData <- data
  fullTrans <- matrix(0,nrow(holdData),ncol(holdData))
  rows <- nrow(holdData)
  istart <- 1
  iend <- min(sizeVal,rows)
  ##to prevent errors in the transformation of the x and y values when geo is a predictor,
  ##extracts the rows with the minimum and maximum x and y values, these rows will be added
  ##onto the "chuck" given to transform, and then immediately removed after the transformation,
  ##this makes sure that the c++ code will always have access to the minimum and maximum
  ##x and y values
  if(geo==TRUE){
    if(dataCheck=="RasterStack" | dataCheck=="RasterLayer" | dataCheck=="RasterBrick"){
      xMaxRow <- holdData[which.max(holdData[,"x"]),]
      xMinRow <- holdData[which.min(holdData[,"x"]),]
      yMaxRow <- holdData[which.max(holdData[,"y"]),]
      yMinRow <- holdData[which.min(holdData[,"y"]),]
    }
  }

  ##transform the data based on the gdm
  ##part of a loop to prevent memory errors
  while(istart < rows){
    ##Call the dll function
    data <- holdData[istart:iend,]
    ##adds coordinate rows to data to be transformed
    if((dataCheck=="RasterStack" | dataCheck=="RasterLayer" | dataCheck=="RasterBrick") & geo==TRUE){
      data <- rbind(xMaxRow, xMinRow, yMaxRow, yMinRow, data)
    }
    transformed <- matrix(0,nrow(data),ncol(data))
    z <- .C( "GDM_TransformFromTable",
             as.integer(nrow(data)),
             as.integer(ncol(data)),
             as.integer(model$geo),
             as.integer(length(model$predictors)),
             as.integer(model$splines),
             as.double(model$knots),
             as.double(model$coefficients),
             as.matrix(data),
             trandata = as.double(transformed),
             PACKAGE = "gdm")

    ## Convert transformed from a vector into a dataframe before returning...
    nRows <- nrow(data)
    nCols <- ncol(data)

    ## z$trandata is the transformed data vector created
    myVec <- z$trandata
    pos <- 1
    ##fills out dataframe with transformed values
    for (i in seq(from = 1, to = nCols, by = 1)) {
      tmp <- myVec[seq(from=pos, to=pos+nRows-1)]
      transformed[,i] <- tmp
      pos <- pos + nRows
    }

    ##remove the coordinate rows before doing anything else
    if((dataCheck=="RasterStack" | dataCheck=="RasterLayer" | dataCheck=="RasterBrick") & geo==TRUE){
      transformed <- transformed[-c(1:4),]
    }

    ##places the transformed values into the readied data frame
    fullTrans[istart:iend,] <- transformed
    istart <- iend + 1
    iend <- min(istart + (sizeVal-1), rows)
  }

  ##if wanted output data as raster, provides maps raster, or output table
  if(dataCheck=="RasterStack" | dataCheck=="RasterLayer" | dataCheck=="RasterBrick"){
    ##maps the transformed data back to the input rasters
    rastLay <- rastDat[[1]]
    rastLay[] <- NA
    outputRasts <- raster::stack()
    for(nn in 1:ncol(fullTrans)){
      #print(nn)
      #nn=1
      holdLay <- rastLay
      holdLay[rastCells] <- fullTrans[,nn]
      #holdLay[rastCells] <- holdData[,nn]

      outputRasts <- raster::stack(outputRasts, holdLay)
    }
    ##renames raster layers to be the same as the input
    if(geo){
      names(outputRasts) <- c("xCoord", "yCoord", names(rastDat))
    } else {
      names(outputRasts) <- names(rastDat)
    }

    ##get the predictors with non-zero sum of coefficients
    splineindex <- 1
    predInd <- NULL
    for(i in 1:length(model$predictors)){
      #i <- 1
      ##only if the sum of the coefficients associated with this predictor is > 0.....
      numsplines <- model$splines[i]
      if(sum(model$coefficients[splineindex:(splineindex+numsplines-1)])>0){
        predInd <- c(predInd, i)
      }
      splineindex <- splineindex + numsplines
    }
    if(geo){
      predInd <- c(1,2,predInd[-1]+1)
    }

    outputRasts <- outputRasts[[predInd]]

    ##returns rasters
    return(outputRasts)
  }else{
    if(is.null(rastDat)){
      ##if not raster data, sends back the transformed data
      colnames(fullTrans) <- colnames(data)
      return(fullTrans)
    }else{
      ##returns only the transformed variable data as a table, and the cells with which to map to
      colnames(fullTrans) <- colnames(data)
      return(list(fullTrans, rastCells))
    }
  }
}

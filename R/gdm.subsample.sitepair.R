#' @title Remove Sites at Random from a Site-Pair Table
#'
#' @description Randomly selects a number of sites from a given site-pair table
#' and removes them from the site-pair table. It will remove all instances of
#' the sites randomly selected to be removed in both s1 and s2 positions.
#'
#' @usage subsample.sitepair(spTable, sampleSites)
#'
#' @param spTable A site-pair table, same as used to fit a \code{\link[gdm]{gdm}}.
#'
#' @param sampleSites The fraction (0-1, though a value of 0 would be silly,
#' wouldn't it?) of \emph{sites to retain} from the full site-pair table. If
#' less than 1, this argument will completely remove a fraction of sites such
#' that they are not used in the permutation routines.
#'
#' @return A site-pair table, such as one created by \code{\link[gdm]{formatsitepair}},
#' ideally smaller than the one given. In the very rare case where the function
#' determines not to remove any sites, or should the sampleSites argument be 1,
#' then the function will return the full site-pair table.
#'
#' @note This function removes sites, not just site-pairs (rows) from the
#' site-pair table. This function is called from several of the other functions
#' within the gdm package, including the \code{\link[gdm]{plotUncertainty}} and
#' \code{\link[gdm]{gdm.varImp}} functions, for the purposes of subsampling the sites
#' in the site-pair table.
#'
#' @seealso \code{\link[gdm]{formatsitepair}}
#'
#' @examples
#' ##set up site-pair table using the southwest data set
#' sppData <- southwest[c(1,2,13,14)]
#' envTab <- southwest[c(2:ncol(southwest))]
#' sitePairTab <- formatsitepair(sppData, 2, XColumn="Long", YColumn="Lat", sppColumn="species",
#'                               siteColumn="site", predData=envTab)
#'
#' subsample.sitepair(sitePairTab, sampleSites=0.7)
#'
#' @keywords gdm
#'
#' @export
subsample.sitepair <- function(spTable, sampleSites){
  #################
  #spTable <- gdmTab             ##sitepair table
  #sampleSites <- 0.8         ##fraction of sites to keep in sitepair table
  #################

  ##adds error checking to gdm function
  ##checks to see if in site-pair format from formatsitepair function
  if(!is(spTable, "gdmData")){
    warning("The spTable object is not of class 'gdmData'. See the formatsitepair function for help.")
  }
  ##checks to makes sure data is a matrix or data frame
  if(!(is(spTable, "gdmData") | is(spTable, "matrix") | is(spTable, "data.frame"))){
    stop("The spTable object needs to be of class 'gdmData', 'matrix', or 'data.frame'.")
  }
  ##makes sure that sampleSites is a number between 0 and 1,
  ##and that it is not equal to 0
  if(is.numeric(sampleSites)==FALSE | sampleSites<0 | sampleSites>1){
    stop("argument sampleSites needs to be a positive number between 0 and 1")
  }
  if(sampleSites==0){
    stop("sampleSites = 0 will remove all sites from the analysis!")
  }

  if(sampleSites>0){
    sortMatX <- sapply(1:nrow(spTable), function(i, spTab){c(spTab[i,3], spTab[i,5])}, spTab=spTable)
    sortMatY <- sapply(1:nrow(spTable), function(i, spTab){c(spTab[i,4], spTab[i,6])}, spTab=spTable)
    sortMatNum <- sapply(1:nrow(spTable), function(i){c(1,2)})
    sortMatRow <- sapply(1:nrow(spTable), function(i){c(i,i)})
    ##adds a column of NA for index to be added to
    fullSortMat <- cbind(as.vector(sortMatX), as.vector(sortMatY), as.vector(sortMatNum), as.vector(sortMatRow), rep(NA, length(sortMatX)))
    ##assigns sites by unique coordinates
    siteByCoords <- as.data.frame(unique(fullSortMat[,1:2]))
    ##number of sites to expect by coordinates
    numSites <- nrow(siteByCoords)
    ##assigns site index based on coordinates
    for(i in 1:numSites){
      fullSortMat[which(fullSortMat[,1]==siteByCoords[i,1] & fullSortMat[,2]==siteByCoords[i,2]),5] <- i
    }

    ##create index table to know where each site is in input site-pair table
    indexTab <- matrix(NA,nrow(spTable),2)
    for(iRow in 1:nrow(fullSortMat)){
      indexTab[fullSortMat[iRow,4],fullSortMat[iRow,3]] <- fullSortMat[iRow,5]
    }

    ##determines the number of sites to remove
    numToRemove <- round(max(indexTab)*(1-sampleSites))
    ##randomly determines the index of sites to remove
    if(numToRemove > 0){
      rmSites <- sample(1:max(indexTab), numToRemove)
      rmIndexCol1 <- which(indexTab[,1] %in% rmSites)
      rmIndexCol2 <- which(indexTab[,2] %in% rmSites)
      ##creates sampled table
      sampTable <- spTable[-c(unique(c(rmIndexCol1, rmIndexCol2))),]
    }else{
      warning("Percentage of sites to be retained is high enough such that no sites were removed.")
      sampTable <- spTable
    }
  }else{
    sampTable <- spTable
  }

  ##returns the sampled table
  return(sampTable)
}

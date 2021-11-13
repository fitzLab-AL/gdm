#' @title Create Site-Pair Table
#'
#' @description Creates a site-pair table from the lower half of a site-by-site distance
#' (dissimilarity) matrix. This function is called from the
#'  \code{\link[gdm]{formatsitepair}} function and not needed by the user.
#'
#' @usage createsitepair(dist, spdata, envInfo, dXCol, dYCol, siteCol, weightsType,
#' custWeights)
#'
#' @param dist The lower half of a site-by-site distance (dissimilarity) matrix,
#' provided by the \code{\link[gdm]{formatsitepair}} function.
#'
#' @param spdata Input species data, the same as the bioData input to the
#' \code{\link[gdm]{formatsitepair}} function.
#'
#' @param envInfo Input environmental data. Only accepts data tables as input.
#' If the environmental data for \code{\link[gdm]{formatsitepair}} are rasters, the
#' data would have been extracted into table format within
#'  \code{\link[gdm]{formatsitepair}}.
#'
#' @param dXCol Input x coordinate, the same as the XColumn input to the
#' \code{\link[gdm]{formatsitepair}} function.
#'
#' @param dYCol Input y coordinate, the same as the YColumn input to the
#' \code{\link[gdm]{formatsitepair}} function.
#'
#' @param siteCol Site column, taken from either the species or environmental
#' tables.
#'
#' @param weightsType The method of determining the site-pair weights used in
#' model fitting.
#'
#' @param custWeights Custom weights, as a vector, if given by the user.
#'
#' @return A site-pair table with appropriate distance (dissimilarity) and
#' weight columns used for fitting GDM.
#'
#' @note This function is called from the \code{\link[gdm]{formatsitepair}} function
#' and not needed by the user.
#'
#' @seealso \code{\link[gdm]{formatsitepair}}
#'
#' @keywords gdm internal
#'
#' @export
createsitepair <- function(dist, spdata, envInfo, dXCol, dYCol, siteCol,
                           weightsType, custWeights){
  ###########################
  ##lines used to quickly test function
  #dist = distData
  #spdata = bioData
  #envInfo = predData
  #dXCol = XColumn
  #dYCol = YColumn
  #siteCol = siteColumn
  #weightsType = weightType
  #custWeights = custWeights
  ###########################

  weightsType <- as.character(weightsType)
  distance <- as.vector(dist)
  ##calculates total richness = the sum of the two most diverse sites
  if(weightsType[1]=="richness"){
    sppOnly <- spdata[, -c(1,2,3)]
    sppSums <- rowSums(sppOnly)
    sppSiteSums <- cbind(spdata[1], sppSums)
    orderedSums <- sppSiteSums[order(-sppSiteSums[,2]),]
    richTotal <- orderedSums[1,2]+orderedSums[2,2]
  }

  ##Builds index needed for site-pair table format
  s1.xCoord <- s1.yCoord <- s2.xCoord <- s2.yCoord <- NULL
  s1 <- s2 <- NULL

  if((siteCol %in% colnames(envInfo))==T){
    count <- seq(length(unique(envInfo[,siteCol]))-1,1,-1)
  }else{
    count <- seq(length(unique(envInfo[,"siteUltimateCoolness"]))-1,1,-1)
  }
  s1 <- unlist(sapply(seq(length(count),1), function(y){c(s1, rep((max(count)-y)+1, times=y))}))
  s2 <- unlist(sapply(seq(length(count),1), function(y){c(s2, (max(count)-y+2):(max(count)+1))}))

  if(length(s1)!=length(distance)){
    stop("The length of distance values is not the same as the expected number of rows in the site-pair table, unable to proceed.")
  }

  if(weightsType[1]=="equal"){
    print("Site weighting type: Equal")
    weights <- rep(1, times=length(distance))
  }else if(weightsType[1]=="custom"){
    print("Site weighting type: Custom")
    weights <- (custWeights[s1, "weights"] + custWeights[s2, "weights"]) / 2
  }else{
    print("Site weighting type: Richness")
    weights <- (sppSiteSums[s1, "sppSums"] + sppSiteSums[s2, "sppSums"]) / richTotal
  }
  gdmTable <- cbind(distance, weights)

  ##from environmental or species table, copy coordinates for site-pair table
  if((dXCol %in% colnames(envInfo))==T){
    if((siteCol %in% colnames(envInfo))==T){
      checkTab <- table(envInfo[siteCol])
    }else{
      checkTab <- table(envInfo["siteUltimateCoolness"])
    }

    if(sum(checkTab>1)>0){
      stop("The function failed because at least one site has more than one set
           of coordinates associated with it. In other words, sites with the same
           name occur in different locations.")
    }
    s1.xCoord <- envInfo[s1, dXCol]
    s2.xCoord <- envInfo[s2, dXCol]
    s1.yCoord <- envInfo[s1, dYCol]
    s2.yCoord <- envInfo[s2, dYCol]
  }else if((dXCol %in% colnames(spdata))==T){
    s1.xCoord <- spdata[s1, dXCol]
    s2.xCoord <- spdata[s2, dXCol]
    s1.yCoord <- spdata[s1, dYCol]
    s2.yCoord <- spdata[s2, dYCol]
  }else{
    stop("The function failed because multiple sites have the same coordinates. In other
    words, sites with different names occur in the same location.")
  }

  ##sets up output table
  gdmForm <- cbind(gdmTable, s1.xCoord, s1.yCoord, s2.xCoord, s2.yCoord)
  xhold <- which(names(envInfo)==dXCol)
  yhold <- which(names(envInfo)==dYCol)
  sitehold <- which(names(envInfo)==siteCol)
  sitehold2 <- which(names(envInfo)=="siteUltimateCoolness")
  envInfo <- envInfo[-c(xhold, yhold, sitehold, sitehold2)]

  ##fills output table
  if(ncol(envInfo)>0){
    gdmTableFill <- cbind(gdmForm, envInfo[s1,1:ncol(envInfo)], envInfo[s2,1:ncol(envInfo)])
    names.s1 <- paste("s1.",names(envInfo[1:ncol(envInfo)]), sep="")
    names.s2 <- paste("s2.",names(envInfo[1:ncol(envInfo)]), sep="")
    colnames(gdmTableFill) <- c(colnames(gdmTableFill)[1:6], names.s1, names.s2)
  }else{
    gdmTableFill <- gdmForm
  }

  ##returns results
  return(gdmTableFill)
}

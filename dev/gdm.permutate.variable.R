#' Permutate the Values of a Site-Pair Table Variable
#'
#' A function which randomizes the values of the given variables within a site-pair table. This function is called from the \code{\link{gdm.varImp}} function and should not need to be called directly by the user.
#'
#' @usage permutateVarSitePair(spTab, siteVarTab, indexTab, vName)
#'
#' @param spTab A given site-pair table.
#'
#' @param siteVarTab A given site x variable table.
#'
#' @param indexTab A table of index values for the given site-pair table.
#'
#' @param vName Vector of variable names in both the site-pair and site x variable tables.
#'
#' @return A site-pair table, with variable values randomized from the original site-pair table.
#'
#' @note This function is called from the \code{\link{gdm.varImp}} function and the user should not need to access directly.
#'
#' @seealso \code{\link[gdm]{gdm.varImp}}
#'
#' @keywords gdm internal
#'
#' @export
permutateVarSitePair <- function(spTab, siteVarTab, indexTab, vName){
  ##only randomizes the values for a particular variable
  #################
  #spTab <- currSitePair    ##site-pair table
  #siteVarTab <- siteData   ##siteXvar table
  #indexTab <- indexTab     ##table of the index of sites
  #vName <- varChar         ##variables names
  #################
  ##randomizes the row order of the given siteXvar table
  randVarTab <- siteVarTab[sample(nrow(siteVarTab), nrow(siteVarTab)), ]

  #site1x <- siteVarTab$xCoord[1]
  #site1y <- siteVarTab$yCoord[1]
  #checkingIn <- siteVarTab[siteVarTab$xCoord==site1x & siteVarTab$yCoord==site1y,]
  #checkX <- siteVarTab[siteVarTab$xCoord==site1x,]
  #checkingRand <- randVarTab[randVarTab$xCoord==site1x & randVarTab$yCoord==site1y,]

  ##identifies variable columns in randVarTab
  randCols <- grep(paste("^", vName, "$", sep=""), colnames(randVarTab))
  ##identifies variable columns in site-pair table
  spCols1 <- grep(paste("^s1.", vName, "$", sep=""), colnames(spTab))
  spCols2 <- grep(paste("^s2.", vName, "$", sep=""), colnames(spTab))

  ##extracts values based on new index position
  s1var <- sapply(1:nrow(spTab), function(i){randVarTab[indexTab[i,1],randCols]})
  s2var <- sapply(1:nrow(spTab), function(i){randVarTab[indexTab[i,2],randCols]})
  ##places values back into site-pair table
  spTab[,spCols1] <- s1var
  spTab[,spCols2] <- s2var

  #getCoords1 <- spTab[spTab$s1.xCoord==site1x,]
  #getCoords2 <- spTab[spTab$s2.xCoord==site1x,]

  return(spTab)
}

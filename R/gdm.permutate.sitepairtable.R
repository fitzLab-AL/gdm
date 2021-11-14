#' @title Permutate Site-Pair Table Rows, Internal Function
#'
#' @description A function which randomizes the rows of a given site-pair table. This function
#' is called from the \code{\link[gdm]{gdm.varImp}} function and not needed by the user.
#' be called directly by the user.
#'
#' @usage permutateSitePair(spTab, siteVarTab, indexTab, vNames)
#'
#' @param spTab A site-pair table.
#'
#' @param siteVarTab A site x variable table.
#'
#' @param indexTab A table of index values for the site-pair table.
#'
#' @param vNames Vector of variable names in both the site-pair and
#' site x variable tables.
#'
#' @return A new site-pair table with rows randomized.
#'
#' @note This function is called from the \code{\link[gdm]{gdm.varImp}} function
#' and not needed by the user.
#'
#' @seealso \code{\link[gdm]{gdm.varImp}}
#'
#' @keywords gdm internal
#'
permutateSitePair <- function(spTab, siteVarTab, indexTab, vNames){

  #################
  #spTab <- currSitePair    ##site-pair table
  #siteVarTab <- siteData   ##siteXvar table
  #indexTab <- indexTab     ##table of the index of sites
  #vNames <- varNames       ##variables names
  #vNames <- c("awcA", "phTotal", "shcA", "solumDepth", "bio5", "bio19")
  #################

  ##randomizes the row order of the given siteXvar table
  randVarTab <- siteVarTab[sample(nrow(siteVarTab), nrow(siteVarTab)), ]

  #site1x <- siteVarTab$xCoord[1]
  #site1y <- siteVarTab$yCoord[1]
  #checkingIn <- siteVarTab[siteVarTab$xCoord==site1x & siteVarTab$yCoord==site1y,]
  #checkX <- siteVarTab[siteVarTab$xCoord==site1x,]
  #checkingRand <- randVarTab[randVarTab$xCoord==site1x & randVarTab$yCoord==site1y,]

  ##sets up the coordinate values for the randomized site-pair table
  s1xCoord <- sapply(1:nrow(spTab), function(i){randVarTab[indexTab[i,1],1]})
  s1yCoord <- sapply(1:nrow(spTab), function(i){randVarTab[indexTab[i,1],2]})
  s2xCoord <- sapply(1:nrow(spTab), function(i){randVarTab[indexTab[i,2],1]})
  s2yCoord <- sapply(1:nrow(spTab), function(i){randVarTab[indexTab[i,2],2]})

  #print(vNames)
  ##extracts values of other variables
  varLists <- lapply(vNames, function(vn, rvTab, spt, inT){if(vn!="Geographic"){
    ###################
    #vn <- vNames[[2]]
    #rvTab=randVarTab
    #spt=spTab
    #inT=indexTab
    ###################
    ##identifies variable columns in randVarTab
    randCols <- grep(paste("^", vn, "$", sep=""), colnames(rvTab))
    #print(randCols)
    ##identifies variable columns in site-pair table
    spCols <- grep(vn, colnames(spt))

    s1var <- sapply(1:nrow(spt), function(i){rvTab[inT[i,1],randCols]})
    s2var <- sapply(1:nrow(spt), function(i){rvTab[inT[i,2],randCols]})

    return(list(s1var, s2var))
  }
  }, rvTab=randVarTab, spt=spTab, inT=indexTab)
  # unravels the varList into a data.frame of the variable portion of a site-pair table
  bySite <- lapply(1:2, function(i,vlist){sapply(vlist, function(vl,k){vl[[k]]},k=i)}, vlist=varLists)

  if(is(bySite[[1]], "list")){
    site1Vars <- do.call("cbind", bySite[[1]])
    site2Vars <- do.call("cbind", bySite[[2]])
  }else{
    site1Vars <- bySite[[1]]
    site2Vars <- bySite[[2]]
  }
  ##sets up new site-pair table
  newSP <- as.data.frame(cbind(spTab$distance, spTab$weights, s1xCoord, s1yCoord, s2xCoord, s2yCoord, site1Vars, site2Vars))
  colnames(newSP) <- colnames(spTab)
  class(newSP) <- c(class(spTab))

  #getCoords1 <- newSP[newSP$s1.xCoord==site1x,]
  #getCoords2 <- newSP[newSP$s2.xCoord==site1x,]

  return(newSP)
}

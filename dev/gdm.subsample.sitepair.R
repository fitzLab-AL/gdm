########################################################################## 
removeSitesFromSitePair <- function(spTable, sampleSites){
  ##a function to remove a random number of sites from a sitepair table
  ##involves assigning an index to each site, picking the indicies to be 
  ##removed, then identifying which site pairs those indices are a part of
  ##and remove those site pairs from the table
  #################
  #spTable <- gdmTab             ##sitepair table
  #sampleSites <- 0.8         ##fraction of sites to keep in sitepair table
  #################
  
  ##adds error checking to gdm function
  ##checks to see if in site-pair format from formatsitepair function
  if(class(spTable)[1] != "gdmData"){
    warning("spTable class does not include type 'gdmData'. Make sure your data is in site-pair format or the gdm model will not fit.")
  }
  ##checks to makes sure data is a matrix or data frame
  if(!(class(spTable)[1]=="gdmData" | class(spTable)[1]=="matrix" | class(spTable)[1]=="data.frame")){
    stop("spTable argument needs to be gdmData, a matrix, or a data frame")
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
##########################################################################
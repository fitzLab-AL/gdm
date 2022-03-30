library(gdm)

## tabular data
# start with the southwest data table
head(southwest)
sppData <- southwest[c(1,2,13,14)]
envTab <- southwest[c(2:ncol(southwest))]

#########table type 1
## site-species table without coordinates
testData1a <- reshape2::dcast(sppData, site~species)
##site-species table with coordinates
coords <- unique(sppData[, 2:ncol(sppData)])
testData1b <- merge(testData1a, coords, by="site")
## site-species, table-table
spTable <- formatsitepair(testData1a, 1, siteColumn="site", XColumn="Long",
                             YColumn="Lat", predData=envTab)

geo=T

##number of variables in the site-pair table, adds 1 if geo=TRUE
nVars <- (ncol(spTable)-6)/2
##collects variable names
varNames <- colnames(spTable[c(7:(6+nVars))])

varNames <- sapply(strsplit(varNames, "s1."), "[[", 2)
if(geo==TRUE){
  nVars <- nVars + 1
  varNames <- c("Geographic", varNames)
}

##First create a spTable to determine the index of each site in the site-pair table
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

## And remove the sorting table and supporting objects to free up memory
rm(fullSortMat)
rm(sortMatX)
rm(sortMatY)
rm(sortMatNum)
rm(sortMatRow)
rm(siteByCoords)

##create siteXvar table, to be able to rebuild site-pair table later in function
exBySite <- lapply(1:numSites, function(i, index, tab){
  rowSites <- which(index[,1] %in% i)
  if(length(rowSites)<1){
    rowSites <- which(index[,2] %in% i)
  }
  exSiteData <- tab[rowSites[1],]
  return(exSiteData)
}, index=indexTab, tab=spTable)
##identifies the one site not in the first column of the index table
outSite <- which(!(1:numSites %in% indexTab[,1]))
#print(outSite)
##sets up siteXvar table, uses for loop to make sure have steps correct
for(i in 1:length(exBySite)){
  #i <- 42
  ##grabs row and identify if should take s1 or s2 by rather or not number appears in outsite
  siteRow <- exBySite[[i]]
  if(i %in% outSite){
    ##extracts the data from the site-pair table by site
    siteRow <- siteRow[grep("s2.", colnames(siteRow))]
    colnames(siteRow) <- sapply(strsplit(colnames(siteRow), "s2."), "[[", 2)
  }else{
    ##extracts the data from the site-pair table by site
    siteRow <- siteRow[grep("s1.", colnames(siteRow))]
    colnames(siteRow) <- sapply(strsplit(colnames(siteRow), "s1."), "[[", 2)
  }
  exBySite[[i]] <- siteRow
}

##transforms data from list to table
siteData <- do.call("rbind", exBySite)

perm <- replicate(10, list(permutateSitePair(spTable,
                                             siteData,
                                             indexTab,
                                             varNames)))

lll <- lapply(1:length(perm), function(x, spt, pSpt, varChar){
  ppp <- pSpt[[x]]
  idx <- grep(varChar, colnames(ppp))
  spt[,idx] <- ppp[,idx]
  return(spt)
}, spt=spTable, pSpt=perm, varChar="bio18")


eee <- lapply(lll, gdm)


permVar <- permutateVarSitePair(spTable, siteData, indexTab, varChar)

library(gdm)
library(pbapply)
library(doParallel)
library(foreach)

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
splines <- knots <- NULL
fullModelOnly <- T
nPerm <- 100
parallel <- TRUE
cores <- 12
sampleSites <- 1
sampleSitePairs <- 1
outFile <- NULL


##assign k to prevent issues with cran checking
k <- NULL

##error checking for input objects
##checks to see if in site-pair format from formatsitepair function
if(!is(spTable, "gdmData")){
  warning("The spTable object is not of class 'gdmData'. See the formatsitepair function for help.")
}
##checks to makes sure data is a matrix or data frame
if(!(is(spTable, "gdmData") | is(spTable, "matrix") | is(spTable, "data.frame"))){
  stop("spTable argument needs to be of class 'gdmData', 'matrix', or 'data frame'")
}

##sanity check on the data table
if(ncol(spTable) < 6){
  stop("spTable object requires at least 6 columns: distance, weights, s1.xCoord, s1.yCoord, s2.xCoord, s2.yCoord")
}
if(nrow(spTable) < 1){
  stop("The spTable object contains zero rows of data.")
}

##checks that geo has either TRUE or FALSE
if(!(geo==TRUE | geo==FALSE)){
  stop("The geo argument must be either TRUE or FALSE.")
}
##makes sure splines is a numeric vector
if(is.null(splines)==FALSE & !is(splines, "numeric")){
  stop("The splines argument needs to be a numeric data type.")
}
##checks knots inputs
if(is.null(knots)==FALSE & !is(knots, "numeric")){
  stop("The knots argument needs to be a numeric data type.")
}
##checks that fullModelOnly has either TRUE or FALSE
if(!(fullModelOnly==TRUE | fullModelOnly==FALSE)){
  stop("The fullModelOnly argument must be either TRUE or FALSE.")
}
##makes sure that nPerm is a positive integer
if((is.null(nPerm)==FALSE & is.numeric(nPerm)==FALSE) | nPerm<1){
  stop("The nPerm argument needs to be a positive integer.")
}
##checks that parallel has either TRUE or FALSE
if(!(parallel==TRUE | parallel==FALSE)){
  stop("The parallel argument must be either TRUE or FALSE.")
}
##makes sure that cores has a value when parallel is true
if(parallel==TRUE & is.null(cores)==TRUE){
  stop("If parallel==TRUE, the number of cores must be specified.")
}
##makes sure that cores is a positive integer
if((is.null(cores)==FALSE & is.numeric(cores)==FALSE) | cores<1){
  stop("The cores argument needs to be a positive integer.")
}
##makes sure that both sampleSites and sampleSitePairs are a number between 0 and 1,
##and that neither is equal to 0
if(is.numeric(sampleSites)==FALSE | sampleSites<0 | sampleSites>1){
  stop("The sampleSites argument needs to be a positive number between 0 and 1.")
}
if(is.numeric(sampleSitePairs)==FALSE | sampleSitePairs<0 | sampleSitePairs>1){
  stop("The sampleSitePairs argument needs to be a positive number between 0 and 1.")
}
if(sampleSites==0){
  stop("A sampleSites value of 0 will remove all sites from the analysis.")
}
if(sampleSitePairs==0){
  stop("A sampleSitePairs value of 0 will remove all sites from the analysis.")
}
##checks to see if the user has requested for an output file to be written, and if so
##makes sure that it is formatted correctly
if(is.null(outFile)==FALSE){
  ##first makes sure outFile is a string
  if(is.character(outFile)==FALSE){
    stop("The outFile argument needs to be a character string of the directory and file name you wish the tables to be written to")
  }
  ##makes sure that text has ".RData" in it, if not, adds it
  outFileChar <- nchar(outFile)
  if(substr(outFile, outFileChar-5, outFileChar)!=".RData"){
    outFile <- paste(outFile, ".RData", sep="")
  }
  ##checks to see if there is a path as well as a file name
  if(length(strsplit(outFile,"/")[[1]])>1){
    splitOutFile <- strsplit(outFile,"/")[[1]][-length(strsplit(outFile,"/")[[1]])]
    dir.create(paste(splitOutFile, collapse="/"))
  }else{
    outFile <- paste("./", outFile, sep="")
  }
}

##double makes sure these values are integers, seems to truncate if not
nPerm <- as.integer(nPerm)
cores <- as.integer(cores)

##removes a user specified number of sites from the site-pair table
if(sampleSites<1){
  spTable <- subsample.sitepair(spTable, sampleSites=sampleSites)
  ##throws warning if sampleSitePairs<1 as well
  if(sampleSitePairs<1){
    warning("You have selected to randomly remove sites and/or site-pairs.")
  }
}
##removes a user specified number of site-pairs from the site-pair table
if(sampleSitePairs<1){
  ##determine which rows to remove
  numRm <- sample(1:nrow(spTable), round(nrow(spTable)*(1-sampleSitePairs)))
  spTable <- spTable[-c(numRm),]
}

##check that the response data is [0..1]
rtmp <- spTable[,1]
if(length(rtmp[rtmp<0]) > 0){
  stop("The spTable contains negative distance values. Must be between 0 - 1.")
}
if (length(rtmp[rtmp>1]) > 0){
  stop("The spTable contains distance values greater than 1. Must be between 0 - 1.")
}

# run initial GDM to see if any vars have zero I-spline coeffs
print(paste0("Fitting initial model with all predictors..."))
fullGDM <- gdm(spTable, geo=geo, splines=splines, knots=knots)

# check for zero coeffs
thiscoeff <- 1
thisquant <- 1
sumCoeff <- NULL
for(i in 1:length(fullGDM$predictors)){
  numsplines <- fullGDM$splines[[i]]
  holdCoeff <- NULL
  for(j in 1:numsplines){
    holdCoeff[j] <- fullGDM$coefficients[[thiscoeff]]
    thiscoeff <- thiscoeff + 1
  }
  sumCoeff[i] <- sum(holdCoeff)
}

zeroSum <- fullGDM$predictors[which(sumCoeff==0)]
if(length(zeroSum)>0){
  print(paste0("Sum of I-spline coefficients for predictors(s) ", zeroSum," = 0"))
  print(paste0("Removing predictor(s) ", zeroSum, " and proceeding with permutation testing..."))

  for(z in zeroSum){
    ##select variable columns to be removed from original site-pair table
    testVarCols1 <- grep(paste("^s1.", z, "$", sep=""), colnames(spTable))
    testVarCols2 <- grep(paste("^s2.", z, "$", sep=""), colnames(spTable))
    spTable <- spTable[,-c(testVarCols1, testVarCols2)]
  }
}

# number of variables in the site-pair table, adds 1 if geo=TRUE
nVars <- (ncol(spTable)-6)/2
# create vector of variable names
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

##create site x predictor table, to be able to rebuild site-pair table later in function
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

##sets up objects to be returned by the function
modelTestValues <- data.frame(matrix(NA,4,nVars,dimnames = list(c("Model deviance",
                                                                  "Percent deviance explained",
                                                                  "Model p-value",
                                                                  "Fitted permutations"),
                                                                c("All predictors",
                                                                  "1-removed",
                                                                  paste(seq(2,nVars-1), "-removed", sep="")))))
##deviance reduction predictor table
varImpTable <- matrix(NA, nVars, nVars-1)
rownames(varImpTable) <- varNames
colnames(varImpTable) <- c("All predictors",
                           "1-removed",
                           paste(seq(2,nVars-2), "-removed", sep=""))
##p value predictor table
pValues <- nModsConverge <- varImpTable

##assigns given site-pair table to new variable, to prevent changing the original input
currSitePair <- spTable
nullGDMFullFit <- 0  ##a variable to track rather or not the fully fitted gdm model returned a NULL object

# create the set up permutated site-pair tables to be used for all
# downstream analyses
print(paste0("Creating ", nPerm, " permutated site-pair tables..."))

permSpt <- pbreplicate(nPerm, list(permutateSitePair(currSitePair,
                                                     siteData,
                                                     indexTab,
                                                     varNames)))
varNames.x <- varNames
for(v in 1:length(varNames)){
  # runs gdm, first time on site-pair table with all variables
  # As variables are eliminated the "full" site-pair table will
  # contain fewer variables
  fullGDM <- gdm(currSitePair, geo=geo, splines=splines, knots=knots)
  print(paste0("Percent deviance explained by the full model =  ", round(fullGDM$explained,3)))

  if(is.null(fullGDM)==TRUE){
    warning(paste("The model did not converge when testing variable: ", varNames[v],
                  ". Terminating analysis and returning output completed up to this point.", sep=""))
    break
  }

  print("Fitting GDMs to the permutated site-pair tables...")
  permGDM <- lapply(permSpt, function(x){
    gdm(x, geo=geo, splines=splines, knots=knots)
  })

  ##extracts deviance of permuted gdms
  permModelDev <- sapply(permGDM, function(mod){mod$gdmdeviance})
  ##if needed, removes nulls from output permModelDev
  modPerms <- length(which(sapply(permModelDev,is.null)==TRUE))
  if(modPerms>0){
    permModelDev <- unlist(permModelDev[-(which(sapply(permModelDev,is.null)==T))])
  }

  if(v>1){
    colnames(modelTestValues)[v] <- colnames(pValues)[v] <- colnames(varImpTable)[v] <- colnames(nModsConverge)[v] <- paste0(names(elimVar),"-removed")
  }

  ##begins to fill in the output table with data from fully fitted model
  modelTestValues[1,v] <- round(fullGDM$gdmdeviance,3)
  modelTestValues[2,v] <- round(fullGDM$explained,3)
  #p-value
  modelTestValues[3,v] <- round(sum(permModelDev<=fullGDM$gdmdeviance)/(nPerm-modPerms),3)
  ##fitted permutations
  modelTestValues[4,v] <- round(nPerm-modPerms,0)

  # now permutate each variable individual and fit models
  if(parallel==TRUE){

    # set up parallel processing
    cl <- makeCluster(cores, outfile="")
    registerDoParallel(cl)

    # foreach function to create site-pair tables with each variable permuted,
    # fit gdms and extract deviance.
    permVarDev <- foreach(k=1:length(varNames.x), .verbose=F, .packages=c("gdm"), .export = c("currSitePair")) %dopar%{
      if(varNames.x[k]!="Geographic"){
        # permute a single variable
        lll <- lapply(permSpt, function(x, spt=currSitePair){
          idx <- grep(varNames.x[k], colnames(x))
          spt[,idx] <- x[,idx]
          return(spt)
        })}

      if(varNames.x[k]=="Geographic"){
        # permute a single variable
        lll <- lapply(permSpt, function(x, spt=currSitePair){
          s1 <- sample(1:nrow(spt), nrow(spt))
          s2 <- sample(1:nrow(spt), nrow(spt))
          s3 <- sample(1:nrow(spt), nrow(spt))
          s4 <- sample(1:nrow(spt), nrow(spt))
          spt[,3] <- spt[s1,3]
          spt[,4] <- spt[s2,4]
          spt[,5] <- spt[s3,5]
          spt[,6] <- spt[s4,6]
          return(spt)
        })}

      gdmPermVar <- lapply(lll, function(x){
        try(gdm(x, geo=geo, splines=splines, knots=knots))
      })

      ##extracts deviance of permuted gdms
      permModelDev <- sapply(gdmPermVar, function(mod){mod$gdmdeviance})
    }

    ##closes cores
    stopCluster(cl)
  }

  names(permVarDev) <- varNames.x
  nullDev <- fullGDM$nulldeviance

  for(var in varNames.x){
    grepper <- grep(var, names(permVarDev))
    varDevTab <- permVarDev[[grepper]]

    # number of perms for which GDM converged
    nConv <- length(which(is.null(varDevTab)))
    nModsConverge[which(rownames(varImpTable) == var),v] <- nPerm-nConv

    # remove NULLs (GDMs that did not converge)
    if(nConv>0){
      varDevTab <- unlist(varDevTab[-(which(sapply(is.null(varDevTab))))])
    }

    # calculate variable importance
    varDevExplained <- 100*(nullDev-varDevTab)/nullDev
    varImpTable[which(rownames(varImpTable) == var),v] <- median(100 * abs((varDevExplained - fullGDM$explained)/fullGDM$explained))

    if(var!="Geographic"){
      ##select variable columns to be removed from original site-pair table
      testVarCols1 <- grep(paste("^s1.", var, "$", sep=""), colnames(currSitePair))
      testVarCols2 <- grep(paste("^s2.", var, "$", sep=""), colnames(currSitePair))
      testSitePair <- currSitePair[,-c(testVarCols1, testVarCols2)]
      ##run gdm for the missing variable
      noVarGDM <- gdm(testSitePair, geo=geo, splines=splines, knots=knots)
    } else {
      noVarGDM <- gdm(currSitePair, geo=F, splines=splines, knots=knots)
    }

    permDevReduct <- noVarGDM$gdmdeviance - varDevTab
    pValues[which(rownames(pValues) == var),v] <- sum(permDevReduct>=(varDevTab - fullGDM$gdmdeviance))/(nPerm-nConv)
  }

  elimVar <- which.min(varImpTable[,v])
  print(paste0("Removing ", names(elimVar), " and proceeding with the next round of permutations."))

  if(names(elimVar)!="Geographic"){
    ##select variable columns to be removed from original site-pair table
    remVarCols1 <- grep(paste("^s1.", names(elimVar), "$", sep=""), colnames(currSitePair))
    remVarCols2 <- grep(paste("^s2.", names(elimVar), "$", sep=""), colnames(currSitePair))
    currSitePair <- currSitePair[,-c(remVarCols1, remVarCols2)]
  } else {
    geo <- F
  }

  permSpt <- lapply(permSpt, function(x){
    x[,-c(remVarCols1, remVarCols2)]
  })

  varNames.x <- varNames.x[-which(varNames.x==names(elimVar))]

  # ends the loop if only 1 variable remains
  if(length(varNames)<2){
    break
  }
}

##lists tables into one object
outObject <- list(round(modelTestValues, 3), round(devReductVars,3), round(pValues,3), numPermsFit)
##if given, writes out files to space on disk
if(is.null(outFile)==FALSE){
  save(outObject, file=outFile)
}
return(outObject)



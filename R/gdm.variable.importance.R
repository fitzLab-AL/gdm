#' @title Quantify Model Significance and Variable Importance in a Fitted
#' Generalized Dissimilarity Model Using Matrix Permutation.
#'
#' @description This function uses matrix permutation to perform model and
#' variable significance testing and to estimate variable importance in a
#' generalized dissimilarity model. The function can be run in parallel on
#' multicore machines to reduce computation time.
#'
#' @usage gdm.varImp(spTable, geo, splines = NULL, knots = NULL,
#' fullModelOnly = FALSE, nPerm = 50, parallel = FALSE, cores = 2,
#' sampleSites = 1, sampleSitePairs = 1, outFile = NULL)
#'
#' @param spTable A site-pair table, same as used to fit a \code{\link[gdm]{gdm}}.
#'
#' @param geo Similar to the \code{\link[gdm]{gdm}} geo argument. The only
#' difference is that the geo argument does not have a default in this function.
#'
#' @param splines Same as the \code{\link[gdm]{gdm}} splines argument.
#'
#' @param knots Same as the \code{\link[gdm]{gdm}} knots argument.
#'
#' @param fullModelOnly Set to TRUE to test only the full variable set. Set to
#' FALSE to estimate model significance and variable importance and significance
#' using matrix permutation and backward elimination. Default is FALSE.
#'
#' @param nPerm Number of permutations to use to estimate p-values. Default is 50.
#'
#' @param parallel Whether or not to run the matrix permutations and model
#' fitting in parallel. Parallel processing is highly recommended when either
#' (i) the nPerms argument is large (>100) or (ii) a large number of site-pairs
#' (and / or variables) are being used in model fitting (note computation demand
#' can be reduced using subsampling - see next arguments). The default is FALSE.
#'
#' @param cores When the parallel argument is set to TRUE, the number of cores
#' to be registered for parallel processing. Must be <= the number of cores in
#' the machine running the function.
#'
#' @param sampleSites The fraction (0-1, though a value of 0 would be silly,
#' wouldn't it?) of \emph{sites to retain} from the full site-pair table. If
#' less than 1, this argument will completely remove a fraction of sites such
#' that they are not used in the permutation routines.
#'
#' @param sampleSitePairs The fraction (0-1) of \emph{site-pairs (i.e., rows)
#' to retain} from the full site-pair table - in other words, all sites will
#' be used in the permutation routines (assuming sampleSites = 1), but not
#' all \emph{site-pair combinations}. In the case where both the sampleSites
#' and the sampleSitePairs argument have values less than 1, sites first will
#' be removed using the sampleSites argument, followed by removal of site-pairs
#' using the sampleSitePairs argument. Note that the number of site-pairs
#' removed is based on the fraction of the resulting site-pair table after
#' sites have been removed, not on the size of the full site-pair table.
#'
#' @param outFile An optional character string to write the object returned by
#' the function to disk as an .RData object (".RData" is not required as part
#' of the file name). The .RData object will contain a single list with the
#' name of "outObject". The default is NULL, meaning that no file will be written.
#'
#' @details To test model significance, first a "full model" is fit using
#' un-permuted environmental data. Next, the environmental data are permuted
#' nPerm times (by randomizing the order of the rows) and a GDM is fit to each
#' permuted table. Model significance is determined by comparing the deviance
#' explained by GDM fit to the un-permuted table to the distribution of deviance
#' explained values from GDM fit to the nPerm permuted tables. To assess variable
#' significance, this process is repeated for each predictor individually (i.e.,
#' only the data for the variable being tested is permuted rather than the entire
#' environmental table). Variable importance is quantified as the percent change
#' in deviance explained between a model fit with and without that variable
#' (technically speaking, with the variable permuted and un-permuted). If
#' fullModelOnly=FALSE, this process continues by next permutating the site-pair
#' table nPerm times, but removing one variable at a time and reassessing
#' variable importance and significance. At each step, the least important
#' variable is dropped (backward elimination) and the process continues until
#' all non-significant predictors are removed.
#'
#' @return A list of four tables. The first table summarizes full model deviance,
#' percent deviance explained by the full model, the p-value of the full model,
#' and the number of permutations used to calculate the statistics for each
#' fitted model (i.e., the full model and each model with variables removed in
#' succession during the backward elimination procedure if fullModelOnly=F). The
#' remaining three tables summarize (1) variable importance, (2) variable
#' significance, and (3) the number of permutations used to calculate the
#' statistics for that model, which is provided because some GDMs may fail
#' to converge for some permutations / variable combinations and you might want to
#' know how many permutations were used when calculating statistics. Or maybe you
#' don't, you decide.
#'
#' Variable importance is measured as the percent decrease in deviance explained
#' between the full model and the deviance explained by a model fit with that variable
#' permuted. Significance is estimated using the bootstrapped p-value when the
#' variable has been permuted. For most cases, the number of permutations will
#' equal the nPerm argument. However, the value may be less should any of the
#' permutations fail to converge.
#'
#' If fullModelOnly=T, the tables will have values only in the first column and
#' NAs elsewhere.
#'
#' NOTE: In some cases, GDM may fail to converge if there is a weak relationship
#' between the response and predictors (e.g., when an important variable is
#' removed). Such cases are indicated by -9999 values in the variable importance,
#' variable significance, and number of permutations tables.
#'
#' @author Karel Mokany and Matt Fitzpatrick
#'
#' @references Ferrier S, Manion G, Elith J, Richardson, K (2007) Using
#' generalized dissimilarity modelling to analyse and predict patterns of
#' beta diversity in regional biodiversity assessment. \emph{Diversity &
#' Distributions} 13, 252-264.
#'
#' Fitzpatrick, MC, Sanders NJ, Ferrier S, Longino JT, Weiser MD, and RR Dunn.
#' 2011. Forecasting the Future of Biodiversity: a Test of Single- and
#' Multi-Species Models for Ants in North America. \emph{Ecography} 34: 836-47.
#'
#' @examples
#' ##fit table environmental data
#' ##set up site-pair table using the southwest data set
#' sppData <- southwest[c(1,2,13,14)]
#' envTab <- southwest[c(2:ncol(southwest))]
#' sitePairTab <- formatsitepair(sppData, 2, XColumn="Long", YColumn="Lat",
#' sppColumn="species", siteColumn="site", predData=envTab)
#'
#' ## not run
#' #modTest <- gdm.varImp(sitePairTab, geo=T, nPerm=50, parallel=T, cores=10)
#' #barplot(sort(modTest[[2]][,1], decreasing=T))
#'
#' @keywords gdm
#'
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @importFrom doParallel registerDoParallel
#'
#' @export
gdm.varImp <- function(spTable, geo, splines=NULL, knots=NULL, fullModelOnly=FALSE,
                       nPerm=50, parallel=FALSE, cores=2, sampleSites=1,
                       sampleSitePairs=1, outFile=NULL){
 #################
  #spTable <- sitePairTab         ##the input site-pair table to subsample from
  #load("M:/UAE/kavyaWorking/Code/GDM/GDMSitepairTable.RData")
  #spTable <- ddd
  #spTable <- gdmTab[-c(samp),]
  #geo <- T              ##rather or not the gdm model takes geography into account, see gdm
  #splines <- NULL       ##splines gdm setting, see gdm
  #knots <- NULL         ##knots gdm setting, see gdm
  #fullModelOnly <- T     ##rather to run the full calculations, or just once on the full model, acceptable values are TRUE and FALSE
  #nPerm <- 10             ##number of permutations
  #parallel <- T          ##rather or not to run in parallel
  #cores <- 7             ##if in parallel, the number of cores to run on
  #sampleSites <- 1   ##the percent of sites to be retained before calculating
  #sampleSitePairs <- 1  ##the percent of site-pairs (rows) to be retained before calculating
  #outFile <- NULL
  #outFile <- "D:/junk/testWrite.RData"
  #outFile <- "testWrite3.RData"
  #################
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

  ##sets up objects to be returned by the function
  modelTestValues <- matrix(NA,4,nVars,dimnames = list(c("Model deviance",
                                                         "Percent deviance explained",
                                                         "Model p-value",
                                                         "Fitted permutations"),
                                                       c("All variables",
                                                         "1-removed",
                                                         paste(seq(2,nVars-1), "-removed", sep=""))))
  ##deviance reduction variable table
  devReductVars <- matrix(NA, nVars, nVars-1)
  rownames(devReductVars) <- varNames
  colnames(devReductVars) <- c("All variables",
                               "1-removed",
                               paste(seq(2,nVars-2), "-removed", sep=""))
  ##p value variable table
  pValues <- numPermsFit <- devReductVars

  ##assigns given site-pair table to new variable, to prevent changing the original input
  currSitePair <- spTable
  nullGDMFullFit <- 0  ##a variable to track rather or not the fully fitted gdm model returned a NULL object

  for(v in 1:nVars){
    #v <- 1
    #print(varNames[v])

    ##runs gdm, first time on full site-pair table
    ##however as variables are removed the "full" site-pair table will have less variables in it
    fullGDM <- gdm(currSitePair, geo=geo, splines=splines, knots=knots)

    if(is.null(fullGDM)==TRUE){
      warning(paste("The model did not fit when testing variable: ", varNames[v],
                    ". Terminating analysis and returning output completed up to this point.", sep=""))
      break
    }

    ##create a series of permutated site-pair tables, randomized site comparisons
    if(parallel==TRUE){
      #require(foreach)
      #require(doParallel)

      ##sets cores
      cl <- makeCluster(cores, outfile="")
      registerDoParallel(cl)

      permSitePairs <- foreach(k=1:nPerm, .verbose=F, .packages=c("gdm"), .export=c("permutateSitePair")) %dopar%
        permutateSitePair(currSitePair, siteData, indexTab, varNames)
      ##runs gdm on the permuted tables
      permGDM <- try(foreach(k=1:length(permSitePairs), .verbose=F, .packages=c("gdm")) %dopar%
                       gdm(permSitePairs[[k]], geo=geo, splines=NULL, knots=NULL))
      ##closes cores
      stopCluster(cl)
    }else{
      ##non-parallel version
      permSitePairs <- lapply(1:nPerm, function(i, csp, sd, it, vn){permutateSitePair(csp,sd,it,vn)},
                              csp=currSitePair, sd=siteData, it=indexTab, vn=varNames)
      ##runs gdm on the permuted tables
      permGDM <- lapply(permSitePairs, gdm, geo=geo, splines=NULL, knots=NULL)
    }

    ##extracts deviance of permuted gdms
    permModelDev <- sapply(permGDM, function(mod){mod$gdmdeviance})
    ##if needed, removes nulls from output permModelDev
    modPerms <- length(which(sapply(permModelDev,is.null)==TRUE))
    if(modPerms>0){
      permModelDev <- unlist(permModelDev[-(which(sapply(permModelDev,is.null)==T))])
    }

    ##begins to fill in the output table with data from fully fitted model
    modelTestValues[1,v] <- fullGDM$gdmdeviance
    modelTestValues[2,v] <- fullGDM$explained
    #p-value
    modelTestValues[3,v] <- sum(permModelDev<=fullGDM$gdmdeviance)/(nPerm-modPerms)
    ##fitted permutations
    modelTestValues[4,v] <- nPerm-modPerms

    ##ends the loop if only 1 variable was used in the model
    if(length(varNames)<2){
      break
    }

    ##begins running tests on variations
    ##runs model without geo if geo was part of the model
    if(geo==TRUE){
      #print("doing Geo")
      noGeoGDM <- gdm(currSitePair, geo=FALSE, splines=NULL, knots=NULL)

      ##create a series of permutated site-pair tables, randomized site comparisons
      if(parallel == TRUE){
        ##sets cores
        cl <- makeCluster(cores, outfile="")
        registerDoParallel(cl)

        permSitePairs <- foreach(k=1:nPerm, .verbose=F, .packages=c("gdm"), .export=c("permutateSitePair")) %dopar%
          permutateSitePair(currSitePair, siteData, indexTab, varNames)

        permGDM <- try(foreach(k=1:length(permSitePairs), .verbose=F, .packages=c("gdm")) %dopar%
                         gdm(permSitePairs[[k]], geo=geo, splines=NULL, knots=NULL))
        ##closes cores
        stopCluster(cl)
      }else{
        ##non-parallel version
        permSitePairs <- lapply(1:nPerm, function(i, csp, sd, it, vn){permutateSitePair(csp,sd,it,vn)},
                                csp=currSitePair, sd=siteData, it=indexTab, vn=varNames)
        ##runs gdm on the permuted tables
        permGDM <- lapply(permSitePairs, gdm, geo=geo, splines=NULL, knots=NULL)
      }

      ##runs gdm on the permuted tables
      permModelDev <- sapply(permGDM, function(mod){mod$gdmdeviance})
      ##if needed, removes nulls from output permModelDev
      modPerms <- length(which(sapply(permModelDev,is.null)==TRUE))
      if(modPerms>0){
        permModelDev <- unlist(permModelDev[-(which(sapply(permModelDev,is.null)==T))])
      }

      ##calculates table values for geographic, and adds them to output objects
      if(is.null(noGeoGDM$gdmdeviance)==TRUE){
        permDevReduct <- -9999
        devReductVars[1,v] <- -9999
        pValues[1,v] <- -9999
      }else{
        permDevReduct <- noGeoGDM$gdmdeviance - permModelDev
        ##change in devience with variable removed
        devReductVars[1,v] <- 100 * abs((noGeoGDM$explained - fullGDM$explained)/fullGDM$explained)  ##new - percent change in devience
        pValues[1,v] <- sum(permDevReduct>=(noGeoGDM$gdmdeviance - fullGDM$gdmdeviance))/(nPerm-modPerms)
      }
      numPermsFit[1,v] <- nPerm-modPerms
    }

    ##now tests all other variables
    for(varChar in varNames){
      #varChar <- varNames[2]
      if(varChar!="Geographic"){
        ##select variable columns to be removed from original site-pair table
        testVarCols1 <- grep(paste("^s1.", varChar, "$", sep=""), colnames(currSitePair))
        testVarCols2 <- grep(paste("^s2.", varChar, "$", sep=""), colnames(currSitePair))
        testSitePair <- currSitePair[,-c(testVarCols1, testVarCols2)]
        ##run gdm for the missing variable
        noVarGDM <- gdm(testSitePair, geo=geo, splines=NULL, knots=NULL)

        ##create a series of permutated site-pair tables, randomized site comparisons
        if(parallel == TRUE){
          ##sets cores
          cl <- makeCluster(cores, outfile="")
          registerDoParallel(cl)

          noVarSitePairs <- foreach(k=1:nPerm, .verbose=F, .packages=c("gdm"), .export=c("permutateVarSitePair")) %dopar%
            permutateVarSitePair(currSitePair, siteData, indexTab, varChar)
          ##runs gdm on the permuted tables
          permGDM <- try(foreach(k=1:length(noVarSitePairs), .verbose=F, .packages=c("gdm")) %dopar%
                           gdm(noVarSitePairs[[k]], geo=geo, splines=NULL, knots=NULL))
          ##closes cores
          stopCluster(cl)
        }else{
          ##non-parallel version
          noVarSitePairs <- lapply(1:nPerm, function(i, csp, sd, it, vn){permutateVarSitePair(csp,sd,it,vn)},
                                   csp=currSitePair, sd=siteData, it=indexTab, vn=varChar)
          ##runs gdm on the permuted tables
          permGDM <- lapply(noVarSitePairs, gdm, geo=geo, splines=NULL, knots=NULL)
        }

        ##extracts deviance of permuted gdms
        permModelDev <- sapply(permGDM, function(mod){mod$gdmdeviance})
        ##if needed, removes nulls from output permModelDev
        modPerms <- length(which(sapply(permModelDev,is.null)==TRUE))
        if(modPerms>0){
          permModelDev <- unlist(permModelDev[-(which(sapply(permModelDev,is.null)==T))])
        }

        ##calculates table values for geographic, and adds them to output objects
        if(is.null(noVarGDM$gdmdeviance)==TRUE){
          permDevReduct <- -9999
          ggg <- which(rownames(devReductVars) %in% varChar)
          devReductVars[ggg,v] <- rep(-9999, times=length(ggg))
          pValues[ggg,v] <- rep(-9999, times=length(ggg))
        }else{
          permDevReduct <- noVarGDM$gdmdeviance - permModelDev
          devReductVars[which(rownames(devReductVars) %in% varChar),v] <- 100 * abs((noVarGDM$explained - fullGDM$explained)/fullGDM$explained)
          pValues[which(rownames(pValues) %in% varChar),v] <- sum(permDevReduct>=(noVarGDM$gdmdeviance - fullGDM$gdmdeviance))/(nPerm-modPerms)
        }
        numPermsFit[which(rownames(numPermsFit) %in% varChar),v] <- nPerm-modPerms
      }
    }

    ##if fullModelOnly == TRUE, breaks script
    if(fullModelOnly==TRUE){
      break
    }

    ##based on the P-value, and then deviance reduction, select the variable to be omitted
    ##from future iterations of the testing
    tempPVals <- as.numeric(pValues[c(1:nVars),v])
    tempDevs <- as.numeric(devReductVars[c(1:nVars),v])
    tempPVals <- tempPVals[!is.na(tempPVals)]
    tempDevs <- tempDevs[!is.na(tempDevs)]
    #print(length(tempPVals))
    #print(length(tempDevs))
    varToOmit <- which.max(tempPVals)

    for(iCheck in 1:length(varNames)){
      if(tempPVals[iCheck] == tempPVals[varToOmit]){
        if(tempDevs[iCheck] < tempDevs[varToOmit]){
          varToOmit <- iCheck
        }
      }
    }

    ##removes variables
    ##if selected, removes geo
    if(varToOmit==1 & geo==TRUE){
      geo <- FALSE
      varNames <- varNames[-1]
    }else{
      ##removes any non-geo variable
      nameToRemove <- varNames[varToOmit]
      #nameToRemove <- "bio1"
      ##remove from variables
      varNames <- varNames[-varToOmit]
      removeFromSitePs1 <- grep(paste("^s1.", nameToRemove, "$", sep=""), colnames(currSitePair))
      removeFromSitePs2 <- grep(paste("^s2.", nameToRemove, "$", sep=""), colnames(currSitePair))
      ##removes variable from important objects
      currSitePair <- currSitePair[,-c(removeFromSitePs1,removeFromSitePs2)]
    }
  }
  ##lists tables into one object
  outObject <- list(round(modelTestValues, 3), round(devReductVars,3), round(pValues,3), numPermsFit)
  ##if given, writes out files to space on disk
  if(is.null(outFile)==FALSE){
    save(outObject, file=outFile)
  }
  return(outObject)
}

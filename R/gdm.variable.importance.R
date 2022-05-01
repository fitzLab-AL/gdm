#' @title Assess Predictor Importance and Quantify Model Significance in a Fitted
#' Generalized Dissimilarity Model.
#'
#' @description This function uses matrix permutation to perform model and
#' predictor significance testing and to estimate predictor importance in a
#' generalized dissimilarity model. The function can be run in parallel on
#' multicore machines to reduce computation time.
#'
#' @usage gdm.varImp(spTable, geo, splines = NULL, knots = NULL,
#' predSelect = FALSE, nPerm = 50, pValue=0.05, parallel = FALSE, cores = 2,
#' sampleSites = 1, sampleSitePairs = 1, outFile = NULL)
#'
#' @param spTable A site-pair table, same as used to fit a \code{\link[gdm]{gdm}}.
#'
#' @param geo Similar to the \code{\link[gdm]{gdm}} geo argument. The only
#' difference is that the geo argument does not have a default in this function.
#'
#' @param splines Same as the \code{\link[gdm]{gdm}} splines argument. Note that
#' the current implementation requires that all predictors have the same number of
#' splines.
#'
#' @param knots Same as the \code{\link[gdm]{gdm}} knots argument.
#'
#' @param predSelect Set to TRUE to perform predictor selection using matrix
#' permutation and backward elimination. Default is FALSE. When predSelect = FALSE
#' results will be returned only for a model fit with all predictors.
#'
#' @param nPerm Number of permutations to use to estimate p-values. Default is 50.
#'
#' @param pValue The p-value to use for predictor selection / elimination. Default is 0.05.
#'
#' @param parallel Whether or not to run the matrix permutations and model
#' fitting in parallel. Parallel processing is highly recommended when either
#' (i) the nPerms argument is large (>100) or (ii) a large number of site-pairs
#' (and / or variables) are being used in model fitting (note computation demand
#' can be reduced using subsampling - see next arguments). The default is FALSE.
#'
#' @param cores When the parallel argument is set to TRUE, the number of cores
#' to be registered for parallel processing. Must be <= the number of cores in
#' the machine running the function. There is no benefit to setting the number of
#' cores greater than the number of predictors in the model.
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
#' @details To test model significance, first a model is fit using all predictors and
#' un-permuted environmental data. Any predictor for which the sum of the I-spline
#' coefficients sum to zero is preemptively removed. Next, the environmental data are permuted
#' nPerm times (by randomizing the order of the rows) and a GDM is fit to each
#' permuted table. Model significance is determined by comparing the deviance
#' explained by GDM fit to the un-permuted table to the distribution of deviance
#' explained values from GDM fit to the nPerm permuted tables. To assess predictor
#' significance, this process is repeated for each predictor individually (i.e.,
#' only the data for the predictor being tested is permuted rather than the entire
#' environmental table). Predictor importance is quantified as the percent change
#' in deviance explained between a model fit with and without that predictor permuted. If
#' predSelect=TRUE, this process continues by next permutating the site-pair
#' table nPerm times, but removing one predictor at a time and reassessing
#' predictor importance and significance. At each step, the least important
#' predictor is dropped (backward elimination) and the process continues until
#' all non-significant predictors are removed, with significance level being set
#' by the user and the pValue argument.
#'
#' @return A list of four tables. The first table summarizes full model deviance,
#' percent deviance explained by the full model, the p-value of the full model,
#' and the number of permutations used to calculate the statistics for each
#' fitted model (i.e., the full model and each model with predictors removed in
#' succession during the backward elimination procedure if predSelect=T). The
#' remaining three tables summarize (1) predictor importance, (2) predictor
#' significance, and (3) the number of permutations used to calculate the
#' statistics for that model, which is provided because some GDMs may fail
#' to converge for some permutations / predictor combinations and you might want to
#' know how many permutations were used when calculating statistics. Or maybe you
#' don't, you decide.
#'
#' Predictor importance is measured as the percent decrease in deviance explained
#' between the full model and the deviance explained by a model fit with that predictor
#' permuted. Significance is estimated using the bootstrapped p-value when the
#' predictor has been permuted. For most cases, the number of permutations will
#' equal the nPerm argument. However, the value may be less should any of the models
#' fit to them permuted tables fail to converge.
#'
#' If predSelect=FALSE, the tables will have values only in the first column.
#'
#' @author Matt Fitzpatrick and Karel Mokany
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
#' #modTest <- gdm.varImp(sitePairTab, geo=T, nPerm=50, parallel=T, cores=10, predSelect=T)
#' #barplot(sort(modTest$`Predictor Importance`[,1], decreasing=T))
#'
#' @keywords gdm
#'
#' @importFrom pbapply pbreplicate
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @importFrom doParallel registerDoParallel
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#'
#' @export
gdm.varImp <- function(spTable, geo, splines=NULL, knots=NULL, predSelect=FALSE,
                       nPerm=50, pValue=0.05, parallel=FALSE, cores=2, sampleSites=1,
                       sampleSitePairs=1, outFile=NULL){

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
  ##checks that predSelect has either TRUE or FALSE
  if(!(predSelect==TRUE | predSelect==FALSE)){
    stop("The predSelect argument must be either TRUE or FALSE.")
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

  # number of variables in the site-pair table, adds 1 if geo=TRUE
  nVars <- (ncol(spTable)-6)/2
  # create vector of variable names
  varNames <- colnames(spTable[c(7:(6+nVars))])
  varNames <- sapply(strsplit(varNames, "s1."), "[[", 2)
  if(geo==TRUE){
    nVars <- nVars + 1
    varNames <- c("Geographic", varNames)
  }

  # run initial GDM to see if any vars have zero I-spline coeffs
  message(paste0("Fitting initial model with all ", nVars,  " predictors..."))
  Sys.sleep(1)
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
    for(p in 1:length(zeroSum)){
      message(paste0("Sum of I-spline coefficients for predictor ", zeroSum[p]," = 0"))
      Sys.sleep(0.5)}
    #message("\n")
    for(p in 1:length(zeroSum)){
      if(zeroSum[p]=="Geographic"){
        message("Setting Geo=FALSE and proceeding with permutation testing...")
      } else {
        message(paste0("Removing ", zeroSum[p], " and proceeding with permutation testing..."))
      }
      Sys.sleep(0.5)
    }

    if(length(grep("Geographic", zeroSum))==1){
      geo <- FALSE
      zeroSum <- zeroSum[-grep("Geographic", zeroSum)]
    }

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

  if(cores>nVars){
    cores <- nVars
  }

  splines <- rep(unique(splines), nVars)

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

  ##sets up siteXvar table, uses for loop to make sure have steps correct
  for(i in 1:length(exBySite)){
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
  modelTestValues <- matrix(NA,4,nVars*5,dimnames = list(c("Model deviance",
                                                                    "Percent deviance explained",
                                                                    "Model p-value",
                                                                    "Fitted permutations"),
                                                                  c("All predictors",
                                                                    "1-removed",
                                                                    paste(seq(2,nVars*5-1), "-removed", sep=""))))
  ##deviance reduction predictor table
  varImpTable <- matrix(NA, nVars, nVars*5-1)
  rownames(varImpTable) <- varNames
  colnames(varImpTable) <- c("All predictors",
                             "1-removed",
                             paste(seq(2,nVars*5-2), "-removed", sep=""))
  ##p value predictor table
  pValues <- nModsConverge <- varImpTable

  ##assigns given site-pair table to new variable, to prevent changing the original input
  currSitePair <- spTable
  nullGDMFullFit <- 0  ##a variable to track rather or not the fully fitted gdm model returned a NULL object

  # create the set up permuted site-pair tables to be used for all
  # downstream analyses
  message(paste0("Creating ", nPerm, " permuted site-pair tables..."))

  if(parallel==F | nPerm <= 25){
    permSpt <- pbreplicate(nPerm, list(permutateSitePair(currSitePair,
                                                       siteData,
                                                       indexTab,
                                                       varNames)))}

  if(parallel==T & nPerm > 25){
    # set up parallel processing
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    permSpt <- foreach(k=1:nPerm,
                             .verbose=F,
                             .packages=c("gdm"),
                             .export=c("permutateSitePair")) %dopar%
      permutateSitePair(currSitePair, siteData, indexTab, varNames)
    stopCluster(cl)}

  varNames.x <- varNames
  message("Starting model assessment...")
  for(v in 1:length(varNames)){
    # ends the loop if only 1 variable remains
    if(length(varNames.x)<2){
      break
    }

    if(v>1){
      message(paste0("Removing ", names(elimVar), " and proceeding with the next round of permutations."))
    }

    if(is.numeric(splines)){
      splines <- rep(unique(splines), length(varNames.x))
    }

    # runs gdm, first time on site-pair table with all variables
    # As variables are eliminated the "full" site-pair table will
    # contain fewer variables
    fullGDM <- gdm(currSitePair, geo=geo, splines=splines, knots=knots)
    message(paste0("Percent deviance explained by the full model =  ", round(fullGDM$explained,3)))

    if(is.null(fullGDM)==TRUE){
      warning(paste("The model did not converge when testing variable: ", varNames.x[v],
                    ". Terminating analysis and returning output completed up to this point.", sep=""))
      break
    }

    message("Fitting GDMs to the permuted site-pair tables...")
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
      if(length(varNames.x)<cores){
        cores <- length(varNames.x)
      }
      # set up parallel processing
      cl <- makeCluster(cores)
      registerDoParallel(cl)

     # foreach function to create site-pair tables with each variable permuted,
      # fit gdms and extract deviance.
      permVarDev <- foreach(k=1:length(varNames.x),
                            .verbose=F,
                            .packages=c("gdm"),
                            .export = c("currSitePair")) %dopar%{
        if(varNames.x[k]!="Geographic"){
          # permute a single variable
          lll <- lapply(permSpt, function(x, spt=currSitePair){
            idx1 <- grep(paste("^s1.", varNames.x[k], "$", sep=""), colnames(x))
            idx2 <- grep(paste("^s2.", varNames.x[k], "$", sep=""), colnames(x))
            spt[,c(idx1, idx2)] <- x[,c(idx1, idx2)]
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
        return(permModelDev)
      }

      ##closes cores
      #close(pb)
      stopCluster(cl)
    }

    if(parallel==FALSE){
      # for-loop to create site-pair tables with each variable permuted,
      # fit gdms and extract deviance.
      permVarDev <- list()
      for(k in 1:length(varNames.x)){
        if(varNames.x[k]!="Geographic"){
          message(paste0("Assessing importance of ", varNames.x[k], "..."))
          # permute a single variable
          lll <- lapply(permSpt, function(x, spt=currSitePair){
            idx1 <- grep(paste("^s1.", varNames.x[k], "$", sep=""), colnames(x))
            idx2 <- grep(paste("^s2.", varNames.x[k], "$", sep=""), colnames(x))
            spt[,c(idx1, idx2)] <- x[,c(idx1, idx2)]
            return(spt)
          })}

        if(varNames.x[k]=="Geographic"){
          message("Assessing importance of geographic distance...")
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
        permVarDev[[k]] <- sapply(gdmPermVar, function(mod){mod$gdmdeviance})
      }
    }

    names(permVarDev) <- varNames.x
    nullDev <- fullGDM$nulldeviance

    for(var in varNames.x){
      grepper <- grep(paste0("^",var,"$"), names(permVarDev))
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
        noVarGDM <- gdm(testSitePair, geo=geo, splines=splines[-1], knots=knots)
      } else {
        noVarGDM <- gdm(currSitePair, geo=F, splines=splines[-1], knots=knots)
      }

      permDevReduct <- noVarGDM$gdmdeviance - varDevTab
      pValues[which(rownames(pValues) == var),v] <- sum(permDevReduct>=(varDevTab - fullGDM$gdmdeviance))/(nPerm-nConv)
    }

    if(max(na.omit(pValues[,v]))<pValue){
      message("All remaining predictors are significant, ceasing assessment.")
      message(paste0("Percent deviance explained by final model = ", round(fullGDM$explained,3)))
      message("Final set of predictors returned: ")
      for(vvv in 1:length(fullGDM$predictors)){
        message(fullGDM$predictors[vvv])
      }
      break
    }

    elimVar <- which.min(varImpTable[,v])

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

    if(v==1 & predSelect==F){
      message("Backwards elimination not selected by user (predSelect=F). Ceasing assessment.")
      message(paste0("Percent deviance explained by final model = ", round(fullGDM$explained,3)))
      message("Final set of predictors returned: ")
      for(vvv in 1:length(fullGDM$predictors)){
        message(fullGDM$predictors[vvv])
      }
      break
    }
  }


  if(v==1 & predSelect==F){
    # Model assessment
    modelTestVals <- data.frame(matrix(round(modelTestValues[,1], 3), ncol=1))
    rownames(modelTestVals) <- rownames(modelTestValues)
    colnames(modelTestVals) <- "All predictors"
    #Variable importance
    varImpTab <- data.frame(matrix(round(varImpTable[,1], 3), ncol=1))
    rownames(varImpTab) <- rownames(varImpTable)
    colnames(varImpTab) <- "All predictors"
    #Variable selection p-values
    pVals <- varImpTab
    pVals[,1] <- round(pValues[,1], 3)
    #Variable selection model convergence
    nModsConv <- varImpTab
    nModsConv[,1] <- round(nModsConverge[,1], 3)
    outObject <- list(modelTestVals, varImpTab, pVals, nModsConv)
    names(outObject) <- c("Model assessment", "Predictor Importance", "Predictor p-values", "Model Convergence")
  } else{
    outObject <- list(round(modelTestValues[,1:v], 3), round(varImpTable[,1:v],3), round(pValues[,1:v],3), nModsConverge[,1:v])
    names(outObject) <- c("Model assessment", "Predictor Importance", "Predictor p-values", "Model Convergence")
  }

  ##if given, writes out files to space on disk
  if(is.null(outFile)==FALSE){
    save(outObject, file=outFile)
  }
  return(outObject)
}

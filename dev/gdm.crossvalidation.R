#'@title GDM cross-validation test
#'
#'@description Undertake a cross-validation assessment of a GDM fit using all the predictors included in the formated GDM input
#' table (spTable). The cross-validation is run using a specified proportion (train.proportion) of the randomly selected sites
#' included in spTable to train the model, with the remaining sites being used to test the performance of the model predictions.
#' The test is repeated a specified number of times (n.crossvalid.tests), with a unique random sample taken each time. Outputs
#' are a number of cross-validation test metrics.
#'
#'@param spTable (dataframe) A dataframe holding the GDM input table for model fitting.
#'@param train.proportion (float) The proportion of sites in 'spTable' to use in training the GDM, with the remaining proportion used to test the model. (default = 0.9)
#'@param n.crossvalid.tests (integer) The number of cross-validation sets to use in testing the GDM. (default = 1)
#'@param geo (boolean) Geographic distance to be used in model fitting (default = FALSE) .
#'@param splines (vector) An optional vector of the number of I-spline basis functions to be used for each predictor in fitting the model.
#'@param knots (vector) An optional vector of knots in units of the predictor variables to be used in the fitting process.
#'
#'@return List, providing cross-validation statistics. These are metrics that describe how well the model fit using the
#' sitepair training table predicts the dissimilarities in the sitepair testing table. Metrics provided include:
#' 'Deviance.Explained' (the deviance explained for the training data);
#' 'Test.Deviance.Explained' (the deviance explained for the test data);
#' 'Mean.Error';
#' 'Mean.Absolute.Error';
#' 'Root.Mean.Squre.Error';
#' 'Obs.Pred.Correlation' (Pearson's correlation coefficient between observed and predicted values);
#' 'Equalised.RMSE' (the average root mean square error across bands of observed dissimilarities (0.05 dissimialrity units));
#' 'Error.by.Observed.Value' (the average root mean square error and number of observations within bands of observed dissimilarities (0.05 dissimialrity units)).
#'
#'@examples output = gdm_crossvalidation(My.GDM.input.table, n.crossvalid.tests=10)
#'
#'@export
gdm.crossvalidation=function(spTable,
                             train.proportion=0.9,
                             n.crossvalid.tests=1,
                             geo=FALSE,
                             splines=NULL,
                             knots=NULL)
{
  ##adds error checking to gdm function
  ##checks to see if in site-pair format from formatsitepair function
  if(class(spTable)[1] != "gdmData"){
    warning("spTable class does not include type 'gdmData'. Make sure your data is in site-pair format or the gdm model will not fit.")
  }
  ##checks to makes sure data is a matrix or data frame
  if(!(class(spTable)[1]=="gdmData" | class(spTable)[1]=="matrix" | class(spTable)[1]=="data.frame")){
    stop("spTable argument needs to be gdmData, a matrix, or a data frame")
  }
  ##makes sure that train.proportion is a number between 0 and 1,
  ##and that it is not equal to 0
  if(is.numeric(train.proportion)==FALSE | train.proportion<=0 | train.proportion>1){
    stop("argument train.proportion needs to be a positive number between 0 and 1")
  }
  ##Check we have at least one cross-validation test to run
  if(n.crossvalid.tests<1){
    stop("set 'n.crossvalid.tests' to 1 or greater")
  }


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
  numTestSites <- round(max(indexTab)*(1-train.proportion))
  numTrainSites <- max(indexTab) - numTestSites

  ##randomly determines the index of sites to remove
  if(numTestSites <= 1)
    {
    stop("train.proportion is too high - no sites are available as test data in the cross-validation.")
    }

  # Set up the catcher for the cross-validation outputs
  Deviance.Explained <- NULL
  Mean.Error <- NULL
  Mean.Absolute.Error <- NULL
  Root.Mean.Squre.Error <- NULL
  Obs.Pred.Correlation <- NULL
  Equalised.RMSE <- NULL
  Error.by.Observed.Value.npairs <- NULL
  Error.by.Observed.Value.RMSE <- NULL

  # Now loop through the cross-validation tests
  for(i.test in 1:n.crossvalid.tests) # turn this into a parallel loop perhaps
    {
    # randomly select the train and test sites
    testSites <- sample(1:max(indexTab), numTestSites)
    trainSites <- c(1:max(indexTab))
    trainSites <- trainSites[-testSites]
    # grab the sitepairs containing only training sites, and only testing sites
    # TEST
    rmIndexCol1 <- which(indexTab[,1] %in% testSites)
    rmIndexCol2 <- which(indexTab[,2] %in% testSites)
    all.test.indices <- c(rmIndexCol1,rmIndexCol2)
    test.pairs <- all.test.indices[duplicated(all.test.indices)]
    sampTableTest <- spTable[test.pairs,]
    # TRAIN
    rmIndexCol1 <- which(indexTab[,1] %in% trainSites)
    rmIndexCol2 <- which(indexTab[,2] %in% trainSites)
    all.test.indices <- c(rmIndexCol1,rmIndexCol2)
    train.pairs <- all.test.indices[duplicated(all.test.indices)]
    sampTableTrain <- spTable[train.pairs,]

    # Test how well a model fit with the training data predicts observed dissimilarities for the test data
    cv.test.out <- gdm.single.crossvalidation(sampTableTrain,
                                                   sampTableTest,
                                                   geo=geo,
                                                   splines=splines,
                                                   knots=knots)
    # Catch the outputs
    Deviance.Explained <- c(Deviance.Explained, cv.test.out$Test.Deviance.Explained)
    Mean.Error <- c(Mean.Error, cv.test.out$Mean.Error)
    Mean.Absolute.Error <- c(Mean.Absolute.Error, cv.test.out$Mean.Absolute.Error)
    Root.Mean.Squre.Error <- c(Root.Mean.Squre.Error, cv.test.out$Root.Mean.Squre.Error)
    Obs.Pred.Correlation <- c(Obs.Pred.Correlation, cv.test.out$Obs.Pred.Correlation)
    Equalised.RMSE <- c(Equalised.RMSE, cv.test.out$Equalised.RMSE)
    Error.by.Observed.Value.npairs <- cbind(Error.by.Observed.Value.npairs, cv.test.out$Error.by.Observed.Value$obs.count)
    Error.by.Observed.Value.RMSE <- cbind(Error.by.Observed.Value.RMSE, cv.test.out$Error.by.Observed.Value$pred.RMSE)

   }# end for i.test
  # Give the observed dissimilarity bands for the error by observed value
  row.names(Error.by.Observed.Value.npairs) <- cv.test.out$Error.by.Observed.Value$obs.dissim
  row.names(Error.by.Observed.Value.RMSE) <- cv.test.out$Error.by.Observed.Value$obs.dissim

  # Create a summary of the RMSE by observed value
  Error.by.Observed.Value <- data.frame('observed.dissimilarity' = cv.test.out$Error.by.Observed.Value$obs.dissim,
                                        'number.obs.sitepairs' = rowMeans(Error.by.Observed.Value.npairs, na.rm=TRUE),
                                        'RMSE' = rowMeans(Error.by.Observed.Value.RMSE, na.rm=TRUE))
  rownames(Error.by.Observed.Value) <- c()

  # Now generate outputs for the cross-validation
  # write the outputs of the function
  list(Deviance.Explained = mean(Deviance.Explained),
       Mean.Error = mean(Mean.Error),
       Mean.Absolute.Error = mean(Mean.Absolute.Error),
       Root.Mean.Squre.Error = mean(Root.Mean.Squre.Error),
       Obs.Pred.Correlation = mean(Obs.Pred.Correlation),
       Equalised.RMSE = mean(Equalised.RMSE),
       Error.by.Observed.Value = Error.by.Observed.Value,
       Full.Crossvalidation.Stats = list(Deviance.Explained = Deviance.Explained,
                                         Mean.Error = Mean.Error,
                                         Mean.Absolute.Error = Mean.Absolute.Error,
                                         Root.Mean.Squre.Error = Root.Mean.Squre.Error,
                                         Obs.Pred.Correlation = Obs.Pred.Correlation,
                                         Equalised.RMSE = Equalised.RMSE,
                                         Error.by.Observed.Value.npairs= Error.by.Observed.Value.npairs,
                                         Error.by.Observed.Value.RMSE = Error.by.Observed.Value.RMSE))

} # end gdm_crossvalidation()


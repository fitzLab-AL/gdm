#'@title Single GDM cross-validation test, core function
#'
#'@description Undertake a cross-validation assessment of a GDM, using a single
#'training and testing dataset.
#'
#' @usage gdm.single.crossvalidation(spTable_train, spTable_test, geo=FALSE,
#' splines=NULL, knots=NULL)
#'
#'@param spTable_train (dataframe) A dataframe holding the GDM input table for
#'model fitting.
#'
#'@param spTable_test (dataframe) A dataframe holding the GDM input table for
#'model testing, having identical column names to 'spTable_train' but using
#'different site-pairs.
#'
#'@param geo (boolean) Geographic distance to be used in model fitting
#'(default = FALSE).
#'
#'@param splines (vector) An optional vector of the number of I-spline basis
#'functions to be used for each predictor in fitting the model.
#'
#'@param knots (vector) An optional vector of knots in units of the predictor
#'variables to be used in the fitting process.
#'
#'@return List, providing cross-validation statistics. These are metrics that
#'describe how well the model fit using the sitepair training table predicts
#'the dissimilarities in the sitepair testing table. Metrics provided include:
#' 'Deviance.Explained' (the deviance explained for the training data);
#' 'Test.Deviance.Explained' (the deviance explained for the test data);
#' 'Mean.Error';
#' 'Mean.Absolute.Error';
#' 'Root.Mean.Squre.Error';
#' 'Obs.Pred.Correlation' (Pearson's correlation coefficient between observed and predicted values);
#' 'Equalised.RMSE' (the average root mean square error across bands of observed dissimilarities (0.05 dissimialrity units));
#' 'Error.by.Observed.Value' (the average root mean square error and number of observations within bands of observed dissimilarities (0.05 dissimialrity units)).
#'
#'@importFrom stats complete.cases
#'@importFrom stats cor
#'@importFrom stats na.omit
#'@importFrom stats predict
#'@importFrom stats sd
#'
#'@export
gdm.single.crossvalidation=function(spTable_train,
                                         spTable_test,
                                         geo=FALSE,
                                         splines=NULL,
                                         knots=NULL)
{
  # omit rows with NAs, just in case
  GDM_Table_Training <- spTable_train[complete.cases(spTable_train),]
  GDM_Table_Testing <- spTable_test[complete.cases(spTable_test),]

  # suppress the warnings GDM often throws up
  oldw <- getOption("warn")
  options(warn = -1)

  # Fit the gdm on the training set of sites
  train.mod <- gdm(GDM_Table_Training,
                        geo=geo,
                        splines=splines,
                        knots=knots)

  # now predict the dissimilarity for the test sites (pairs)
#NOTE - have to use predict.gdm for testing. Not sure this is appropriate or viable for package function
  pred.test <- predict(train.mod,
                            GDM_Table_Testing)

  # and fit a gdm to the test data to get the null deviance
  test.mod <- gdm(GDM_Table_Testing[,c(1:6)],
                       geo=TRUE)

  # reset the warnings
  options(warn = oldw)

  # Calculate deviance explained for the test data
  test.data.gdm.deviance <- calculate.gdm.deviance(pred.test, GDM_Table_Testing$distance)
  test.D2 <- ((test.mod$nulldeviance - test.data.gdm.deviance) / test.mod$nulldeviance)*100

  # calculate an array of test statistics using observed & predicted dissimilarities in the
  # testing data
  # Error
  error = GDM_Table_Testing[,1] - pred.test
  # Mean Error
  ME = mean(error)
  # Mean Absolute Error
  MAE = mean(abs(error))
  # Root Mean Square Error
  RMSE = sqrt(mean(error^2))
  # correlation, obs vs exp
  ObsPred.R <- cor(pred.test , GDM_Table_Testing[,1] , method = "pearson")
  # something new breaking error up over the range of observed values
  obs.breaks<-seq(from=0, to=1, by=0.05)
  obs.value.err <- data.frame("obs.dissim"=seq(from=0.025, to=0.975, by=0.05),
                              "obs.count"=rep(0, times=20),
                              "pred.RMSE"=rep(0, times=20))
  obs.err<-cbind(GDM_Table_Testing[,1], error)
  for(i.cat in 1:nrow(obs.value.err))
  {
    # find the observed values within the specified range of dissimilarities
    if(i.cat == 1){
      cat.obs.err<-obs.err[(obs.err[,1]>=obs.breaks[i.cat] & obs.err[,1]<=obs.breaks[(i.cat+1)]),2]
    }else{
      cat.obs.err<-obs.err[(obs.err[,1]>obs.breaks[i.cat] & obs.err[,1]<=obs.breaks[(i.cat+1)]),2]
    }
    # check if there are any observations in this category
    if(length(cat.obs.err)>0) {
      obs.value.err[i.cat,2]<-length(cat.obs.err)
      obs.value.err[i.cat,3]<-sqrt(mean(cat.obs.err^2))
    }else{
      obs.value.err[i.cat,3]<-NA
    }
  } #end for i.cat

  # summarise the prediction error across dissimilarity bands
  equ.RMSE<-mean(obs.value.err[,3], na.rm=TRUE)

  # write the outputs of the function
  list(Deviance.Explained = train.mod$explained,
       Test.Deviance.Explained = test.D2,
       Mean.Error = ME,
       Mean.Absolute.Error = MAE,
       Root.Mean.Squre.Error = RMSE,
       Obs.Pred.Correlation = ObsPred.R,
       Equalised.RMSE = equ.RMSE,
       Error.by.Observed.Value = obs.value.err)

} # end gdm.single.crossvalidation function

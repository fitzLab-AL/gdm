#' @title Calculate GDM Deviance for Observed & Predicted Dissimilarities
#'
#' @description Calculate GDM deviance for observed & predicted dissimilarities.
#'Can be used for assessing cross-validation data. Translated from the c++
#'function CalcGDMDevianceDouble() in the file NNLS_Double.cpp from the
#'GDM R package.
#'
#' @usage calculate.gdm.deviance(predDiss, obsDiss)
#'
#' @param predDiss (float) A vector of predicted dissimilarity values, of same
#'length as obsDiss.
#' @param obsDiss (float) A vector of observed dissimilarity values, of same
#'length as predDiss.
#'
#' @return A single value (float) being the deviance.
#'
#' @export
calculate.gdm.deviance <- function(predDiss, obsDiss){
  dTotal = 0
  for(i in 1:length(predDiss)){
    # calculate the first part (t1)
    if(predDiss[i] == 0){t1 <- obsDiss[i]} # end if
    else{
      if(obsDiss[i] == 0){t1 <- 0}
      # end if
      else{t1 <- obsDiss[i] * log(obsDiss[i] / predDiss[i])} # end else
    } # end else
    # calculate the second part (t2)
    if(predDiss[i] == 1){t2 <- 1 - obsDiss[i]} # end if
    else{
      if(obsDiss[i] == 1){t2 <- 1 - obsDiss[i]} # end if
      else{
        t2 <- (1.0 - obsDiss[i]) * log(( 1.0 - obsDiss[i]) / (1.0 - predDiss[i]))
      }# end else
    }# end else
    # accumulate the running sum
    dTotal <- dTotal+(t1 + t2)
  }# end for i.row
  dTotal <- dTotal * 2
  return(dTotal)
} # end calculate.gdm.deviance

#' @title Summarize a Fitted Generalized Dissimilarity Model
#'
#' @description This function summarizes the gdm model object returned from \code{\link[gdm]{gdm}}.
#'
#' @usage \method{summary}{gdm}(object, ...)
#'
#' @param object A gdm model object resulting from a call to \code{\link[gdm]{gdm}}.
#'
#' @param ... Ignored.
#'
#' @return summary prints its output to the R Console window and returns no value.
#'
#' @seealso \code{\link[gdm]{gdm}}
#'
#' @examples
#' ##set up site-pair table using the southwest data set
#' sppData <- southwest[, c(1,2,14,13)]
#' envTab <- southwest[, c(2:ncol(southwest))]
#' sitePairTab <- formatsitepair(sppData, 2, XColumn="Long", YColumn="Lat", sppColumn="species",
#'                               siteColumn="site", predData=envTab)
#'
#' ##create GDM
#' gdmMod <- gdm(sitePairTab, geo=TRUE)
#'
#' ##summary of GDM
#' summary(gdmMod)
#'
#' @keywords gdm
#'
#' @export
##function to print a summary of a gdm object
summary.gdm <- function (object, ...){
  print( "", quote=F )
  print( "", quote=F )
  print( "GDM Modelling Summary", quote=F );
  print( paste( "Creation Date: ", object$creationdate ), quote=F );
  print( "", quote=F )
  ##        call <- match.call()
  m <- match.call(expand.dots = F)
  print( paste( "Name: ", m[[2]] ), quote=F )
  print( "", quote=F )
  print( paste( "Data: ", object$dataname ), quote=F )
  print( "", quote=F )
  print( paste( "Samples: ", object$sample ), quote=F )
  print( "", quote=F )
  print( paste( "Geographical distance used in model fitting? ", object$geo ), quote=F )
  print( "", quote=F )
  print( paste( "NULL Deviance: ", round(object$nulldeviance,3) ), quote=F )
  print( paste( "GDM Deviance: ",  round(object$gdmdeviance,3) ), quote=F )
  print( paste( "Percent Deviance Explained: ", round(object$explained, 3)),
         quote=F )
  print( "", quote=F )
  print( paste( "Intercept: ", round(object$intercept, 3) ), quote=F )
  print( "", quote=F )
  print("PREDICTOR ORDER BY SUM OF I-SPLINE COEFFICIENTS:", quote=F)
  print( "", quote=F )

  #########################
  # reorder the predictors, splines, coeffs in order of
  # importance based on sum(coeffs)
  thiscoeff <- 1
  thisquant <- 1
  sumCoeff <- NULL
  for(i in 1:length(object$predictors)){
    numsplines <- object$splines[[i]]
    holdCoeff <- NULL
    for(j in 1:numsplines){
      holdCoeff[j] <- object$coefficients[[thiscoeff]]
      thiscoeff <- thiscoeff + 1
    }
    sumCoeff[i] <- sum(holdCoeff)
  }

  lateralus <- NULL
  schism <- NULL
  orderPreds <- order(sumCoeff, decreasing = T)
  for(op in orderPreds){
    parabol <- 1+(cumsum(object$splines)[op]-object$splines[op])#1+(op*object$splines[op]-(object$splines[op]))
    parabola <- cumsum(object$splines)[op]#op*object$splines[op]
    lateralus <- c(lateralus, object$coefficients[parabol:parabola])
    schism <- c(schism, object$knots[parabol:parabola])
  }

  object$predictors <- object$predictors[orderPreds]
  object$splines <- object$splines[orderPreds]
  object$coefficients <- lateralus
  object$knots <- schism
  object$sumCoeff <- sumCoeff[orderPreds]
  #########################

  for(i in 1:length(object$predictors)){
    print( paste( "Predictor ",i,": ",object$predictors[[i]], sep="" ), quote=F )
    print( paste( "Splines: ",object$splines[[i]], sep="" ), quote=F )
    numsplines <- object$splines[[i]]
    for(j in 1:numsplines){
      if ( j == 1 ) print( paste( "Min Knot: ", round(object$knots[[thisquant]], 3), sep="" ), quote=F )
      else if ( j == numsplines ) print( paste( "Max Knot: ", round(object$knots[[thisquant]], 3), sep="" ), quote=F )
      else print( paste( round(100/(numsplines-1),digits=2),"% Knot: ", round(object$knots[[thisquant]], 3), sep="" ), quote=F )
      thisquant <- thisquant + 1
    }
    #for(j in 1:numsplines){
    #  print( paste( "Coefficient[",j,"]: ", round(object$coefficients[[thiscoeff]], 3), sep="" ), quote=F )
    #  thiscoeff <- thiscoeff + 1
    #}
    #print( "", quote=F )
    print(paste0("Sum of coefficients for ", object$predictors[[i]], ": ",
                 round(object$sumCoeff[i], 3)), quote=F)
    print( "", quote=F )
  }
}

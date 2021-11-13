#' Summarize a Fitted Generalized Dissimilarity Model
#'
#' This function summarizes the gdm model object returned from \code{\link[gdm]{gdm}}.
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
  print( paste( "NULL Deviance: ", object$nulldeviance ), quote=F )
  print( paste( "GDM Deviance: ", object$gdmdeviance ), quote=F )
  print( paste( "Percent Deviance Explained: ", object$explained ), quote=F )
  print( "", quote=F )
  print( paste( "Intercept: ", object$intercept ), quote=F )
  print( "", quote=F )
  thiscoeff <- 1
  thisquant <- 1
  for(i in 1:length(object$predictors)){
    print( paste( "Predictor ",i,": ",object$predictors[[i]], sep="" ), quote=F )
    print( paste( "Splines: ",object$splines[[i]], sep="" ), quote=F )
    numsplines <- object$splines[[i]]
    for(j in 1:numsplines){
      if ( j == 1 ) print( paste( "Min Knot: ",object$knots[[thisquant]], sep="" ), quote=F )
      else if ( j == numsplines ) print( paste( "Max Knot: ",object$knots[[thisquant]], sep="" ), quote=F )
      else print( paste( round(100/(numsplines-1),digits=2),"% Knot: ",object$knots[[thisquant]], sep="" ), quote=F )
      thisquant <- thisquant + 1
    }
    for(j in 1:numsplines){
      print( paste( "Coefficient[",j,"]: ",object$coefficients[[thiscoeff]], sep="" ), quote=F )
      thiscoeff <- thiscoeff + 1
    }
    print( "", quote=F )
  }
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gdm.summary.R
\name{summary.gdm}
\alias{summary.gdm}
\title{Summarize a Fitted Generalized Dissimilarity Model}
\usage{
\method{summary}{gdm}(object, ...)
}
\arguments{
\item{object}{A gdm model object resulting from a call to \code{\link[gdm]{gdm}}.}

\item{...}{Ignored.}
}
\value{
summary prints its output to the R Console window and returns no value.
}
\description{
This function summarizes the gdm model object returned from \code{\link[gdm]{gdm}}.
}
\examples{
##set up site-pair table using the southwest data set
sppData <- southwest[, c(1,2,14,13)]
envTab <- southwest[, c(2:ncol(southwest))]
sitePairTab <- formatsitepair(sppData, 2, XColumn="Long", YColumn="Lat", sppColumn="species",
                              siteColumn="site", predData=envTab)

##create GDM
gdmMod <- gdm(sitePairTab, geo=TRUE)

##summary of GDM
summary(gdmMod)

}
\seealso{
\code{\link[gdm]{gdm}}
}
\keyword{gdm}

\name{createsitepair}
\alias{createsitepair}
\title{
Creates Site-Pair Table
}
\description{
Takes the lower half of a site-by-site distance (dissimilarity) matrix and creates the site-pair table. This function is called from the \code{\link{formatsitepair}} function and should not need to be called directly by the user.
}
\usage{
createsitepair(dist, spdata, envInfo, dXCol, dYCol, siteCol, weightsType, custWeights)
}
\arguments{
  \item{dist}{
The lower half of a site-by-site distance (dissimilarity) matrix, provided by the \code{\link{formatsitepair}} function.
}
  \item{spdata}{
Input species data, the same as the bioData input to the \code{\link{formatsitepair}} function.
}
  \item{envInfo}{
Input environmental data. Only accepts data tables as input. If the environmetal data for \code{\link{formatsitepair}} are rasters, the data would have been extracted into table format within \code{\link{formatsitepair}}.
}
  \item{dXCol}{
Input x coordinate, the same as the XColumn input to the \code{\link{formatsitepair}} function.
}
  \item{dYCol}{
Input y coordinate, the same as the YColumn input to the \code{\link{formatsitepair}} function.
}
  \item{siteCol}{
Site column, taken from either the species or environmental tables.
}
  \item{weightsType}{
The method of determining the site-pair weights.
}
  \item{custWeights}{
Custom wieghts, as a vector, if given by the user.
}
}
\value{
A site-pair table with appropriate distance (dissimilarity) and weight columns used for fitting GDM.
}
\note{
This function is called from the \code{\link{formatsitepair}} function and the user should not need to access directly.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link[gdm]{formatsitepair}}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ gdm }
\keyword{ internal }

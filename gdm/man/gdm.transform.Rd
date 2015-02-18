\name{gdm.transform}
\alias{gdm.transform}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transform Predictor Data Using A Generalized Dissimilarity Model
}
\description{
This function transforms a table of geographic and environmental predictors using (1) the fitted functions from a model object returned from \code{\link{gdm}} and (2) a data frame containing predictor data for a set of sites. Future versions will incoporate direct transformation of rasters stacks for mapping biological variation.
}
\usage{
gdm.transform(model, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
A gdm model object resulting from a call to \code{\link{gdm}}.
}
  \item{data}{
A data frame containing values for each predictor variable in the model, formatted as follows: X, Y, var1, var2, var3, ..., varN. If geographic distance was not used as a predictor in model fitting, the first two columns (X, Y) need to be removed. Note that the order of the predictor columns MUST be the same as the order of the predictors in the site-pair table used in model fitting. There is currently no checking to ensure that the order of the variables to be transformed are the same as those in the site-pair table used in model fitting.
}
}
\value{
gdm.transform returns a matrix with the same number of rows as the input data frame. If the model uses geographic distance as a predictor the output table will have the first and second columns set to the transformed X and Y values for each site in the input table. The transformed predictor data will be in the remaining columns. Values of the predictors can be taken from a set of rasters, that once transformed, can be mapped to reveal spatial biological variation. 
}
\references{
Ferrier S, Manion G, Elith J, Richardson, K (2007) Using generalized dissimilarity modelling to analyse and predict patterns of beta diversity in regional biodiversity assessment. \emph{Diversity & Distributions} 13, 252-264.

Fitzpatrick MC, Keller SR (2015) Ecological genomics meets community-level modeling of biodiversity: Mapping the genomic landscape of current and future environmental adaptation. \emph{Ecology Letters} 18: 1-16
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link[gdm]{predict.gdm}}}

\examples{
##sets up site-pair table
load(system.file("./data/gdm.RData", package="gdm"))
sppData <- gdmExpData[, c(1,2,13,14)]
envTab <- gdmExpData[, c(2:ncol(gdmExpData))]
sitePairTab <- formatsitepair(sppData, 2, XColumn="Long", YColumn="Lat", sppColumn="species", siteColumn="site", predData=envTab)

##create GDM
gdmMod <- gdm(sitePairTab, geo=TRUE)

##reordering environmental data to create table for transform function
envTrans <- envTab[, c(13,12,2:11)] # same order as sitePairTab

##tranform geo & env to biological importance
bioImp <- gdm.transform(gdmMod, envTrans)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ gdm }

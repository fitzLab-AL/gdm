\name{gdm.transform}
\alias{gdm.transform}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transform Environmental Data Using a Generalized Dissimilarity Model
}
\description{
This function transforms geographic and environmental predictors using (1) the fitted functions from a model object returned from \code{\link{gdm}} and (2) a data frame or raster stack containing predictor data for a set of sites.
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
Either (i) a data frame containing values for each predictor variable in the model, formatted as follows: X, Y, var1, var2, var3, ..., varN or (ii) a raster stack with one layer per predictor variable used in the model, excluding X and Y (rasters for x- and y-coordinates are built automatically from the input rasters if the model was fit with geo=T). The order of the columns (data frame) or raster layers (raster stack) MUST be the same as the order of the predictors in the site-pair table used in model fitting. There is currently no checking to ensure that the order of the variables to be transformed are the same as those in the site-pair table used in model fitting. If geographic distance was not used as a predictor in model fitting, the x- and y-columns need to be removed from the data to be transformed. Output is provided in the same format as the input data. 
}
}
\value{
gdm.transform returns either a data frame with the same number of rows as the input data frame or a raster stack, depending on the format of the input data. If the model uses geographic distance as a predictor the output object will contain columns or layers for the transformed X and Y values for each site. The transformed environmental data will be in the remaining columns or layers. 
}
\references{
Ferrier S, Manion G, Elith J, Richardson, K (2007) Using generalized dissimilarity modelling to analyse and predict patterns of beta diversity in regional biodiversity assessment. \emph{Diversity & Distributions} 13, 252-264.

Fitzpatrick MC, Keller SR (2015) Ecological genomics meets community-level modeling of biodiversity: Mapping the genomic landscape of current and future environmental adaptation. \emph{Ecology Letters} 18: 1-16
}

\examples{
load(system.file("./data/gdm.RData", package="gdm"))

# grab the columns with xy, site ID, and species data 
sppTab <- gdmExpData[, c("species", "site", "Lat", "Long")]

##fit gdm using rasters
rastFile <- system.file("./extdata/stackedVars.grd", package="gdm")
envRast <- stack(rastFile)
sitePairRast <- formatsitepair(sppTab, 2, XColumn="Long", YColumn="Lat", sppColumn="species",
	siteColumn="site", predData=envRast)
##remove NA values
sitePairRast <- na.omit(sitePairRast)

##fit raster GDM
gdmRastMod <- gdm(sitePairRast, geo=TRUE)

##raster input, raster output
transRasts <- gdm.transform(gdmRastMod, envRast)

# map biological patterns
rastDat <- sampleRandom(transRasts, 10000)
pcaSamp <- prcomp(rastDat)
 
# note the use of the 'index' argument
pcaRast <- predict(transRasts, pcaSamp, index=1:3)

# scale rasters
pcaRast[[1]] <- (pcaRast[[1]]-pcaRast[[1]]@data@min) /
  (pcaRast[[1]]@data@max-pcaRast[[1]]@data@min)*255
pcaRast[[2]] <- (pcaRast[[2]]-pcaRast[[2]]@data@min) /
  (pcaRast[[2]]@data@max-pcaRast[[2]]@data@min)*255
pcaRast[[3]] <- (pcaRast[[3]]-pcaRast[[3]]@data@min) /
  (pcaRast[[3]]@data@max-pcaRast[[3]]@data@min)*255

plotRGB(pcaRast, r=1, g=2, b=3)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ gdm }

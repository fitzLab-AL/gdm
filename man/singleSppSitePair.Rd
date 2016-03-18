\name{singleSppSitePair}
\alias{singleSppSitePair}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Converts Common Biological (Response) and Environmental (Predictor) Data Formats to a GDM Site-Pair Table for a Single Species
}
\description{
This function takes input biological and environmental data in various formats and builds a site-pair formatted table for a single species required for fitting a Generalized Dissimilarity Model using the \code{\link{gdm}} function. NOTE: Sample site x-y coordinates MUST be present in either the biological or the environmental data.
}
\usage{
singleSppSitePair(bioData, bioFormat, siteColumn = NULL, XColumn, YColumn, 
	sppColumn = NULL, sppID, abundColumn = NULL, predData, distPreds = NULL, 
	weightType = "equal", custWeightVect = NULL, samples = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bioData}{
The input biological (response) data table. Similar to the the bioData arugment in \code{\link{formatsitepair}}, but only accepts either a site-species matrix or a x, y, species list (table formats 1 and 2).
}
  \item{bioFormat}{
An integer code specifying the format of bioData. Acceptable values are 1 or 2 (see \code{\link{formatsitepair}}).
}
  \item{siteColumn}{
Same as the siteColumn argument of \code{\link{formatsitepair}}.
}
  \item{XColumn}{
Same as the XColumn argument of \code{\link{formatsitepair}}.
}
  \item{YColumn}{
Same as the YColumn argument of \code{\link{formatsitepair}}.
}
  \item{sppColumn}{
Same as the sppColumn argument of \code{\link{formatsitepair}}.
}
  \item{sppID}{
The name of the species for which the site-pair table will be created for.
}
  \item{abundColumn}{
Same as the abundColumn argument of \code{\link{formatsitepair}}.
}
  \item{predData}{
Same as the predData argument of \code{\link{formatsitepair}}.
}
  \item{distPreds}{
Same as the distPreds argument of \code{\link{formatsitepair}}.
}
  \item{weightType}{
Same as the abundColumn argument of \code{\link{formatsitepair}}.
}
  \item{custWeightVect}{
Similar to the the custWeightVect arugment in \code{\link{formatsitepair}}, but only accepts either "equal" or "custom" as inputs.
}
  \item{samples}{
Same as the samples argument of \code{\link{formatsitepair}}.
}
}
\details{
Similar to the site-pair output of \code{\link{formatsitepair}}. The main differance is that a site-pair table created from this function is specifically set to model a single species. The main differance between the dissimularity messurement of a multi-species and a single species site-pair matrix is that the single species matrix values are binary. A dissimularity value of 1 means that the species is present in both of the sites. A dissimularity value of 0 means that the species is present in one site and absent in the other. Site-pairs in which the species is absent in both sites are removed.
}
\value{
A site-pair formatted table containing the response (biological distance or dissimilarity), predictors, and weights as required for fitting Generalized Dissimilarity Models.
}

\seealso{\code{\link[gdm]{formatsitepair}}}
\examples{

##table data, species and environmental
load(system.file("./data/gdm.RData", package="gdm"))
sppData <- gdmExpData[c(1,2,13,14)]
envTab <- gdmExpData[c(2:ncol(gdmExpData))]

##environmental raster data
##commented out to reduce example run time
#rastFile <- system.file("./extdata/stackedVars.grd", package="gdm")
#envRast <- stack(rastFile)

#########table type 1
##site-species table without coordinates
testData1a <- reshape2::dcast(sppData, site~species)
##site-species table with coordinates
coords <- unique(sppData[2:ncol(sppData)])
testData1b <- merge(testData1a, coords, by="site")
##site-species, table-table
exFormat1a <- singleSppSitePair(testData1a, 1, siteColumn="site", XColumn="Long", YColumn="Lat",
	sppID="spp258", predData=envTab)
##site-species, table-raster
##commented out to reduce example run time
#exFormat1b <- singleSppSitePair(testData1b, 1, siteColumn="site", XColumn="Long", YColumn="Lat",
#	sppID="spp258", predData=envRast)

#########table type 2
##site xy spp list, table-table
exFormat2a <- singleSppSitePair(sppData, 2, XColumn="Long", YColumn="Lat", sppColumn="species",
	sppID="spp258", siteColumn="site", predData=envTab)
##site xy spp list, table-raster
##commented out to reduce example run time
#exFormat2b <- singleSppSitePair(sppData, 2, XColumn="Long", YColumn="Lat", sppColumn="species",
#	sppID="spp258", siteColumn="site", predData=envRast)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ gdm }


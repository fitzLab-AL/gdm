\name{removeSitesFromSitePair}
\alias{removeSitesFromSitePair}
\title{
Remove Random Sites from Site-Pair Table
}
\description{
Randomly selects a number of sites from a given site-pair table and removes them from the site-pair table. It will remove all instances of the sites randomly selected to be removed in both s1 and s2 positions.
}
\usage{
removeSitesFromSitePair(spTable, sampleSites)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{spTable}{
A site-pair table, same as used to fit a \code{\link{gdm}}.
}
  \item{sampleSites}{
The fraction (0-1, though a value of 0 would be silly, wouldn't it?) of \emph{sites to retain} from the full site-pair table. If less than 1, this argument will completely remove a fraction of sites such that they are not used in the permutation routines.
}
}
\value{
A site-pair table, such as one created by \code{\link{formatsitepair}}, ideally smaller than the one given. In the very rare case where the function determines not to remove any sites, or should the sampleSites argument be 1, then the function will return the full site-pair table.
}
\note{
This function removes sites, not just site-pairs (rows) from the site-pair table. This function is called from several of the other functions within the gdm package, including the \code{\link{plotUncertainty}} and \code{\link{gdm.varImp}} functions, for the purposes of subsampling the sites in the site-pair table.  
}

\seealso{\code{\link{formatsitepair}}}

\examples{
##fit table environmental data
##sets up site-pair table, environmental tabular data
load(system.file("./data/gdm.RData", package="gdm"))
sppData <- gdmExpData[c(1,2,13,14)]
envTab <- gdmExpData[c(2:ncol(gdmExpData))]
sitePairTab <- formatsitepair(sppData, 2, XColumn="Long", YColumn="Lat", sppColumn="species", 
	siteColumn="site", predData=envTab)

removeSitesFromSitePair(sitePairTab, sampleSites=0.7)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ gdm }


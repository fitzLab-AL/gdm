\name{plotUncertainty}
\alias{plotUncertainty}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot I-splines with error bands using bootstrapping.
}
\description{
This function estimates uncertainty in the fitted I-splines using bootstrapping. The function can run in parallel on multicore machines to reduce computation time (recommended for large number of iterations). I-spline plots with error bands (+/- one standard deviation) are produced showing (1) the variance of I-spline coefficients and (2) a rug plot indicating how sites used in model fitting are distributed along each gradient.
}
\usage{
plotUncertainty(spTable, sampleSites, bsIters, geo=FALSE, splines=NULL, knots=NULL, 
                splineCol="blue", errCol="grey80", plot.linewidth=2.0, plot.layout=c(2,2),
                parallel=FALSE, cores=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{spTable}{
A site-pair table, same as used to fit a \code{\link{gdm}}.
}
  \item{sampleSites}{
The fraction (0-1) of sites to retain from the full site-pair table when subsampling.
}
  \item{bsIters}{
The number of bootstrap iterations to perform.
}
  \item{geo}{
Same as the \code{\link{gdm}} geo argument.
}
  \item{splines}{
Same as the \code{\link{gdm}} splines argument.
}
  \item{knots}{
Same as the \code{\link{gdm}} knots argument.
}
  \item{splineCol}{
The color of the plotted mean spline. The default is "blue".
}
  \item{errCol}{
The color of shading for the error bands (+/- one standard deviation around the mean line). The default is "grey80".
}
  \item{plot.linewidth}{
The line width of the plotted mean spline line. The default is 2.
}
  \item{plot.layout}{
Same as the \code{\link{plot.gdm}} plot.layout argument.
}
  \item{parallel}{
Whether or not to run the uncertainty estimation in parallel. The parallel processing is highly recommended when the bsIters argument is hundreds or more. The default is FALSE (no parallel processing). 
}
  \item{cores}{
When the parallel argument is set to TRUE, the number of cores to be registered for the foreach loop. Must be <= the number of cores in the machine running the function.
}
}

\value{
plotUncertainty returns NULL.
}
\references{
Shryock, D. F., C. A. Havrilla, L. A. DeFalco, T. C. Esque, N. A. Custer, and T. E. Wood. 2015. Landscape genomics of \emph{Sphaeralcea ambigua} in the Mojave Desert: a multivariate, spatially-explicit approach to guide ecological restoration. \emph{Conservation Genetics} 16:1303-1317.
}

\seealso{\code{\link[gdm]{plot.gdm}, \link[gdm]{formatsitepair}, \link[gdm]{removeSitesFromSitePair}}}

\examples{
##sets up site-pair table
load(system.file("./data/gdm.RData", package="gdm"))
sppData <- gdmExpData[c(1,2,13,14)]
envTab <- gdmExpData[c(2:ncol(gdmExpData))]
sitePairTab <- formatsitepair(sppData, 2, XColumn="Long", YColumn="Lat", 
                              sppColumn="species", siteColumn="site", predData=envTab)

##plot GDM uncertainty using one core
#not run
#plotUncertainty(sitePairTab, sampleSites=0.70, bsIters=5, geo=TRUE, plot.layout=c(3,3))

##plot GDM uncertainty in parallel
#not run
#plotUncertainty(sitePairTab, sampleSites=0.70, bsIters=50, geo=TRUE, plot.layout=c(3,3),
#parallel=T, cores=10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ gdm }

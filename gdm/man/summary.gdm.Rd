\name{summary}
\alias{summary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarize a Fitted Generalized Dissimilarity Model
}
\description{
This function summarizes the model object returned from \code{\link{gdm}}.
}
\usage{
summary.gdm(model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
A gdm model object resulting from a call to \code{\link{gdm}}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
summary prints its output to the R Console window and returns no value.
}
\references{
%% ~put references to the literature/web site here ~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link[gdm]{gdm}}}
\examples{
##sets up site-pair table
load(system.file("./data/gdm.RData", package="gdm"))
sppData <- gdmExpData[, c(1,2,14,13)]
envTab <- gdmExpData[, c(2:ncol(gdmExpData))]
sitePairTab <- formatsitepair(sppData, 2, XColumn="Long", YColumn="Lat", sppColumn="species", siteColumn="site", predData=envTab)

##create GDM
gdmMod <- gdm(sitePairTab, geo=T)

##summary of GDM
summary(gdmMod)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ gdm }

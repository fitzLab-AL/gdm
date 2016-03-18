\name{removeSitesFromSitePair}
\alias{removeSitesFromSitePair}
\title{
Remove Sites from Site-Pair Table
}
\description{
Randomly selects a number of sites from a given site-pair table and removes them from the site-pair table. This function is called from the \code{\link{plotUncertainty}} function and should not need to be called directly by the user.
}
\usage{
removeSitesFromSitePair(n, tab, rmPer)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
An integer, not actually used, but needed to preform lapply
}
  \item{tab}{
A site-pair table, same as used to fit a \code{\link{gdm}}.
}
  \item{rmPer}{
The percent of sites to remove from the full site-pair table.
}
}
\value{
A site-pair table, smaller than the one given.
}
\note{
This function is called from the \code{\link{plotUncertainty}} function and the user should not need to access directly.
}

\seealso{\code{\link[gdm]{plotUncertainty}}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ gdm }
\keyword{ internal }

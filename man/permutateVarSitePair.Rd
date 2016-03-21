\name{permutateVarSitePair}
\alias{permutateVarSitePair}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Permutate the Values of a Site-Pair Table Variable
}
\description{
A function which randomizes the values of the given variables within a site-pair table. This function is called from the \code{\link{gdm.varImp}} function and should not need to be called directly by the user.
}
\usage{
permutateVarSitePair(spTab, siteVarTab, indexTab, vName)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{spTab}{
A given site-pair table.
}
  \item{siteVarTab}{
A given site x variable table. 
}
  \item{indexTab}{
A table of index values for the given site-pair table.
}
  \item{vName}{
Vector of variable names in both the site-pair and site x variable tables.
}
}
\value{
A site-pair table, with variable values randomized from the original site-pair table.
}
\note{
This function is called from the \code{\link{gdm.varImp}} function and the user should not need to access directly.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link[gdm]{gdm.varImp}}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ gdm }
\keyword{ internal }

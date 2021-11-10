\name{permutateSitePair}
\alias{permutateSitePair}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Permutate Site-Pair Table Rows
}
\description{
A function which randomizes the rows of a given site-pair table. This function is called from the \code{\link{gdm.varImp}} function and should not need to be called directly by the user.
}
\usage{
permutateSitePair(spTab, siteVarTab, indexTab, vNames)
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
  \item{vNames}{
Vector of variable names in both the site-pair and site x variable tables.
}
}
\value{
A site-pair table, with rows randomized from the original site-pair table.
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

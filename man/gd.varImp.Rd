\name{gdm.varImp}
\alias{gdm.varImp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate variable importance/significance in gdm using matrix permutation.
}
\description{
This function uses matrix permutation to perform variable significance testing in a generalized dissimilarity model. The function can run in parallel on multicore machines to reduce computation time (recommended). 
}
\usage{
gdm.varImp(spTable, geo, splines = NULL, knots = NULL, 
fullModelOnly = FALSE, nPerm = 100, parallel = FALSE, cores = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{spTable}{
A site-pair table, same as used to fit a \code{\link{gdm}}
}
  \item{geo}{
Similar to the \code{\link{gdm}} geo argument. The only difference is that the geo argument does not have a default in this function.
}
  \item{splines}{
Same as the \code{\link{gdm}} splines argument.
}
  \item{knots}{
Same as the \code{\link{gdm}} knots argument.
}
  \item{fullModelOnly}{
Set to TRUE to test only the full variable set. Set to false to estiatme variable importance and significance using matrix permutation and backward elimination. Default is FALSE.
}
  \item{nPerm}{
Number of permutations to use to estimate p-values. Default is 100.
}
  \item{parallel}{
Whether or not to run the matrix permutations and model fitting in parallel. The parallel processing is done using a foreach loop and it is highly recommended when the nPerms argument is hundreds or more. When is argument is set to FALSE, the processes are completed using lapply. The default is FALSE. 
}
  \item{cores}{
When the parallel argument is set to TRUE, the number of cores to be registered for the foreach loop. Must be <= the number of cores in the machine running the function.
}
}
\details{
This function implements matrix permutation to test variable significance in gdm as described in Ferrier et al. (2007) and Fitzpatrick et al. (2011). The function first fits a "full model" using all predictors in the site-pair table. Next, it permutates the site-pair table nPerm times by randomizing the order of the rows. A new set of gdm's are fit to these permutated site-pair tables to estimate an overall p-value for model significance. If fullModelOnly=F, the function then permutates the site-pair table nPerm times by randomizing the order of the rows, but one variable at a time. At each step, the least important variable is dropped (backward elimination) and the process continues until all variables have been tested.    
}
\value{
A list of three matrices. The first summarizes model deviance, percent deviance explained, and p-value for each fitted model (i.e., the full model and each model with variables removed in succession). If fullModelOnly=T, this table will have values only in the first column. The remaining two tables summarize variable importance and significance respectively. Variable importnace is measured as the amount full model deviacne is reduced when that variable is removed.   
}
\Author(s){
Karel Mokany, Matthew Lisk, and Matt Fitzpatrick
}
\references{
Ferrier S, Manion G, Elith J, Richardson, K (2007) Using generalized dissimilarity modelling to analyse and predict patterns of beta diversity in regional biodiversity assessment. \emph{Diversity & Distributions} 13, 252-264.

Fitzpatrick, MC, Sanders NJ, Ferrier S, Longino JT, Weiser MD, and RR Dunn. 2011. Forecasting the Future of Biodiversity: a Test of Single- and Multi-Species Models for Ants in North America. \emph{Ecography} 34: 836–47.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
##fit table environmental data
##sets up site-pair table, environmental tabular data
load(system.file("./data/gdm.RData", package="gdm"))
sppData <- gdmExpData[c(1,2,13,14)]
envTab <- gdmExpData[c(2:ncol(gdmExpData))]
sitePairTab <- formatsitepair(sppData, 2, XColumn="Long", YColumn="Lat", sppColumn="species", 
	siteColumn="site", predData=envTab)

## not run
#modTest <- gdm.varImp(sitePairTab, geo=T, nPerm=50, parallel=T, cores=10)
#barplot(modTest[[2]][,1])

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ gdm }

\name{gdm.varImp}
\alias{gdm.varImp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Quantify model significance and variable importance/significance in gdm using matrix permutation.
}
\description{
This function uses matrix permutation to perform model and variable significance testing and to estimate variable importance in a generalized dissimilarity model. The function can be run in parallel on multicore machines to reduce computation time (recommended until we learn to program in C++). 
}
\usage{
gdm.varImp(spTable, geo, splines = NULL, knots = NULL, fullModelOnly = FALSE, 
nPerm = 50, parallel = FALSE, cores = 2, sampleSites = 1, sampleSitePairs = 1,
outFile = NULL)
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
Set to TRUE to test only the full variable set. Set to FALSE to estimate model significance and variable importance and significance using matrix permutation and backward elimination. Default is FALSE.
}
  \item{nPerm}{
Number of permutations to use to estimate p-values. Default is 50.
}
  \item{parallel}{
Whether or not to run the matrix permutations and model fitting in parallel. Parallel processing is highly recommended when either (i) the nPerms argument is large (>100) or (ii) a large number of site-pairs (and or variables) are being used in model fitting (note computation demand can be reduced using subsampling - see next arguments). The default is FALSE. 
}
  \item{cores}{
When the parallel argument is set to TRUE, the number of cores to be registered for parallel processing. Must be <= the number of cores in the machine running the function.
}
  \item{sampleSites}{
The fraction (0-1, though a value of 0 would be silly, wouldn't it?) of \emph{sites to retain} from the full site-pair table. If less than 1, this argument will completely remove a fraction of sites such that they are not used in the permutation routines.
}
  \item{sampleSitePairs}{
The fraction (0-1) of \emph{site-pairs (i.e., rows) to retain} from the full site-pair table - in other words, all sites will be used in the permutation routines (assuming sampleSites = 1), but not all \emph{site-pair combinations}. In the case where both the sampleSites and the sampleSitePairs argument have values less than 1, sites first will be removed using the sampleSites argument, followed by removal of site-pairs using the sampleSitePairs argument. Note that the number of site-pairs removed is based on the fraction of the resulting site-pair table after sites have been removed, not on the size of the full site-pair table. 
}
  \item{outFile}{
An optional character string to write the object returned by the function to disk as an .RData object (".RData"" is not required as part of the file name). The .RData object will contain a single list with the name of "outObject". The default is NULL, meaning that no file will be written.
}
}
\details{
To test model significance, first a "full model" is fit using un-permuted environmental data. Next, the environmental data are permuted nPerm times (by randomizing the order of the rows) and a GDM is fit to each permuted table. Model significance is determined by comparing the deviance explained by GDM fit to the un-permuted table to the distribution of deviance explained values from GDM fit to the nPerm permuted tables. To assess variable significance, this process is repeated for each predictor individually (i.e., only the data for the variable being tested is permuted rather than the entire environmental table). Variable importance is quantified as the percent change in deviance explained between a model fit with and without that variable (technically speaking, with the variable permuted and un-permuted). If fullModelOnly=FALSE, this process continues by then permutating the site-pair table nPerm times, but removing one variable at a time and reassessing variable importance and significance. At each step, the least important variable is dropped (backward elimination) and the process continues until all non-significant predictors are removed. 
}
\value{
A list of four tables. The first table summarizes full model deviance, percent deviance explained by the full model, the p-value of the full model, and the number of permutations used to calculate the statistics for each fitted model (i.e., the full model and each model with variables removed in succession during the backward elimination procedure if fullModelOnly=F). The remaining three tables summarize (1) variable importance, (2) variable significance, and (3) the number of permutations used to calculate the statistics for that model, which is provided because some GDMs may fail to fit for some permutations / variable combinations and you might want to know how many permutations were used when clacuating statistics. Or maybe you don't, you decide. 

Variable importance is measured as the percent change in deviance explained by the full model and the deviance explained by a model fit with that variable permuted. Significance is estimated using the bootstrapped p-value when the variable has been permuted. For most cases, the number of permutations will equal the nPerm argument. However, the value may be less should any of the permutations fail to fit. 

If fullModelOnly=T, the tables will have values only in the first column and NAs elsewhere.    

NOTE: In some cases, GDM may fail to fit if there is a weak relationship between the response and predictors (e.g., when an important variable is removed). Such cases are indicated by -9999 values in the variable importance, variable significance, and number of permutations tables.
}
\author{
Karel Mokany, Matthew Lisk, and Matt Fitzpatrick
}
\references{
Ferrier S, Manion G, Elith J, Richardson, K (2007) Using generalized dissimilarity modelling to analyse and predict patterns of beta diversity in regional biodiversity assessment. \emph{Diversity & Distributions} 13, 252-264.

Fitzpatrick, MC, Sanders NJ, Ferrier S, Longino JT, Weiser MD, and RR Dunn. 2011. Forecasting the Future of Biodiversity: a Test of Single- and Multi-Species Models for Ants in North America. \emph{Ecography} 34: 836-47.
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
#barplot(sort(modTest[[2]][,1], decreasing=T))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ gdm }

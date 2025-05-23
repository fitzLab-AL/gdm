% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gdm.variable.importance.R
\name{gdm.varImp}
\alias{gdm.varImp}
\title{Assess Predictor Importance and Quantify Model Significance in a Fitted
Generalized Dissimilarity Model.}
\usage{
gdm.varImp(spTable, geo, splines = NULL, knots = NULL,
predSelect = FALSE, nPerm = 50, pValue=0.05, parallel = FALSE, cores = 2,
sampleSites = 1, sampleSitePairs = 1, outFile = NULL)
}
\arguments{
\item{spTable}{A site-pair table, same as used to fit a \code{\link[gdm]{gdm}}.}

\item{geo}{Similar to the \code{\link[gdm]{gdm}} geo argument. The only
difference is that the geo argument does not have a default in this function.}

\item{splines}{Same as the \code{\link[gdm]{gdm}} splines argument. Note that
the current implementation requires that all predictors have the same number of
splines.}

\item{knots}{Same as the \code{\link[gdm]{gdm}} knots argument.}

\item{predSelect}{Set to TRUE to perform predictor selection using matrix
permutation and backward elimination. Default is FALSE. When predSelect = FALSE
results will be returned only for a model fit with all predictors.}

\item{nPerm}{Number of permutations to use to estimate p-values. Default is 50.}

\item{pValue}{The p-value to use for predictor selection / elimination. Default is 0.05.}

\item{parallel}{Whether or not to run the matrix permutations and model
fitting in parallel. Parallel processing is highly recommended when either
(i) the nPerms argument is large (>100) or (ii) a large number of site-pairs
(and / or variables) are being used in model fitting (note computation demand
can be reduced using subsampling - see next arguments). The default is FALSE.}

\item{cores}{When the parallel argument is set to TRUE, the number of cores
to be registered for parallel processing. Must be <= the number of cores in
the machine running the function. There is no benefit to setting the number of
cores greater than the number of predictors in the model.}

\item{sampleSites}{The fraction (0-1, though a value of 0 would be silly,
wouldn't it?) of \emph{sites to retain} from the full site-pair table. If
less than 1, this argument will completely remove a fraction of sites such
that they are not used in the permutation routines.}

\item{sampleSitePairs}{The fraction (0-1) of \emph{site-pairs (i.e., rows)
to retain} from the full site-pair table - in other words, all sites will
be used in the permutation routines (assuming sampleSites = 1), but not
all \emph{site-pair combinations}. In the case where both the sampleSites
and the sampleSitePairs argument have values less than 1, sites first will
be removed using the sampleSites argument, followed by removal of site-pairs
using the sampleSitePairs argument. Note that the number of site-pairs
removed is based on the fraction of the resulting site-pair table after
sites have been removed, not on the size of the full site-pair table.}

\item{outFile}{An optional character string to write the object returned by
the function to disk as an .RData object (".RData" is not required as part
of the file name). The .RData object will contain a single list with the
name of "outObject". The default is NULL, meaning that no file will be written.}
}
\value{
A list of four tables. The first table summarizes full model deviance,
percent deviance explained by the full model, the p-value of the full model,
and the number of permutations used to calculate the statistics for each
fitted model (i.e., the full model and each model with predictors removed in
succession during the backward elimination procedure if predSelect=T). The
remaining three tables summarize (1) predictor importance, (2) predictor
significance, and (3) the number of permutations used to calculate the
statistics for that model, which is provided because some GDMs may fail
to converge for some permutations / predictor combinations and you might want to
know how many permutations were used when calculating statistics. Or maybe you
don't, you decide.

Predictor importance is measured as the percent decrease in deviance explained
between the full model and the deviance explained by a model fit with that predictor
permuted. Significance is estimated using the bootstrapped p-value when the
predictor has been permuted. For most cases, the number of permutations will
equal the nPerm argument. However, the value may be less should any of the models
fit to them permuted tables fail to converge.

If predSelect=FALSE, the tables will have values only in the first column.
}
\description{
This function uses matrix permutation to perform model and
predictor significance testing and to estimate predictor importance in a
generalized dissimilarity model. The function can be run in parallel on
multicore machines to reduce computation time.
}
\details{
To test model significance, first a model is fit using all predictors and
un-permuted environmental data. Any predictor for which the sum of the I-spline
coefficients sum to zero is preemptively removed. Next, the environmental data are permuted
nPerm times (by randomizing the order of the rows) and a GDM is fit to each
permuted table. Model significance is determined by comparing the deviance
explained by GDM fit to the un-permuted table to the distribution of deviance
explained values from GDM fit to the nPerm permuted tables. To assess predictor
significance, this process is repeated for each predictor individually (i.e.,
only the data for the predictor being tested is permuted rather than the entire
environmental table). Predictor importance is quantified as the percent change
in deviance explained between a model fit with and without that predictor permuted. If
predSelect=TRUE, this process continues by next permutating the site-pair
table nPerm times, but removing one predictor at a time and reassessing
predictor importance and significance. At each step, the least important
predictor is dropped (backward elimination) and the process continues until
all non-significant predictors are removed, with significance level being set
by the user and the pValue argument.
}
\examples{
##fit table environmental data
##set up site-pair table using the southwest data set
sppData <- southwest[, c(1,2,13,14)]
envTab <- southwest[, c(2:ncol(southwest))]
sitePairTab <- formatsitepair(sppData, 2, XColumn="Long", YColumn="Lat",
sppColumn="species", siteColumn="site", predData=envTab)

## not run
#modTest <- gdm.varImp(sitePairTab, geo=T, nPerm=50, parallel=T, cores=10, predSelect=T)
#barplot(sort(modTest$`Predictor Importance`[,1], decreasing=T))

}
\references{
Ferrier S, Manion G, Elith J, Richardson, K (2007) Using
generalized dissimilarity modelling to analyse and predict patterns of
beta diversity in regional biodiversity assessment. \emph{Diversity &
Distributions} 13, 252-264.

Fitzpatrick, MC, Sanders NJ, Ferrier S, Longino JT, Weiser MD, and RR Dunn.
2011. Forecasting the Future of Biodiversity: a Test of Single- and
Multi-Species Models for Ants in North America. \emph{Ecography} 34: 836-47.
}
\author{
Matt Fitzpatrick and Karel Mokany
}
\keyword{gdm}

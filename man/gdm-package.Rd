\name{gdm-package}

\alias{gdm-package}

\docType{package}

\title{
Overview of the functions in the gdm package
}

\description{
Generalized Dissimilarity Modeling is a statistical technique for modelling variation in biodiversity between pairs of geographical locations or through time. The \pkg{gdm} package provides functions to fit, evaluate, summarize, and plot Generalized Dissimilarity Models and to make predictions (across space and/or through time) and map biological patterns by transforming environmental predictor variables.
}
\details{
The functions in the \pkg{gdm} package provide the tools necessary for fitting GDMs, including functions to prepare biodiversity and environmental data. Major functionality includes:
  \itemize{
    \item Formatting various types of biodiversity and environmental data to \pkg{gdm}'s site-pair format used in model fitting
    \item Fitting GDMs using geographic and environmental distances between sites
    \item Plotting fitted functions & extracting I-spline values
    \item Estimating predictor importance using matrix permutation and predictor contributions using deviance paritioning
    \item Using cross-validation to evaluate models
    \item Predicting pairwise dissimiliarites between sites or times and transforming environmental predictors to biological importance and mapping these patterns.
  }

To see the preferable citation of the package, type \code{citation("gdm")}.
}

\section{I. Formatting input data}{

GDM fits biological distances to pairwise site geographical and environmental distances. Most users will need to first format their data to \pkg{gdm}'s site-pair table format:
  \tabular{ll}{
    \code{\link{formatsitepair}}\tab To convert biodiversity and environmental data to site-pair format \cr
}
}

\section{II. Model fitting, evaluation, and summary}{
  \tabular{ll}{
    \code{\link{gdm}}\tab To fit a GDM\cr
    \code{\link{gdm.crossvalidation}}\tab To evaluate a GDM\cr
    \code{\link{gdm.partition.deviance}}\tab To asses predictor contributions to deviance explained\cr
    \code{\link{gdm.varImp}}\tab To asses model significance and predictor importance\cr
    \code{\link{summary}}\tab To summarize a GDM\cr
}
}

\section{III. Model prediction and transformation of environmental data}{
  \tabular{ll}{
    \code{\link{predict}}\tab To predict biological dissimilarities between sites in space or between time periods\cr
    \code{\link{gdm.transform}}\tab To transform each environmental predictor to biological importance\cr
}
}

\section{IV. Plotting model output and fitted functions}{
  \tabular{ll}{
    \code{\link{plot}}\tab To plot model fit and I-splines\cr
    \code{\link{isplineExtract}}\tab To extract I-spline values to allow for custom plotting\cr
    \code{\link{plotUncertainty}}\tab To estimate and plot model sensitivity using bootstrapping\cr
}
}

\author{The \pkg{gdm} development team is Matt Fitzpatrick and Karel Mokany. The R package is based on code originally developed by Glenn Manion under the direction of Simon Ferrier. Where others have contributed to individual functions, credits are provided in function help pages.

The maintainer of the R version of \pkg{gdm} is Matt Fitzpatrick <mfitzpatrick@umces.edu>.
}

\name{predict}
\alias{predict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predict Biological Dissimilarities Between Sites Using a Generalized Dissimilarity Model
}
\description{
This function predicts response values (site-pair bioligcal dissimilarities) using a model object returned from \code{\link{gdm}} and a data frame containing the values of predictors for pairs of sites (these can be any sites, real or hypothetical, and are not restricted to those used to fit the model) and formatted as follows: Response, Weights, X0, Y0, X1, Y1, Pred1SiteA, Pred2SiteA, ..., PredNSiteA, Pred1SiteB, Pred2SiteB, ..., PredNSiteB.

Note that the first two columns (Response and Weights) are not actually used in generating predictions with the predict.gdm function, but is required and can be filled with dummy data (e.g. all zeros) if observed dissimilarities are not available. The second and third columns, X0 and Y0 represent the coordinates of the first site from a site pair. The fourth and fifth columns, X1 and Y1 represent the coordinates of the second site from the pair. Note that these columns are REQUIRED, even if geographic distance is not to be used a predictor. These columns can be loaded with dummy data if the actual coordinates are unknown. The next columns are for N predictors for SiteA and followed by N predictors for Site B. The columns must be in the same order and named the same as in the site-pair table used in model fitting. 
}
\usage{
predict(model, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
A gdm model object resulting from a call to gdm.
}
  \item{data}{
A data frame containing the values of predictors for pairs of sites, in the same format and structure as used to fit the model using \code{\link{gdm}}. The first two columns – Response and Weights – are ignored by predict.gdm and can therefore be filled with dummy data (e.g. all zeros). If geo is TRUE, then the X0, Y0 and Y0, Y1 columns will be used for calculating the geographic distance between each site in each site pair for inclusion of the geographic predictor term into the GDM model. If geo is FALSE, the default, then the X0, Y0, X1 and Y1 data columns are ignored (but these columns are still REQUIRED and can be filled with dummy data as needed). The predictor data for Site A and the predictor data for Site B follow.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
predict returns a response vector with the same length as the number of rows in the input data frame.
}
\references{
%% ~put references to the literature/web site here ~
}

\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link[gdm]{gdm.transform}}}

\examples{
##sets up site-pair table
load(system.file("./data/gdm.RData", package="gdm"))
sppData <- gdmExpData[, c(1,2,14,13)]
envTab <- gdmExpData[, c(2:ncol(gdmExpData))]
sitePairTab <- formatsitepair(sppData, 2, XColumn="Long", YColumn="Lat", sppColumn="species", siteColumn="site", predData=envTab)

##create GDM
gdmMod <- gdm(sitePairTab, geo=T)

##predict GDM
predDiss <- predict(gdmMod, sitePairTab)

##plot observed vs. predicted dissimilarity
plot(sitePairTab$distance, predDiss, xlim=c(0,1), ylim=c(0,1), col=rgb(0,0,1,0.5), pch=20,
xlab="Observed", ylab="Predicted")
lines(-1:2, -1:2, lwd=3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ gdm }

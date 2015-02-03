\name{gdm}
\alias{gdm}

\title{
Fit a Generalized Dissimilarity Model to Tabular Site-Pair Data
}
\description{

For an overview of the functions in the gdm package have a look here: \code{\link{gdm-package}}.

gdm is used to fit a generalized dissimilarity model to tabular site-pair data formatted as follows using the \code{\link{formatsitepair}} function: Response, Weights, X0, Y0, X1, Y1, Pred1SiteA, Pred2SiteA, ..., PredNSiteA, Pred1SiteB, Pred2SiteB, ..., PredNSiteB. The first column (Response) must be any ratio based dissimilarity (distance) measure between Site A and Site B. The second column defines any weighting to be applied during fitting of the model. If equal weighting is required, then all entries in this column should be set to 1.0 (default). The third and fourth columns, X0 and Y0, represent the spatial coordinates of the first site in the site pair (Site A). The fifth and sixth columns, X1 and Y1, represent the coordinates of the second site (Site B). Note that the first six columns are REQUIRED, even if you do not intend to use geographic distance as a predictor (in which case these columns can be loaded with dummy data if the actual coordinates are unknown). The next N*2 columns contain values for N predictors for Site A, followed by values for the same N predictors for Site B. 

The following is an example of a GDM input table header: 

Response, Weights, X0, Y0, X1, Y1, S1_Temp, S1_Rain, S1_Bedrock, S2_Temp, S2_Rain, S2_Bedrock
}
\usage{
gdm(data, geo=FALSE, splines=NULL, knots=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame containing the site pairs to be used to fit the GDM (typically obtained using the \code{\link{formatsitepair}} function). The observed response data must be located in the first column. The weights to be applied to each site pair must be located in the second column. If geo is TRUE, then the X0, Y0 and Y0, Y1 columns will be used to calculate the geographic distance between site pairs for inclusion as the geographic predictor term in the model. If geo is FALSE (default), then the X0, Y0, X1 and Y1 data columns must still be included, but are ignored in fitting the model. Columns containing the predictor data for Site A, and the predictor data for Site B, follow.
}
  \item{geo}{
Set to TRUE if geographic distance between sites is to be included as a model term. Set to FALSE if geographic distance is to be omitted from the model. Default is FALSE.
}
  \item{splines}{
An optional vector of the number of I-spline basis functions to be used for each predictor in fitting the model. If supplied, it must have the same length as the number of predictors (including geographic distance if geo is TRUE). If this vector is not provided (splines=NULL), then a default of 3 basis functions is used for all predictors.
}
  \item{knots}{
An optional vector of knots in  \emph{units of the predictor variables} to be used in the fitting process. If knots are supplied and splines=NULL, then the knots argument must have the same length as the number of predictors * 3. If both knots and the number of splines are supplied, then the length of the knots argument must be the same as the sum of the values in the splines vector. Note that the default values for knots when the default three I-spline basis functions are 0 (minimum), 50 (median), and 100 (maximum) quantiles.  
}
}
\details{
}
\value{
gdm returns a gdm model object. The function \code{\link{summary.gdm}} can be used to obtain or print a summary of the results. A gdm model object is a list containing at least the following components:
\item{}{}

\item{dataname}{ The name of the table used as the data argument to the model.}
\item{geo}{ Whether geographic distance was used as a predictor in the model.}
\item{gdmdeviance}{ The deviance of the fitted GDM model.}
\item{nulldeviance}{ The deviance of the null model.}
\item{explained}{ The percentage of null deviance explained by the fitted GDM model.}
\item{intercept}{ The fitted value for the intercept term in the model.}
\item{predictors}{ A list of the names of the predictors that were used to fit the model.}
\item{coefficients}{ A list of the coefficients for each spline for each of the predictors considered in model fitting.}
\item{knots}{ A vector of the knots derived from the x data (or user defined), for each predictor.}
\item{splines}{ A vector of the number of I-spline basis functions used for each predictor.}
\item{creationdate}{ The date and time of model creation.}
\item{observed}{ The observed response for each site pair (from data column 1).}
\item{predicted}{ The predicted response for each site pair, from the fitted model (after applying the link function).}
\item{ecological}{ The linear predictor (ecological distance) for each site pair, from the fitted model (before applying the link function).}
}
\references{
Ferrier S, Manion G, Elith J, Richardson, K (2007) Using generalized dissimilarity modelling to analyse and predict patterns of beta diversity in regional biodiversity assessment. \emph{Diversity & Distributions} 13, 252–264.
}


\seealso{\code{\link[gdm]{formatsitepair}, \link[gdm]{summary.gdm}, 
\link[gdm]{plot.gdm}, \link[gdm]{predict.gdm}, \link[gdm]{gdm.transform}} }

\examples{
##sets up site-pair table
load(system.file("./data/gdm.RData", package="gdm"))
sppData <- gdmExpData[c(1,2,13,14)]
envTab <- gdmExpData[c(2:ncol(gdmExpData))]
sitePairTab <- formatsitepair(sppData, 2, XColumn="Long", YColumn="Lat", sppColumn="species", siteColumn="site", predData=envTab)

##fit GDM
gdmMod <- gdm(sitePairTab, geo=TRUE)
summary(gdmMod)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ gdm }

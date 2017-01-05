\name{gdm}
\alias{gdm}

\title{
Fit a Generalized Dissimilarity Model to Tabular Site-Pair Data
}
\description{

For an overview of the functions in the gdm package have a look here: \code{\link{gdm-package}}.

The gdm function is used to fit a generalized dissimilarity model to tabular site-pair data formatted as follows using the \code{\link{formatsitepair}} function: distance, weights, s1.xCoord, s1.yCoord, s2.xCoord, s2.yCoord, s1.Pred1, s1.Pred2, ..., s1.PredN, s2.Pred1, s2.Pred2, ..., s2.PredN. The first column (Response) must be any ratio based dissimilarity (distance) measure between Site 1 and Site 2. The second column defines any weighting to be applied during fitting of the model. If equal weighting is required, then all entries in this column should be set to 1.0 (default). The third and fourth columns, s1.X and s1.Y, represent the spatial coordinates of the first site in the site pair (s1). The fifth and sixth columns, s2.X and s2.Y, represent the coordinates of the second site (s2). Note that the first six columns are REQUIRED, even if you do not intend to use geographic distance as a predictor (in which case these columns can be loaded with dummy data if the actual coordinates are unknown - though that would be weird, no?). The next N*2 columns contain values for N predictors for Site 1, followed by values for the same N predictors for Site 2. 

The following is an example of a GDM input table header with three environmental predictors (Temp, Rain, Bedrock): 

distance, weights, s1.xCoord, s1.yCoord, s2.xCoord, s2.yCoord, s1.Temp, s1.Rain, s1.Bedrock, s2.Temp, s2.Rain, s2.Bedrock
}
\usage{
gdm(data, geo=FALSE, splines=NULL, knots=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame containing the site pairs to be used to fit the GDM (typically obtained using the \code{\link{formatsitepair}} function). The observed response data must be located in the first column. The weights to be applied to each site pair must be located in the second column. If geo is TRUE, then the s1.xCoord, s1.yCoord and s2.xCoord, s2.yCoord columns will be used to calculate the geographic distance between site pairs for inclusion as the geographic predictor term in the model. If geo is FALSE (default), then the s1.xCoord, s1.yCoord, s2.xCoord and s2.yCoord data columns must still be included, but are ignored in fitting the model. Columns containing the predictor data for Site 1, and the predictor data for Site 2, follow.
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
\value{
gdm returns a gdm model object. The function \code{\link{summary.gdm}} can be used to obtain or print a synopsis of the results. A gdm model object is a list containing at least the following components:
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
Ferrier S, Manion G, Elith J, Richardson, K (2007) Using generalized dissimilarity modelling to analyse and predict patterns of beta diversity in regional biodiversity assessment. \emph{Diversity & Distributions} 13, 252-264.
}


\seealso{\code{\link[gdm]{formatsitepair}, \link[gdm]{summary.gdm}, 
\link[gdm]{plot.gdm}, \link[gdm]{predict.gdm}, \link[gdm]{gdm.transform}} }

\examples{
##fit table environmental data
##sets up site-pair table, environmental tabular data
load(system.file("./data/gdm.RData", package="gdm"))
sppData <- gdmExpData[c(1,2,13,14)]
envTab <- gdmExpData[c(2:ncol(gdmExpData))]
sitePairTab <- formatsitepair(sppData, 2, XColumn="Long", YColumn="Lat", sppColumn="species", 
	siteColumn="site", predData=envTab)

##fit table GDM
gdmTabMod <- gdm(sitePairTab, geo=TRUE)
summary(gdmTabMod)

##fit raster environmental data
##sets up site-pair table
rastFile <- system.file("./extdata/stackedVars.grd", package="gdm")
envRast <- stack(rastFile)

##environmental raster data
sitePairRast <- formatsitepair(sppData, 2, XColumn="Long", YColumn="Lat", sppColumn="species",
	siteColumn="site", predData=envRast)
##sometimes raster data returns NA in the site-pair table, these rows will have to be removed 
##before fitting gdm
sitePairRast <- na.omit(sitePairRast)

##fit raster GDM
gdmRastMod <- gdm(sitePairRast, geo=TRUE)
summary(gdmRastMod)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ gdm }

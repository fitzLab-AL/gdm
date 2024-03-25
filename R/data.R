#' @title Species and Environmental Data from Southwestern Australia.
#'
#' @description A data set containing species occurrence and associated environmental
#' data at 94 sites in southwestern Australia.
#'
#' @format A data frame with 29364 rows and 14 variables:
#' \describe{
#'   \item{species}{species name}
#'   \item{site}{site name}
#'   \item{awcA}{plant-available water capacity in soil horizon A}
#'   \item{phTotal}{soil pH}
#'   \item{sandA}{percent sand content in soil horizon A}
#'   \item{shcA}{saturated hydraulic conductivity in soil horizon A}
#'   \item{solumDepth}{soil depth to unweathered parent material}
#'   \item{bio5}{maximum temperature of the coldest month}
#'   \item{bio6}{minimum temperature of the coldest month}
#'   \item{bio15}{precipitation seasonality}
#'   \item{bio18}{precipitation of warmest quarter}
#'   \item{bio19}{precipitation of coldest quarter}
#'   \item{Lat}{latitude}
#'   \item{Long}{longitude}
#' }
"southwest"

#' An example biological dissimilarity matrix
#'
#' Pairwise Bray-Curtis dissimilarity calculated using the species
#' occurrence data from the \code{\link[gdm]{southwest}} data set.
#'
#' @format A data frame with 94 rows and 95 columns (extra column holds site IDs):
"gdmDissim"

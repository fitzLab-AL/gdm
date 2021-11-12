#' Species and environmental data from southwestern Australia.
#'
#' A data set containing species occurrence and associated environmental
#' data at 94 sites in southwestern Australia.
#'
#' @format A data frame with 29364 rows and 14 variables:
#' \describe{
#'   \item{species}{species name}
#'   \item{site}{site name}
#'   \item{awcA}{weight of the diamond, in carats}
#'   \item{phTotal}{weight of the diamond, in carats}
#'   \item{sandA}{weight of the diamond, in carats}
#'   \item{shcA}{weight of the diamond, in carats}
#'   \item{solumDepth}{weight of the diamond, in carats}
#'   \item{bio5}{weight of the diamond, in carats}
#'   \item{bio6}{weight of the diamond, in carats}
#'   \item{bio15}{weight of the diamond, in carats}
#'   \item{bio18}{weight of the diamond, in carats}
#'   \item{bio19}{weight of the diamond, in carats}
#'   \item{Lat}{weight of the diamond, in carats}
#'   \item{Long}{weight of the diamond, in carats}
#' }
"southwest"

#' An example biological dissimilarity matrix
#'
#' Pairwise Bray-Curtis dissimilarity calculated using the species
#' occurrence data from the \code{\link[gdm]{southwest}} data set.
#'
#' @format A data frame with 94 rows and 94 columns:
"gdmDissim"

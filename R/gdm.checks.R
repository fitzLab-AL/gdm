# Date : March 2024
# Version 0.1
# Licence GPL v3

# is it a raster object
.is_raster <- function(x){
  z <- class(x)
  return(
    z %in% c("SpatRaster", "RasterStack", "RasterLayer", "RasterBrick", "stars")
  )
}

# check for r
.check_rast <- function(r, name = "r"){
  if(!methods::is(r, "SpatRaster")){
    tryCatch(
      {
        r <- terra::rast(r)
      },
      error = function(cond) {
        message(sprintf("'%s' is not convertible to a terra SpatRaster object!", name))
        message(sprintf("'%s' must be a SpatRaster, stars, or Raster* object.", name))
      }
    )
  }
  return(r)
}


# check for required packages
.check_pkgs <- function(pkg){
  pkgna <- names(which(sapply(sapply(pkg, find.package, quiet = TRUE), length) == 0))
  if(length(pkgna) > 0){
    nm <- paste(pkgna, collapse = ", ")
    message("This function requires these packages: ", nm, "\nWould you like to install them now?\n1: yes\n2: no")
    user <- readline(prompt = paste0("Selection: "))
    if(tolower(user) %in% c("1", "yes", "y")){
      utils::install.packages(pkgna)
    } else{
      stop("Please install these packages for function to work: ", nm)
    }
  }
}

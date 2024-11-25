library(gdm)
library(terra)

sppTab <- southwest[, c("species", "site", "Lat", "Long")]

rastFile <- system.file("./extdata/swBioclims.grd", package="gdm")
envRast <- terra::rast(rastFile)


# test formatsitepair function --------------------------------------------

# expect warning because of using raster inputs
expect_warning({

  sitePairRast <- formatsitepair(
    sppTab,
    2,
    XColumn = "Long",
    YColumn = "Lat",
    sppColumn = "species",
    siteColumn = "site",
    predData = envRast
  )

})

sitePairRast <- na.omit(sitePairRast)

# build the model ---------------------------------------------------------

expect_silent({

  gdmRastMod <- gdm(sitePairRast, geo = TRUE)

})

# test predict function ---------------------------------------------------
futureRast <- envRast
futureRast[[3]] <- futureRast[[3]] * 0.75

# check prediction is done silently i.e. no errors/warnings
expect_silent({

  pred_raster <- gdm:::predict.gdm(
    object = gdmRastMod,
    data = envRast,
    time = TRUE,
    predRasts = futureRast,
    filename = tempfile(fileext = ".tif")
  )

})

# read the prediction created by gdm version 1.5.x
pred_v1.5 <- terra::rast(
  system.file("./extdata/test_data/pred_3_75.tif", package="gdm")
)

# round as file compression might lose some precision
diff <- round(terra::global(pred_raster - pred_v1.5, "mean", na.rm = TRUE)[1,1], 5)

expect_true(
  diff == 0
)


# test transform function -------------------------------------------------

# check gdm.transform is done silently i.e. no errors/warnings
expect_silent({

  trans_raster <- gdm.transform(
    model = gdmRastMod,
    data = envRast,
    filename = tempfile(fileext = ".tif")
  )

})

# read the transformed vars created by gdm version 1.5.x
trans_v1.5 <- terra::rast(
  system.file("./extdata/test_data/gdm_trans_vars.grd", package="gdm")
)

# round as file compression might lose some precision
diff <- round(terra::global(trans_raster - trans_v1.5, "mean", na.rm = TRUE)[, 1], 5)

# expect all difference be zero
expect_true(
  all(diff == 0)
)


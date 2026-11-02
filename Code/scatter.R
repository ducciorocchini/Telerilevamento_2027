
#### Example use ####

# Required packages

library(terra)
library(imageRy)

# Data import

ndvi_series <- im.import("Sentinel2_NDVI")

# Set names

names(ndvi_series) <- c("Feb2020", "May2020", "Aug2020", "Nov2020")

# Create classification for time series

ndvi_class <- im.classify(
  ndvi_series,
  num_clusters = 4,
  seed = 42,
  do_plot = TRUE
)

# Create scatterplot

ndvi_scatter <- im.scatter(
  band_x = ndvi_series[["Feb2020"]],
  band_y = ndvi_series[["Aug2020"]],
  classified_image = ndvi_class,
  n_sample = 1000,
  seed = 42,
  marginal = TRUE,
  legend = TRUE,
  title = "Seasonal NDVI variation by cluster",
  xlab = "Winter NDVI",
  ylab = "Summer NDVI"
)

ndvi_scatter

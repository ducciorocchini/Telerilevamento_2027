# Parallel sets

library(terra)

# Import data
fires <- rast("https://raw.githubusercontent.com/ducciorocchini/Telerilevamento_2027/main/Data/morrone_26aug2017.tif")
pre_fires <- rast("https://raw.githubusercontent.com/ducciorocchini/Telerilevamento_2027/main/Data/morrone_6aug2017.tif")
post_fires <- rast("https://raw.githubusercontent.com/ducciorocchini/Telerilevamento_2027/main/Data/morrone_15sept2017.tif")

# Create list of rasters
fire_list <- list(pre_fires, fires, post_fires)
names(fire_list) <- c("pre_fires", "fires", "post_fires")

# Create classification for time series
fire_class <- im.classify.timeseries(
  input_images = fire_list,
  num_clusters = 4,
  seed = 42,
  do_plot = TRUE
)

fire_class

# Run the function
im.parallel(fire_class)

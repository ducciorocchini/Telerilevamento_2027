# ## PARALLEL SETS ##

#### Required packages ####

library(terra)
library(imageRy)
library(tidyverse)
library(ggforce)

#### Data loading and cropping ####
# We are going to load three images of the Monte Morrone area, which was affected 
# by a large fire in August 2017

# pre_files <- list.files("./images/", pattern = "20170806", full.names = TRUE)
# files <- list.files("./images/", pattern = "20170826", full.names = TRUE)
# post_files <- list.files("./images/", pattern = "20170915", full.names = TRUE)
# 
# pre_fires <- rast(pre_files)
# fires <- rast(files)
# post_fires <- rast(post_files)
# 
# aoi <- ext(
#   13.955727, 14.061985,
#   42.033994, 42.110958
# )
# 
# aoi <- as.polygons(aoi)
# crs(aoi) <- "EPSG:4326"
# 
# # Project AOI to Sentinel-2 CRS
# aoi <- project(aoi, crs(pre_fires))
# 
# # Crop to aoi
# 
# fire_list <- list(pre_fires, fires, post_fires)
# 
# # Select bands and crop to AOI
# 
# bands <- c("B02", "B03", "B04", "B8A", "B11", "B12")
# 
# fire_list <- fire_list |>
#   purrr::map(
#     \(x) {
# 
#       x <- x[[
#         grep(
#           paste0("_(", paste(bands, collapse = "|"), ")_20m$"),
#           names(x)
#         )
#       ]]
# 
#       terra::crop(x, aoi)
#     }
#   )
# 
# names(fire_list) <- c("pre_fires", "fires", "post_fires")
# 
# list2env(fire_list, envir = .GlobalEnv)
# 
# # Save cropped rasters
# 
# writeRaster(pre_fires, "morrone_6aug2017.tif")
# writeRaster(fires, "morrone_26aug2017.tif")
# writeRaster(post_fires, "morrone_15sept2017.tif")

#### Loading cropped files ####

pre_fires <- rast("morrone_6aug2017.tif")
fires <- rast("morrone_26aug2017.tif")
post_fires <- rast("morrone_15sept2017.tif")

fire_list <- list(pre_fires, fires, post_fires)
names(fire_list) <- c("pre_fires", "fires", "post_fires")

#### Image classification ####

pre_c <- im.classify(pre_fires, num_clusters = 4, seed = 42)
fires_c <- im.classify(fires, num_clusters = 4, seed = 42)
post_c <- im.classify(post_fires, num_clusters = 4, seed = 42)

im.multiframe(1, 3)

plot(pre_c)
plot(fires_c)
plot(post_c)

#### BUILD FUNCTION FOR TIME SERIES CLASSIFICATION ####

im.classify.timeseries <- function(input_images,
                                   num_clusters = 3,
                                   seed = NULL,
                                   do_plot = TRUE,
                                   custom_colors = NULL,
                                   num_colors = 100) {
  
  # Check input images
  if (!is.list(input_images)) {
    stop("input_images should be a list of SpatRaster objects.")
  }
  
  if (length(input_images) < 2) {
    stop("input_images should contain at least two SpatRaster objects.")
  }
  
  if (!all(vapply(
    input_images,
    inherits,
    logical(1),
    what = "SpatRaster"
  ))) {
    stop("All elements of input_images should be SpatRaster objects.")
  }
  
  # Check number of layers
  n_layers <- vapply(
    input_images,
    terra::nlyr,
    numeric(1)
  )
  
  if (length(unique(n_layers)) != 1) {
    stop("All input images should have the same number of layers.")
  }
  
  # Check geometry
  same_geometry <- vapply(
    input_images[-1],
    function(x) {
      terra::compareGeom(
        input_images[[1]],
        x,
        stopOnError = FALSE
      )
    },
    logical(1)
  )
  
  if (!all(same_geometry)) {
    stop("All input images should have the same geometry.")
  }
  
  # Check number of clusters
  if (!is.numeric(num_clusters) ||
      length(num_clusters) != 1 ||
      is.na(num_clusters) ||
      num_clusters < 2 ||
      num_clusters %% 1 != 0) {
    
    stop("num_clusters must be an integer greater than or equal to 2.")
  }
  
  num_clusters <- as.integer(num_clusters)
  
  # Check seed
  if (!is.null(seed)) {
    
    if (!is.numeric(seed) ||
        length(seed) != 1 ||
        is.na(seed) ||
        seed %% 1 != 0) {
      
      stop("seed must be NULL or a single integer.")
    }
    
    seed <- as.integer(seed)
  }
  
  # Check do_plot
  if (!is.logical(do_plot) ||
      length(do_plot) != 1 ||
      is.na(do_plot)) {
    
    stop("do_plot must be either TRUE or FALSE.")
  }
  
  # Check custom colors
  if (!is.null(custom_colors)) {
    
    if (!is.character(custom_colors)) {
      stop("custom_colors must be NULL or a character vector of colors.")
    }
    
    if (length(custom_colors) < 2) {
      stop("custom_colors should contain at least two colors.")
    }
  }
  
  # Check number of colors
  if (!is.numeric(num_colors) ||
      length(num_colors) != 1 ||
      is.na(num_colors) ||
      num_colors < 1 ||
      num_colors %% 1 != 0) {
    
    stop("num_colors must be a positive integer.")
  }
  
  num_colors <- as.integer(num_colors)
  
  
  # Set colors
  if (is.null(custom_colors)) {
    class_colors <- viridisLite::viridis(num_colors)
  } else {
    class_colors <- grDevices::colorRampPalette(custom_colors)(num_colors)
  }
  
  # Extract raster values
  image_values <- lapply(
    input_images,
    terra::values,
    mat = TRUE
  )
  
  # Identify complete pixels
  complete_pixels <- lapply(
    image_values,
    stats::complete.cases
  )
  
  # Remove incomplete pixels
  image_values_complete <- Map(
    function(values, complete) {
      values[complete, , drop = FALSE]
    },
    image_values,
    complete_pixels
  )
  
  # Count complete pixels
  n_complete <- vapply(
    image_values_complete,
    nrow,
    integer(1)
  )
  
  # Check complete pixels
  if (any(n_complete == 0)) {
    stop("Each input raster must contain at least one complete pixel observation.")
  }
  
  # Combine pixel values from all images
  all_values_complete <- do.call(
    rbind,
    image_values_complete
  )
  
  if (nrow(all_values_complete) < num_clusters) {
    stop(
      paste(
        "The number of complete raster cells must be",
        "greater than or equal to num_clusters."
      )
    )
  }
  
  # Set seed
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Perform common k-means classification
  kmeans_result <- stats::kmeans(
    x = all_values_complete,
    centers = num_clusters
  )
  
  # Split cluster assignments among images
  cluster_groups <- rep(
    seq_along(n_complete),
    times = n_complete
  )
  
  clusters <- split(
    kmeans_result$cluster,
    cluster_groups
  )
  
  # Build classified images
  classified_images <- Map(
    function(image, complete, cluster) {
      
      classified_image <- image[[1]]
      
      classified_values <- rep(
        NA_integer_,
        terra::ncell(classified_image)
      )
      
      classified_values[complete] <- cluster
      
      terra::values(classified_image) <- classified_values
      
      names(classified_image) <- "class"
      
      classified_image
    },
    input_images,
    complete_pixels,
    clusters
  )
  
  # Preserve input names
  names(classified_images) <- names(input_images)
  
  # Plot classified images
  if (isTRUE(do_plot)) {
    
    for (classified_image in classified_images) {
      
      terra::plot(
        classified_image,
        col = class_colors,
        zlim = c(1, num_clusters),
        axes = FALSE
      )
    }
  }
  
  # Build output
  output <- list(
    result = classified_images,
    centers = kmeans_result$centers
  )
  
  return(output)
}

#### Test function ####

fire_class <- im.classify.timeseries(
  input_images = fire_list,
  num_clusters = 4,
  seed = 42,
  do_plot = TRUE
)

#### Plot parallel sets ####

# Create stack containing the three classifications

class_stack <- do.call(
  c,
  fire_class$result
)

# Convert to data frame

class_df <- terra::as.data.frame(
  class_stack,
  na.rm = TRUE
)

names(class_df) <- names(class_stack)

# Count temporal trajectories

transition_counts <- class_df |>
  dplyr::count(
    pre_fires,
    fires,
    post_fires,
    name = "count"
  )

# Convert classes to factors

transition_counts <- transition_counts |>
  mutate(
    across(
      c(pre_fires, fires, post_fires),
      ~ factor(
        .x,
        levels = 1:4
      )
    )
  )

# Prepare data for parallel sets

parallel_data <- ggforce::gather_set_data(
  transition_counts,
  1:3
)

# Order temporal axes

parallel_data$x <- factor(
  parallel_data$x,
  levels = c(
    "pre_fires",
    "fires",
    "post_fires"
  ),
  labels = c(
    "6 Aug 2017",
    "26 Aug 2017",
    "15 Sep 2017"
  )
)

# Set class colors

color_palette <- viridisLite::viridis(4)

# Plot parallel sets

ggplot(
  parallel_data,
  aes(
    x,
    id = id,
    split = y,
    value = count
  )
) +
  ggforce::geom_parallel_sets(
    aes(fill = pre_fires),
    alpha = 0.7,
    axis.width = 0.13
  ) +
  ggforce::geom_parallel_sets_axes(
    axis.width = 0.1,
    fill = "grey80",
    color = "grey80"
  ) +
  ggforce::geom_parallel_sets_labels(
    color = "black",
    size = 3,
    angle = 90
  ) +
  scale_x_discrete(
    name = NULL,
    expand = c(0, 0.2)
  ) +
  scale_y_continuous(
    breaks = NULL,
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = color_palette,
    guide = "none"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(
      color = "black",
      size = 10
    ),
    plot.margin = margin(14, 1.5, 2, 1.5)
  )

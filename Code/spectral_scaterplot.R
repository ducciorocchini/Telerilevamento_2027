# Load required packages
library(terra)
library(ggplot2)
library(dplyr)
library(imageRy)
library(cowplot)


# Import images and assign band names
sent <- im.import("ntinel2_NDVI_2020-")

names(sent) <- c(
  "Blue",
  "Green",
  "Red",
  "NIR"
)


# Classify pixels into three clusters
set.seed(42)

classified <- imageRy::im.classify(
  sent,
  num_clusters = 3,
  seed = 42,
  do_plot = FALSE
)


# Extract blue, nir, and cluster values for each pixel
pixel_data <- data.frame(
  Blue    = as.numeric(values(sent[["Blue"]], mat = FALSE)),
  NIR     = as.numeric(values(sent[["NIR"]], mat = FALSE)),
  Cluster = as.numeric(values(classified, mat = FALSE))
)


# Remove pixels with missing values
pixel_data <- pixel_data |>
  filter(
    !is.na(Blue),
    !is.na(NIR),
    !is.na(Cluster)
  )


# Assign labels a, b, and c to the three clusters
cluster_sizes <- pixel_data |>
  count(Cluster, sort = TRUE)

cluster_labels <- cluster_sizes$Cluster

pixel_data <- pixel_data |>
  mutate(
    Cluster_Label = case_when(
      Cluster == cluster_labels[1] ~ "A",
      Cluster == cluster_labels[2] ~ "B",
      Cluster == cluster_labels[3] ~ "C"
    )
  )


# Randomly sample 300 pixels for the plot
set.seed(1)

plot_data <- pixel_data |>
  slice_sample(n = 300)


# Define colours and point shapes
cluster_colours <- c(
  A = "#011959",
  B = "#D68A00",
  C = "#A3311D"
)

cluster_shapes <- c(
  A = 16,
  B = 17,
  C = 15
)


# Create the main plot of blue and nir pixel values
main_plot <- ggplot(
  plot_data,
  aes(
    x = Blue,
    y = NIR
  )
) +
  geom_point(
    aes(
      colour = Cluster_Label,
      shape = Cluster_Label
    ),
    size = 1.8,
    alpha = 0.6
  ) +
  scale_colour_manual(
    values = cluster_colours
  ) +
  scale_shape_manual(
    values = cluster_shapes
  ) +
  labs(
    title = "Spectral Scatterplot of Classified Sentinel-2 Pixels",
    subtitle = "300 sampled pixels with marginal cluster distributions",
    x = "Blue reflectance (B02)",
    y = "Near-infrared reflectance (B08)"
  ) +
  theme_minimal(
    base_size = 13
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 11)
  )


# Calculate the distribution of blue values for each cluster
blue_dens <- pixel_data |>
  group_by(Cluster_Label) |>
  group_modify(~ {
    d <- density(.x$Blue, na.rm = TRUE)
    data.frame(Blue = d$x, density = d$y)
  }) |>
  ungroup()


# Calculate the distribution of nir values for each cluster
nir_dens <- pixel_data |>
  group_by(Cluster_Label) |>
  group_modify(~ {
    d <- density(.x$NIR, na.rm = TRUE)
    data.frame(NIR = d$x, density = d$y)
  }) |>
  ungroup()


# Find the peak of each blue value distribution
blue_max <- blue_dens |>
  group_by(Cluster_Label) |>
  filter(density == max(density)) |>
  ungroup()


# Prepare labels for the blue value density curves
blue_labels <- blue_max |>
  mutate(
    label = Cluster_Label,
    hjust = case_when(
      Cluster_Label == "A" ~ 0,
      Cluster_Label == "B" ~ 0.5,
      Cluster_Label == "C" ~ 0
    ),
    vjust = case_when(
      Cluster_Label == "A" ~ 1,
      Cluster_Label == "B" ~ 0,
      Cluster_Label == "C" ~ 1
    ),
    nudge_x = case_when(
      Cluster_Label == "A" ~ 0.003,
      Cluster_Label == "B" ~ 0,
      Cluster_Label == "C" ~ 0.003
    ),
    nudge_y = case_when(
      Cluster_Label == "A" ~ 3,
      Cluster_Label == "B" ~ 5,
      Cluster_Label == "C" ~ 3
    )
  )


# Create the blue density plot to place above the main plot
xdens <- axis_canvas(
  main_plot,
  axis = "x"
) +
  geom_ribbon(
    data = blue_dens,
    aes(
      x = Blue,
      ymin = 0,
      ymax = density,
      fill = Cluster_Label,
      group = Cluster_Label
    ),
    alpha = 0.35
  ) +
  geom_line(
    data = blue_dens,
    aes(
      x = Blue,
      y = density,
      colour = Cluster_Label,
      group = Cluster_Label
    ),
    linewidth = 0.5
  ) +
  geom_text(
    data = blue_labels,
    aes(
      x = Blue + nudge_x,
      y = density + nudge_y,
      label = label,
      hjust = hjust,
      vjust = vjust,
      colour = Cluster_Label
    ),
    size = 13 / .pt,
    fontface = "bold",
    inherit.aes = FALSE
  ) +
  scale_colour_manual(
    values = cluster_colours,
    guide = "none"
  ) +
  scale_fill_manual(
    values = cluster_colours,
    guide = "none"
  ) +
  scale_y_continuous(
    limits = c(0, max(blue_dens$density) * 1.10),
    expand = c(0, 0)
  ) +
  coord_cartesian(clip = "off")


# Create the nir density plot to place on the right
ydens <- axis_canvas(
  main_plot,
  axis = "y",
  coord_flip = TRUE
) +
  geom_ribbon(
    data = nir_dens,
    aes(
      x = NIR,
      ymin = 0,
      ymax = density,
      fill = Cluster_Label,
      group = Cluster_Label
    ),
    alpha = 0.35
  ) +
  geom_line(
    data = nir_dens,
    aes(
      x = NIR,
      y = density,
      colour = Cluster_Label,
      group = Cluster_Label
    ),
    linewidth = 0.5
  ) +
  scale_colour_manual(
    values = cluster_colours,
    guide = "none"
  ) +
  scale_fill_manual(
    values = cluster_colours,
    guide = "none"
  ) +
  scale_y_continuous(
    limits = c(0, max(nir_dens$density) * 1.05),
    expand = c(0, 0)
  ) +
  coord_flip()


# Add the blue density plot above the main plot
p1 <- insert_xaxis_grob(
  main_plot + theme(legend.position = "none"),
  xdens,
  grid::unit(3 * 14, "pt"),
  position = "top",
  clip = "off"
)


# Add the nir density plot to the right
p2 <- insert_yaxis_grob(
  p1,
  ydens,
  grid::unit(3 * 14, "pt"),
  position = "right",
  clip = "off"
)


# Combine and display the final plot
marginal_plot <- ggdraw(p2)

marginal_plot

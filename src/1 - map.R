rm(list = ls())

library(tidyverse); library(sf); library(terra); library(tidyterra); source("src/colores.R")

### Load map layers

raster_file <- "data/map/altitud.tif"
vector_file <- "data/map/IberoAtlantic.shp"

r <- terra::rast(raster_file)
v <- sf::st_read(vector_file)

terra::crs(r)
sf::st_crs(v)

### Reproject CRS and check

if (sf::st_crs(v)$wkt != terra::crs(r)) {
  v <- sf::st_transform(v, terra::crs(r))
}

plot(r)
plot(v$geometry, add = TRUE, col = NA, border = "red", lwd = 2)

## Crop raster

v_vect <- terra::vect(v)
r_crop <- terra::crop(r, v_vect)
r_mask <- terra::mask(r_crop, v_vect)

### Plot limits

bb <- sf::st_bbox(v)

### Load data

read.csv("data/header.csv", fileEncoding = "latin1") -> header

ggplot() +
  ggthemes::theme_tufte() +
  geom_spatraster(data = r_mask, show.legend = FALSE) +
  scale_fill_gradient(
    low = "grey90",
    high = "grey10",
    na.value = "transparent",
    name = "Altitude") +
  geom_sf(data = v, fill = NA, color = "black", linewidth = 0.5) +
  coord_sf(
    xlim = c(bb["xmin"], bb["xmax"]),
    ylim = c(bb["ymin"], bb["ymax"])) +
  ggnewscale::new_scale_fill() +
  geom_point(data = header, aes(longitude, latitude, fill = habitat), shape = 21, size = 2.4) +
  xlab("Longitude (º)") + ylab("Latitude (º)") +
  theme(text = element_text(family = "sans"),
        strip.background = element_blank(),
        legend.position = c(.90,.30), 
        legend.direction = "vertical",
        legend.title = element_blank(),
        legend.margin = margin(0, 0, 0, 0),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"),
        legend.text = element_text(size = 10, face = "italic"), 
        panel.background = element_rect(color = "black", fill = NULL),
        panel.spacing = unit(0.05, "lines"),
        strip.text = element_text(size = 7, face = "italic"),
        # strip.text = element_blank(),
        plot.title = element_text(size = 12),
        axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_fill_manual(name = "Dominant invasive",
                     values = colores) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) -> Fig2; Fig2

ggsave(Fig2, file = "results/figures/F2 - map.png", bg = "white", 
       path = NULL, scale = 1, width = 180, units = "mm", dpi = 600)


library(sf)
library(rnaturalearth)
library(ggplot2)
library(cowplot)

# Europe base map
europe <- ne_countries(
  scale = "medium",
  continent = "Europe",
  returnclass = "sf"
)

# Make sure CRS match
ecoregion_wgs84 <- st_transform(v, 4326)

# Inset map
inset <- ggplot() +
  geom_sf(data = europe, fill = "grey90", colour = "grey65", linewidth = 0.2) +
  geom_sf(data = ecoregion_wgs84, fill = "grey40", colour = "black", linewidth = 0.3) +
  coord_sf(
    xlim = c(-12, 20),
    ylim = c(35, 60),
    expand = FALSE
  ) +
  theme_void() +
  theme(
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.4),
    plot.background = element_rect(fill = "white", colour = NA)
  )

inset

# Main map + inset
p_inset <- ggdraw() +
  draw_plot(Fig2) +
  draw_plot(
    inset,
    x = 0.50,   # horizontal position
    y = 0.18,   # vertical position
    width = 0.40,
    height = 0.40
  )

p_inset

ggsave(p_inset, file = "results/figures/Fig1.png", bg = "white", 
       path = NULL, scale = 1, width = 180, units = "mm", dpi = 600)

library(tidyverse); library(vegan); library(indicspecies); library(vegclust)
library(twinspanR); library(cluster); library(permute) 

rm(list = ls())

read.csv("data/header.csv", fileEncoding = "latin1") -> header
read.csv("data/species.csv", fileEncoding = "latin1") -> species

### Import SIVIM relevés

read.csv("data/header-legacy.csv", fileEncoding = "latin1") -> legacyplots
read.csv("data/species-legacy.csv", fileEncoding = "latin1") -> legacyspecies

### Species frequency by habitat in legacy plots

legacyplots %>%
  group_by(habitat) %>%
  tally %>% 
  rename(N = n) -> n
 
legacyplots %>%
  select(plotID, habitat) %>%
  merge(legacyspecies) %>%
  group_by(taxon, habitat) %>%
  tally %>% 
  merge(n) %>%
  mutate(f = n/N) %>%
  select(-c(n,N)) %>%
  filter(habitat != "Noise") %>%
  mutate(eunis1 =   str_sub(habitat, 1, 1)) %>% 
  mutate(eunis2 =   str_sub(habitat, 1, 2)) %>%
  mutate(eunis1 = ifelse(eunis1 == "M", "MA", eunis1)) %>%
  mutate(eunis2 = ifelse(eunis2 == "MA", "MA2", eunis2)) %>%
  rename(eunis3 = habitat) %>%
  select(taxon, f, eunis1, eunis2, eunis3) -> f

### Network, EUNIS level 1

header %>%
  select(plotID, habitat) %>%
  merge(species) %>%
  select(habitat, taxon) %>%
  unique %>%
  merge(f)  %>%
  group_by(habitat, eunis1) %>%
  summarise(f = sum(f)) %>%
  spread(eunis1, f, fill = 0) %>%
  column_to_rownames("habitat") %>%
  as.matrix()-> web

# 2) Ordenar hábitats alfabéticamente (columnas)
web <- web[, sort(colnames(web)), drop = FALSE]

# (Opcional) ordenar invasoras alfabéticamente (filas)
web <- web[sort(rownames(web)), , drop = FALSE]

# 3) Colores por invasora (filas = lower)
inv_names <- rownames(web)
inv_cols <- grDevices::hcl.colors(length(inv_names), palette = "Dark 3")
names(inv_cols) <- inv_names

# 4) Plot bipartito:
# - lower_color: colores de invasoras (filas)
# - higher_color: color fijo para hábitats (columnas) para no competir
# - link_color="lower": enlaces coloreados por invasora
# - srt=90: rota etiquetas (si quieres)

# Abrir dispositivo PNG
png(
  filename = "bipartite_invasive_habitat.png",
  width = 2000,      # en píxeles (ajusta según calidad deseada)
  height = 1600,
  res = 300          # resolución (300 dpi = calidad publicación)
)

# Dibujar el gráfico
bipartite::plotweb(
  web,
  sorting = "normal",
  lower_color = inv_cols,
  higher_color = "grey85",
  link_color = "lower",
  link_alpha = 0.8,
  srt = 90
)

# Cerrar dispositivo (MUY IMPORTANTE)
dev.off()  

### Network, EUNIS level 2

header %>%
  select(plotID, habitat) %>%
  merge(species) %>%
  select(habitat, taxon) %>%
  unique %>%
  merge(f)  %>%
  group_by(habitat, eunis2) %>%
  summarise(f = sum(f)) %>%
  spread(eunis2, f, fill = 0) %>%
  column_to_rownames("habitat") %>%
  as.matrix()-> web

# 2) Ordenar hábitats alfabéticamente (columnas)
web <- web[, sort(colnames(web)), drop = FALSE]

# (Opcional) ordenar invasoras alfabéticamente (filas)
web <- web[sort(rownames(web)), , drop = FALSE]

# 3) Colores por invasora (filas = lower)
inv_names <- rownames(web)
inv_cols <- grDevices::hcl.colors(length(inv_names), palette = "Dark 3")
names(inv_cols) <- inv_names

# 4) Plot bipartito:
# - lower_color: colores de invasoras (filas)
# - higher_color: color fijo para hábitats (columnas) para no competir
# - link_color="lower": enlaces coloreados por invasora
# - srt=90: rota etiquetas (si quieres)

# Abrir dispositivo PNG
png(
  filename = "bipartite_invasive_habitat.png",
  width = 2000,      # en píxeles (ajusta según calidad deseada)
  height = 1600,
  res = 300          # resolución (300 dpi = calidad publicación)
)

# Dibujar el gráfico
bipartite::plotweb(
  web,
  sorting = "normal",
  lower_color = inv_cols,
  higher_color = "grey85",
  link_color = "lower",
  link_alpha = 0.8,
  srt = 90
)

# Cerrar dispositivo (MUY IMPORTANTE)
dev.off()  

### Network, EUNIS level 3

header %>%
  select(plotID, habitat) %>%
  merge(species) %>%
  select(habitat, taxon) %>%
  unique %>%
  merge(f)  %>%
  group_by(habitat, eunis3) %>%
  summarise(f = sum(f)) %>%
  spread(eunis3, f, fill = 0) %>%
  column_to_rownames("habitat") %>%
  as.matrix()-> web

# 2) Ordenar hábitats alfabéticamente (columnas)
web <- web[, sort(colnames(web)), drop = FALSE]

# (Opcional) ordenar invasoras alfabéticamente (filas)
web <- web[sort(rownames(web)), , drop = FALSE]

# 3) Colores por invasora (filas = lower)
inv_names <- rownames(web)
inv_cols <- grDevices::hcl.colors(length(inv_names), palette = "Dark 3")
names(inv_cols) <- inv_names

# 4) Plot bipartito:
# - lower_color: colores de invasoras (filas)
# - higher_color: color fijo para hábitats (columnas) para no competir
# - link_color="lower": enlaces coloreados por invasora
# - srt=90: rota etiquetas (si quieres)

# Abrir dispositivo PNG
png(
  filename = "bipartite_invasive_habitat.png",
  width = 2000,      # en píxeles (ajusta según calidad deseada)
  height = 1600,
  res = 300          # resolución (300 dpi = calidad publicación)
)

# Dibujar el gráfico
bipartite::plotweb(
  web,
  sorting = "normal",
  lower_color = inv_cols,
  higher_color = "grey85",
  link_color = "lower",
  link_alpha = 0.8,
  srt = 90
)

# Cerrar dispositivo (MUY IMPORTANTE)
dev.off()  


### CHORD

# install.packages(c("tidyverse", "circlize"))
library(tidyverse)
library(circlize)

# df1 debe tener columnas: invasive, habitat, n
header %>%
  select(plotID, habitat) %>%
  merge(species) %>%
  select(habitat, taxon) %>%
  unique %>%
  merge(f)  %>%
  group_by(habitat, eunis1) %>%
  summarise(f = sum(f)) %>%
  rename(invasive = habitat,
         habitat = eunis1,
         n = f)  -> df1

links <- df1 %>%
  mutate(
    invasive = as.character(invasive),
    habitat  = as.character(habitat),
    n        = as.numeric(n)
  ) %>%
  filter(!is.na(invasive), !is.na(habitat), !is.na(n), n > 0)

# Niveles / orden
inv_levels <- sort(unique(links$invasive))                # invasoras (puedes cambiar a tu orden preferido)
hab_levels <- sort(unique(links$habitat))                 # hábitats alfabético

# Convertir a matriz invasoras x hábitats (como antes)
header %>%
  select(plotID, habitat) %>%
  merge(species) %>%
  select(habitat, taxon) %>%
  unique %>%
  merge(f)  %>%
  group_by(habitat, eunis1) %>%
  summarise(f = sum(f)) %>%
  spread(eunis1, f, fill = 0) %>%
  column_to_rownames("habitat") %>%
  as.matrix()-> web

# Colores: invasoras con paleta; hábitats en gris
inv_cols <- grDevices::hcl.colors(length(inv_levels), palette = "Dark 3")
names(inv_cols) <- inv_levels

hab_cols <- rep("grey85", length(hab_levels))
names(hab_cols) <- hab_levels

grid_cols <- c(inv_cols, hab_cols)

# Orden final de sectores: invasoras (arriba) + hábitats alfabéticos
sector_order <- c(inv_levels, hab_levels)

# (Recomendado) limpiar cualquier circlize previo
circos.clear()

# Parámetros del círculo (ajusta gap para separar grupos)
circos.par(
  start.degree = 90,
  gap.after = c(rep(2, length(inv_levels) - 1), 8, rep(2, length(hab_levels) - 1), 8)
)

# Dibujar chord: las cuerdas se colorean por invasora (origen)
chordDiagram(
  x = web,
  order = sector_order,
  grid.col = grid_cols,
  transparency = 0.2,
  directional = 1,                 # dirección desde invasora -> hábitat
  direction.type = c("arrows"),
  diffHeight = 0.02,
  link.arr.type = "big.arrow",
  link.sort = TRUE,
  link.largest.ontop = TRUE
)

# Añadir etiquetas (nombres de sectores)
circos.trackPlotRegion(
  track.index = 1, bg.border = NA,
  panel.fun = function(x, y) {
    sector_name <- get.cell.meta.data("sector.index")
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")
    circos.text(
      x = mean(xlim), y = ylim[1] + .1,
      labels = sector_name,
      facing = "clockwise", niceFacing = TRUE,
      adj = c(0, 0.5), cex = 0.7
    )
  }
)

# IMPORTANTE: limpiar al final
circos.clear()

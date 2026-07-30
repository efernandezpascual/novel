rm(list = ls())

library(tidyverse); library(bipartite); source("src/colores.R")

# =========================================================
# 1. CONSTRUIR MATRIZ (como ya hacías antes)
# =========================================================

read.csv("data/header.csv", fileEncoding = "latin1") -> header
read.csv("data/species.csv", fileEncoding = "latin1") -> species
read.csv("data/header-legacy.csv", fileEncoding = "latin1") -> legacyplots
read.csv("data/species-legacy.csv", fileEncoding = "latin1") -> legacyspecies

# Número de relevés por hábitat
legacyplots %>%
  group_by(habitat) %>%
  tally() %>%
  rename(N = n) -> n

# Frecuencia especies por hábitat
legacyplots %>%
  select(plotID, habitat) %>%
  merge(legacyspecies) %>%
  group_by(taxon, habitat) %>%
  tally() %>%
  merge(n) %>%
  mutate(f = n / N) %>%
  select(-c(n, N)) %>%
  filter(habitat != "Noise") %>%
  mutate(eunis1 = str_sub(habitat, 1, 1)) %>%
  mutate(eunis1 = ifelse(eunis1 == "M", "MA", eunis1)) %>%
  rename(eunis3 = habitat) %>%
  select(taxon, f, eunis1) %>%
  mutate(eunis1 = ifelse(eunis1 == "J", "V", eunis1)) -> f

# =========================================================
# 2. NETWORK EUNIS LEVEL 1
# =========================================================

header %>%
  select(plotID, habitat) %>%
  merge(species) %>%
  select(habitat, taxon) %>%
  unique() %>%
  merge(f) %>%
  group_by(habitat, eunis1) %>%
  summarise(f = sum(f), .groups = "drop") %>%
  spread(eunis1, f, fill = 0) %>%
  column_to_rownames("habitat") %>%
  as.matrix() -> web

# =========================================================
# 1) ASEGURAR ESTRUCTURA CORRECTA
# filas = especies
# columnas = EUNIS
# =========================================================

web_plot <- web

# Si tus especies NO están en filas, descomenta:
# web_plot <- t(web_plot)

# Verificación opcional:
# rownames(web_plot)  # deben ser especies
# colnames(web_plot)  # deben ser EUNIS


# =========================================================
# 2) ELIMINAR FLUJOS PEQUEÑOS
# =========================================================

threshold <- 5
web_plot[web_plot < threshold] <- 0

web_plot <- web_plot[rowSums(web_plot) > 0, , drop = FALSE]
web_plot <- web_plot[, colSums(web_plot) > 0, drop = FALSE]


# =========================================================
# 3) ORDENAR EUNIS POR APORTE TOTAL
# =========================================================

web_plot <- web_plot[, order(colSums(web_plot), decreasing = TRUE), drop = FALSE]

# ordenar especies alfabéticamente (opcional)
web_plot <- web_plot[sort(rownames(web_plot)), , drop = FALSE]

colnames(web_plot) <- c("Forests", "Grasslands", "Anthropogenic", "Shrublands", "Wetlands")


# =========================================================
# 4) COLORES POR ESPECIE
# =========================================================

sp_names <- rownames(web_plot)

# Asegurar que el vector colores esté nombrado correctamente
if(is.null(names(colores))){
  names(colores) <- sp_names
} else {
  colores <- colores[sp_names]
}

# =========================================================
# 5) EXPORTAR FIGURA
# =========================================================

png(
  filename = "results/figures/Fig3.png",
  width = 2200,
  height = 1600,
  res = 300
)

bipartite::plotweb(
  web_plot,
  horizontal = TRUE,
  sorting = "normal",
  
  # IZQUIERDA = EUNIS (gris)
  higher_color = "grey80",
  higher_border = "grey80",
  
  # DERECHA = especies (coloreadas)
  lower_color = colores,
  
  # enlaces coloreados por especie
  link_color = "lower",
  link_alpha = .8,
  
  # cursiva solo en especies
  lower_italic = TRUE,
  higher_italic = FALSE,
  
  spacing = 0.34,
  srt = 0
)

dev.off()

####### Now, Figure EUNIS level 2

# =========================================================
# 1. CONSTRUIR MATRIZ (como ya hacías antes)
# =========================================================

read.csv("data/header.csv", fileEncoding = "latin1") -> header
read.csv("data/species.csv", fileEncoding = "latin1") -> species
read.csv("data/header-legacy.csv", fileEncoding = "latin1") -> legacyplots
read.csv("data/species-legacy.csv", fileEncoding = "latin1") -> legacyspecies

# Número de relevés por hábitat
legacyplots %>%
  group_by(habitat) %>%
  tally() %>%
  rename(N = n) -> n

# Frecuencia especies por hábitat
legacyplots %>%
  select(plotID, habitat) %>%
  merge(legacyspecies) %>%
  group_by(taxon, habitat) %>%
  tally() %>%
  merge(n) %>%
  mutate(f = n / N) %>%
  select(-c(n, N)) %>%
  filter(habitat != "Noise") %>%
  mutate(eunis1 = ifelse(str_sub(habitat, 1, 2) == "MA",
                         str_sub(habitat, 1, 3),
                         str_sub(habitat, 1, 2))) %>%
  select(taxon, f, eunis1) -> f

# =========================================================
# 2. NETWORK EUNIS LEVEL 1
# =========================================================

header %>%
  select(plotID, habitat) %>%
  merge(species) %>%
  select(habitat, taxon) %>%
  unique() %>%
  merge(f) %>%
  group_by(habitat, eunis1) %>%
  summarise(f = sum(f), .groups = "drop") %>%
  spread(eunis1, f, fill = 0) %>%
  column_to_rownames("habitat") %>%
  as.matrix() -> web

# =========================================================
# 1) ASEGURAR ESTRUCTURA CORRECTA
# filas = especies
# columnas = EUNIS
# =========================================================

web_plot <- web

# Si tus especies NO están en filas, descomenta:
# web_plot <- t(web_plot)

# Verificación opcional:
# rownames(web_plot)  # deben ser especies
# colnames(web_plot)  # deben ser EUNIS



# =========================================================
# 4) COLORES POR ESPECIE
# =========================================================

sp_names <- rownames(web_plot)

# Asegurar que el vector colores esté nombrado correctamente
if(is.null(names(colores))){
  names(colores) <- sp_names
} else {
  colores <- colores[sp_names]
}

# =========================================================
# 5) EXPORTAR FIGURA
# =========================================================

png(
  filename = "results/figures/FigS2.png",
  width = 2200,
  height = 1600,
  res = 300
)

bipartite::plotweb(
  web_plot,
  horizontal = TRUE,
  sorting = "normal",
  
  # IZQUIERDA = EUNIS (gris)
  higher_color = "grey80",
  higher_border = "grey80",
  
  # DERECHA = especies (coloreadas)
  lower_color = colores,
  
  # enlaces coloreados por especie
  link_color = "lower",
  link_alpha = .8,
  
  # cursiva solo en especies
  lower_italic = TRUE,
  higher_italic = FALSE,
  
  spacing = 0.34,
  srt = 0
)

dev.off()

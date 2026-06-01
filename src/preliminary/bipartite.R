read.csv("results/invasives-x-habitats.csv") -> df1

library(tidyverse)
library(bipartite)

# 1) Matriz invasoras (filas) x hábitats (columnas)
web <- df1 %>%
  mutate(
    invasive = as.character(invasive),
    habitat  = as.character(habitat),
    n        = as.numeric(n)
  ) %>%
  filter(!is.na(invasive), !is.na(habitat), !is.na(n), n > 0) %>%
  pivot_wider(names_from = habitat, values_from = n, values_fill = 0) %>%
  column_to_rownames("invasive") %>%
  as.matrix()

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



rm(list = ls())

library(tidyverse); library(vegan); source("src/colores.R")

read.csv("data/species.csv", fileEncoding = "latin1") -> species
read.csv("data/header.csv", fileEncoding = "latin1") -> header

### NMDS plant communities

header %>%
  select(habitat, plotID) %>%
  merge(species) %>%
  select(plotID, taxon, cover) %>%
  spread(taxon, cover, fill = 0) %>%
  column_to_rownames(var = "plotID") -> matrix 

matrix %>%
  vegan::metaMDS(trymax = 100000, k = 2) -> 
  nmds # Ordination output

### Site scores

vegan::scores(nmds, "sites") %>%
  data.frame() %>%
  rownames_to_column("plotID") %>%
  rename(Dim.1 = NMDS1, Dim.2 = NMDS2) %>%
  merge(header) %>%
  dplyr::select(plotID, habitat, Dim.1, Dim.2) -> sitescores

### Species scores

vegan::scores(nmds, "species") %>%
  data.frame() %>%
  rownames_to_column("Taxon") %>%
  rename(Dim.1 = NMDS1, Dim.2 = NMDS2) %>%
  mutate(Species = Taxon) %>%
  separate(Species, into = c("G", "S")) %>%
  mutate(Species = paste(substr(G, 1, 3), substr(S, 1, 3), sep = "")) %>%
  dplyr::select(Taxon, Species, Dim.1, Dim.2) -> sppscores

### Plot NMDS sites

ggplot(sitescores, aes(x = Dim.1, y = Dim.2)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(aes(fill = habitat), color = "black", size = 3.5, shape = 21) +
  ggrepel::geom_text_repel(data = sppscores, aes(label = Species), color = "grey33", fontface = "italic", size = 2, segment.color = "transparent") +
  ggthemes::theme_tufte() +
  guides(fill = guide_legend(nrow = 1)) + 
  theme(text = element_text(family = "sans"),
        strip.background = element_blank(),
        legend.position = "bottom", 
        #legend.direction = "vertical",
        legend.title = element_blank(),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"),
        legend.text = element_text(size = 9, face = "italic"), 
        panel.background = element_rect(color = "black", fill = NULL),
        strip.text = element_text(size = 9.5, hjust = 0, margin = margin(l = 0, b = 4)),
        #strip.text = element_blank(),
        plot.title = element_text(size = 10),
        axis.title = element_text(size = 10),
        #axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_x_continuous(name = "NMDS1", labels = scales::number_format(accuracy = 0.1),
                     breaks = seq(-2, 2, by = 0.5)) + 
  scale_y_continuous(name = "NMDS2", labels = scales::number_format(accuracy = 0.1),
                     breaks = seq(-2, 2, by = 0.5)) +
  scale_color_manual(values = colores) +
  scale_fill_manual(values = colores) -> NMDS1;NMDS1

ggsave(NMDS1, file = "results/figures/Fig2.png", bg = "white", 
       path = NULL, scale = 1, width = 180, height = 2*75, units = "mm", dpi = 600)

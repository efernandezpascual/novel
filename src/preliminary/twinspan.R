library(tidyverse); library(vegan); library(indicspecies)

rm(list = ls())

colores <- c("#f8fda2", # Bidens
             "#a75bb3", # Buddleia
             "#8ea2a2", # Cortaderia
             "#fed81d", # Oenothera
             "#8bac39", # Paspalum
             "#780300", # Reynoutria
             "#fffff2") # Robinia

colores <- c(
  "#E6C84F",  # Bidens – amarillo flor (no pastel)
  "#8F4E9C",  # Buddleja – violeta
  "#6C9E8E",  # Cortaderia – verde azulado
  "#F2A541",  # Oenothera – naranja flor
  "#5A8F2E",  # Paspalum – verde gramínea
  "#8B1E1E",  # Reynoutria – rojo oscuro vegetal
  "#C7D07A"   # Robinia – verde claro hoja
)

read.csv("data/species.csv", fileEncoding = "latin1") -> species
read.csv("data/header.csv", fileEncoding = "latin1") -> header

species %>%
  select(plotID, taxon, cover) %>%
  spread(taxon, cover, fill = 0) %>%
  column_to_rownames(var = "plotID") -> df1

### Twinspan 

twinspanR::twinspan(
  df1,
  modif = TRUE,
  cut.levels = c(0, 15, 25),
  min.group.size = 10,
  clusters = 5, #Optimal number chosen with n-groups script
  diss = "multi.sorensen",
  mean.median = "mean",
  show.output.on.console = FALSE,
  quiet = TRUE) %>%
  cut() %>%
  cbind(df1) %>%
  select(".") %>%
  rename("Cluster" = ".") %>%
  rownames_to_column(var = "plotID") -> kclusters

header %>% merge(kclusters, by = "plotID") -> header4

### DCA

decorana(df1) -> dca1

vegan::scores(dca1) %>%
  data.frame() %>%
  rownames_to_column("plotID") %>%
  merge(header4) %>%
  mutate(Cluster = as.factor(Cluster)) -> df2

cent <- aggregate(cbind(DCA1, DCA2) ~ Cluster, data = df2, FUN = mean)
segs <- merge(df2, setNames(cent, c("Cluster", "oDCA1", "oDCA2")), by = "Cluster", sort = FALSE)

df2 %>%
  ggplot(aes(x = DCA1, y = DCA2)) + 
  geom_segment(data = segs, mapping = aes(xend = oDCA1, yend = oDCA2, color = Cluster), show.legend = F) +
  geom_point(data = cent, shape = 21, size = 5, aes(fill = Cluster), show.legend = T)

### NMDS plant communities

header4 %>%
  select(Cluster, plotID) %>%
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
  merge(header4) %>%
  dplyr::select(plotID, Cluster, Dim.1, Dim.2) -> sitescores

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
  geom_point(aes(fill = as.factor(Cluster)), color = "black", size = 3.5, shape = 21) +
  ggrepel::geom_text_repel(data = sppscores, aes(label = Species), color = "grey33", fontface = "italic", size = 2, segment.color = "transparent") +
  ggthemes::theme_tufte() +
  theme(text = element_text(family = "sans"),
        strip.background = element_blank(),
        legend.position = "right", 
        #legend.direction = "vertical",
        legend.title = element_text(size = 10),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"),
        legend.text = element_text(size = 10), 
        panel.background = element_rect(color = "black", fill = NULL),
        strip.text = element_text(size = 9.5, hjust = 0, margin = margin(l = 0, b = 4)),
        #strip.text = element_blank(),
        plot.title = element_text(size = 10),
        axis.title = element_blank(),
        #axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        plot.margin = unit(c(0.25,0,0,0), "cm")) +
  scale_x_continuous(name = "NMDS1", labels = scales::number_format(accuracy = 0.1),
                     breaks = seq(-2, 2, by = 0.5)) + 
  scale_y_continuous(name = "NMDS2", labels = scales::number_format(accuracy = 0.1),
                     breaks = seq(-2, 2, by = 0.5)) -> NMDS2;NMDS2

ggsave(NMDS2, file = "results/figures/NMDS2.png", bg = "white", 
       path = NULL, scale = 1, width = 180, height = 2*75, units = "mm", dpi = 600)

### Habitat

header4 %>%
  group_by(Cluster, habitat) %>%
  tally %>%
  arrange(Cluster, -n) %>%
  data.frame

### Dominant

species %>%
  merge(header4, by = "plotID") %>%
  group_by(Cluster, taxon) %>%
  summarise(D = sum(cover)) %>%
  group_by(Cluster) %>%
  slice_max(D, n = 5) %>%
  data.frame

### Indicator spp

species %>%
  merge(header4, by = "plotID") %>%
  select(Cluster, plotID, taxon, cover) %>%
  spread(taxon, cover, fill = 0) -> plots

plots %>%
  pull(Cluster) -> groups

comm <- plots %>%
  select(-c(Cluster, plotID)) %>%
  as.data.frame() %>%
  mutate(across(everything(), ~ as.numeric(as.character(.)))) %>%  # fuerza numérico
  as.matrix()

comm[is.na(comm)] <- 0

grp <- factor(groups)

set.seed(1)
iv <- multipatt(comm, grp, func = "IndVal.g",
                control = how(nperm = 999))

summary(iv)

iv_sig <- data.frame(iv$sign) |>
  tibble::rownames_to_column("species") |>
  dplyr::filter(p.value <= 0.05) |>
  dplyr::arrange(p.value)

iv_sig
library(tidyverse); library(vegan); library(indicspecies); library(vegclust)

rm(list = ls())

read.csv("data/species.csv", fileEncoding = "latin1") -> species
read.csv("data/header.csv", fileEncoding = "latin1") -> header

### Import SIVIM relevés

read.csv("data/species-legacy.csv", fileEncoding = "latin1")  %>%
  filter(taxon %in% c("Bidens aurea",
                               "Buddleja davidii",
                               "Cortaderia selloana",
                               "Reynoutria japonica",
                               "Oenothera glazioviana",
                               "Robinia pseudoacacia",
                               "Paspalum dilatatum")) %>%
  filter(cover > 30)
  pull(plotID) %>%
  unique -> ids

read.csv("data/header-legacy.csv", fileEncoding = "latin1") %>%
  filter(plotID %in% ids) -> legacyplots
read.csv("data/species-legacy.csv", fileEncoding = "latin1") %>%
  filter(plotID %in% ids) -> legacyspecies

rbind(species, legacyspecies) -> allspecies
rbind(header, legacyplots) -> allheader

allspecies %>%
  select(plotID, taxon, cover) %>%
  spread(taxon, cover, fill = 0) %>%
  column_to_rownames(var = "plotID") -> df1

### Twinspan 

twinspanR::twinspan(
  df1,
  modif = TRUE,
  cut.levels = c(0, 15, 25),
  min.group.size = 10,
  clusters = 10, #Optimal number chosen with n-groups script
  diss = "multi.sorensen",
  mean.median = "mean",
  show.output.on.console = FALSE,
  quiet = TRUE) %>%
  cut() %>%
  cbind(df1) %>%
  select(".") %>%
  rename("Cluster" = ".") %>%
  rownames_to_column(var = "plotID") -> kclusters

allheader %>% merge(kclusters, by = "plotID") -> header4

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


### Habitat

header4 %>%
  group_by(Cluster, habitat) %>%
  tally %>%
  arrange(Cluster, -n) %>%
  data.frame

### Dominant

allspecies %>%
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
library(tidyverse); library(DarkDiv)

rm(list = ls())

read.csv("data/header.csv", fileEncoding = "latin1") -> header
read.csv("data/species.csv", fileEncoding = "latin1") -> species
read.csv("data/header-legacy.csv", fileEncoding = "latin1") -> legacyplots
read.csv("data/species-legacy.csv", fileEncoding = "latin1") -> legacyspecies

# species %>%
#   rbind(legacyspecies) %>%
#   spread(taxon,cover, fill = 0) %>%
#   column_to_rownames(var = "plotID") -> x
# 
# DarkDiv(x, r = x, method = "Hypergeometric", limit = "min",
#         const = 0.01, removeAbsent = T, wa = F, weights = NULL) -> dd
# 
# save(dd, file = "results/darkdiversity.RData")

load(file = "results/darkdiversity.RData")

### Habitat specific species pools

dd$Pool %>%
  data.frame %>%
  rownames_to_column(var = "plotID") %>%
  gather(taxon, p, -plotID) %>%
  merge(header) -> pools

### Present species in the invaded habitats

merge(header, species) %>%
  # select(habitat, taxon) %>%
  select(taxon) %>%
  unique %>%
  mutate(alpha = 1) %>%
  mutate(taxon = gsub(" ", ".", taxon)) %>%
  mutate(taxon = gsub("-", ".", taxon)) -> alphas

### Missing species per invaded habitat

pools %>%
  group_by(habitat, taxon) %>%
  summarise(p = sum(p))  %>%
  merge(alphas, all.x = TRUE) %>%
  mutate(alpha = ifelse(is.na(alpha), 0, alpha)) %>%
  group_by(habitat, alpha) %>%
  # filter(alpha == 0) %>%
  arrange(habitat, -p) %>%
  filter(p >= 14) %>%
  group_by(taxon) %>%
  spread(habitat, p, fill = 0) %>%
  mutate(taxon = gsub("\\.", " ", taxon)) %>%
  write.csv("results/species-pools.csv", row.names = FALSE, fileEncoding = "latin1")

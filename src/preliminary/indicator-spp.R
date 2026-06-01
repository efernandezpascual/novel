library(tidyverse); library(vegan); library(indicspecies)

rm(list = ls())

read.csv("data/species.csv", fileEncoding = "latin1") -> species
read.csv("data/header.csv", fileEncoding = "latin1") -> header

### Indicator spp

species %>%
  merge(header, by = "plotID") %>%
  select(habitat, plotID, taxon, cover) %>%
  spread(taxon, cover, fill = 0) -> plots

plots %>%
  pull(habitat) -> groups

comm <- plots %>%
  select(-c(habitat, plotID)) %>%
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

### Constant

header %>%
  group_by(habitat) %>%
  tally() -> ns

species %>%
  merge(header, by = "plotID") %>%
  group_by(habitat, taxon) %>%
  summarise(F = length(taxon)) %>%
  merge(ns) %>%
  mutate(Frequency = F/n) %>%
  filter(Frequency >= .50) %>%
  mutate(Species = taxon,
         Type = "Constant",
         Value = Frequency*100) %>%
  select(Species, habitat, Type, Value) %>%
  arrange(habitat, -Value)

### Dominant

species %>%
  merge(header, by = "plotID") %>%
  filter(cover > 25) %>%
  group_by(habitat, taxon) %>%
  summarise(F = length(taxon)) %>%
  merge(ns) %>%
  mutate(Frequency = F/n) %>%
  filter(Frequency >= .05) %>%
  mutate(Species = taxon,
         Type = "Dominant",
         Value = Frequency*100) %>%
  select(Species, habitat, Type, Value) %>%
  arrange(habitat, -Value) 

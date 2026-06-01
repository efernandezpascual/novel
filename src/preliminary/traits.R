library(tidyverse)

rm(list = ls())

read.csv("data/sivim-eivs.csv", fileEncoding = "latin1") %>%
  mutate(taxon = gsub("-", " ", taxon)) -> traits
read.csv("results/species-pools.csv", fileEncoding = "latin1") %>%
  mutate(taxon = gsub("aggr ", "aggr.", taxon)) %>%
  mutate(taxon = gsub("sect ", "sect.", taxon)) %>%
  mutate(taxon = gsub("Salix x fragilis f  fragilis", "Salix x fragilis f. fragilis", taxon)) -> pools

merge(pools, traits, all.x = TRUE) %>%
  select(taxon, alpha, Disturbance.Frequency:Soil.Disturbance) %>%
  na.omit %>%
  gather(Trait, Value, Disturbance.Frequency:Soil.Disturbance) %>%
  ggplot(aes(alpha, Value, color = as.factor(alpha))) +
  geom_boxplot() +
  facet_wrap(~Trait, scales = "free_y")

  
  
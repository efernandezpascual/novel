### Con invasoras

rm(list = ls())

library(tidyverse)
library(vegan)
library(pairwiseAdonis)
set.seed(1)

species <- read.csv("data/species.csv", fileEncoding = "latin1")
header  <- read.csv("data/header.csv", fileEncoding = "latin1")

comm <- species %>%
  select(plotID, taxon, cover) %>%
  pivot_wider(names_from = taxon, values_from = cover, values_fill = 0) %>%
  as.data.frame()

rownames(comm) <- comm$plotID
comm$plotID <- NULL

meta <- header[match(rownames(comm), header$plotID), ]
meta$habitat <- factor(meta$habitat)

dist_mat <- vegdist(comm, method = "bray")

perma <- adonis2(dist_mat ~ habitat, data = meta, permutations = 9999)
perma

pairperma <- pairwise.adonis2(comm ~ habitat, data = meta, method = "bray",
                              permutations = 9999, nperm = 9999)
pairperma

tibble(par = names(pairperma)[-1],
       R2  = sapply(pairperma[-1], function(x) x$R2[1]),
       F   = sapply(pairperma[-1], function(x) x$F[1]),
       p   = sapply(pairperma[-1], function(x) x[["Pr(>F)"]][1])) %>%
  mutate(p_holm = p.adjust(p, "holm")) %>%
  arrange(p) %>%
  print(n = Inf)

betadis <- betadisper(dist_mat, meta$habitat)
anova(betadis)
permutest(betadis, permutations = 9999)

### Sin invasoras

rm(list = ls())

library(tidyverse)
library(vegan)
library(pairwiseAdonis)
set.seed(1)

species <- read.csv("data/species.csv", fileEncoding = "latin1") %>%
  filter(! taxon %in% c("Bidens aureus",
                      "Buddleja davidii",
                      "Cortaderia selloana",
                      "Reynoutria japonica",
                      "Oenothera glazioviana",
                      "Robinia pseudoacacia",
                      "Paspalum dilatatum"))
header  <- read.csv("data/header.csv", fileEncoding = "latin1")

comm <- species %>%
  select(plotID, taxon, cover) %>%
  pivot_wider(names_from = taxon, values_from = cover, values_fill = 0) %>%
  as.data.frame()

rownames(comm) <- comm$plotID
comm$plotID <- NULL

meta <- header[match(rownames(comm), header$plotID), ]
meta$habitat <- factor(meta$habitat)

dist_mat <- vegdist(comm, method = "bray")

perma <- adonis2(dist_mat ~ habitat, data = meta, permutations = 9999)
perma

pairperma <- pairwise.adonis2(comm ~ habitat, data = meta, method = "bray",
                              permutations = 9999, nperm = 9999)
pairperma

tibble(par = names(pairperma)[-1],
       R2  = sapply(pairperma[-1], function(x) x$R2[1]),
       F   = sapply(pairperma[-1], function(x) x$F[1]),
       p   = sapply(pairperma[-1], function(x) x[["Pr(>F)"]][1])) %>%
  mutate(p_holm = p.adjust(p, "holm")) %>%
  arrange(p) %>%
  print(n = Inf)

betadis <- betadisper(dist_mat, meta$habitat)
anova(betadis)
permutest(betadis, permutations = 9999)


library(tidyverse); library(DarkDiv)

rm(list = ls())

read.csv("data/header.csv", fileEncoding = "latin1") -> header
read.csv("data/species.csv", fileEncoding = "latin1") -> species
read.csv("data/header-legacy.csv", fileEncoding = "latin1") -> legacyplots
read.csv("data/species-legacy.csv", fileEncoding = "latin1") -> legacyspecies

legacyspecies %>% group_by(taxon) %>% tally %>% filter(n < 3) %>% pull(taxon) -> rarespp

legacyspecies %>% filter(! taxon %in% rarespp) -> legacyspecies

species %>%
  rbind(legacyspecies) %>%
  group_by(plotID, taxon) %>%
  filter(cover == max(cover)) %>%
  sample_n(1) %>%
  group_by() %>%
  spread(taxon,cover, fill = 0) %>%
  column_to_rownames(var = "plotID") -> x

DarkDiv(x, r = x, method = "Hypergeometric", limit = "min",
        const = 0.01, removeAbsent = T, wa = F, weights = NULL) -> dd

save(dd, file = "results/darkdiversity.RData")

load(file = "results/darkdiversity.RData")

### Plot specific species pools

dd$Pool %>%
  data.frame %>%
  rownames_to_column(var = "plotID") %>%
  gather(taxon, p, -plotID) %>%
  merge(header) %>%
  select(plotID:habitat) -> pools

### Observed species in the invaded habitats

merge(header, species) %>%
  select(habitat, taxon) %>%
  unique %>%
  mutate(observed = 1) %>%
  mutate(taxon = gsub(" ", ".", taxon)) %>%
  mutate(taxon = gsub("-", ".", taxon)) -> observed

### Native species per habitat

merge(header, species) %>%
  select(habitat, taxon) %>%
  unique %>%
  merge(read.csv("data/origin.csv", fileEncoding = "latin1")) %>%
  filter(origin == "native") %>%
  group_by(habitat) %>%
  tally %>%
  rename(threshold_n = n) -> threshold

### Habitat-specific species pools

pools %>%
  group_by(habitat, taxon) %>%
  summarise(p = sum(p))  %>%
  merge(observed, all.x = TRUE) %>%
  mutate(observed = ifelse(is.na(observed), 0, observed)) %>%
  mutate(taxon = gsub("\\.", " ", taxon))   %>%
  mutate(taxon = gsub("aggr ", "aggr.", taxon)) %>%
  mutate(taxon = gsub("sect ", "sect.", taxon)) %>%
  mutate(taxon = gsub("Salix x fragilis f  fragilis", "Salix x fragilis f. fragilis", taxon)) -> habitatpools

### Threshold (t = observed native species per habitat)

read.csv("data/origin.csv", fileEncoding = "latin1") %>%
  mutate(taxon = gsub("-", " ", taxon)) -> origin

habitatpools %>%
  filter(observed == 0) %>%
  merge(origin, all.x = TRUE) %>% 
  filter(origin == "native") %>%
  merge(threshold) %>%
  arrange(habitat, -p) %>%
  group_by(habitat)  %>%
  mutate(rank = row_number()) %>% 
  filter(rank <= threshold_n) %>%
  select(habitat, taxon, p, observed) -> darkdiv

merge(header, species) %>%
  select(habitat, taxon) %>%
  unique %>%
  merge(read.csv("data/origin.csv", fileEncoding = "latin1")) %>%
  filter(origin == "native") %>%
  select(habitat, taxon) %>%
  mutate(p = 15, observed = "1") %>%
  rbind(darkdiv) -> finalpools

finalpools %>%
  group_by(habitat, observed) %>%
  tally() %>%
  spread(observed, n)

finalpools %>% filter(observed == "1") %>% pull(taxon) %>% unique -> obsspp; length(obsspp)
finalpools %>% filter(observed == "0") %>% pull(taxon) %>% unique %>% length
finalpools %>% filter(observed == "0") %>% filter(! taxon %in% obsspp) %>% pull(taxon) %>% unique %>% length

finalpools %>% filter(observed == "0") %>% filter(! taxon %in% obsspp) %>%
  select(habitat, taxon) %>%
  unique %>%
  group_by(taxon) %>%
  tally %>% filter(n == 2)

finalpools %>% filter(observed == "0") %>% filter(! taxon %in% obsspp) %>%
  select(habitat, taxon) %>%
  unique %>%
  group_by(taxon) %>%
  tally %>%
  rename(n_habitats = n) %>%
  mutate(taxon = gsub("Dryopteris filix mas", "Dryopteris filix-mas", taxon)) %>%
  mutate(taxon = gsub("Silene flos cuculi", "Silene flos-cuculi", taxon)) %>%
  merge(legacyspecies, all.x = TRUE) %>%
  group_by(taxon, n_habitats) %>% 
  tally() %>%
  arrange(-n)

finalpools %>%
  write.csv("results/darkdiversity/species-pools-native.csv", row.names = FALSE, fileEncoding = "latin1")

read.csv("results/darkdiversity/species-pools-native.csv", fileEncoding = "latin1") %>%
  select(-p) %>%
  spread(habitat, observed, fill = "-") %>%
  write.csv("results/supplementary/Supplementary-Information-S3-v2.csv", fileEncoding = "latin1", row.names = FALSE)

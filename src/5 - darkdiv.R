library(tidyverse); library(DarkDiv)

rm(list = ls())

read.csv("data/header.csv", fileEncoding = "latin1") -> header
read.csv("data/species.csv", fileEncoding = "latin1") -> species
read.csv("data/species-legacy.csv", fileEncoding = "latin1") -> legacyspecies
read.csv("data/origin.csv", fileEncoding = "latin1") -> origin

legacyspecies %>% group_by(taxon) %>% tally %>% filter(n < 3) %>% pull(taxon) -> rarespp

legacyspecies %>% filter(! taxon %in% rarespp) -> legacyspecies

species %>%
  rbind(legacyspecies) %>%
  group_by(plotID, taxon) %>%
  filter(cover == max(cover)) %>%
  slice(1) %>%
  group_by() %>%
  spread(taxon,cover, fill = 0) %>%
  column_to_rownames(var = "plotID") -> x

# DarkDiv, hypergeometric method. 'limit' and 'const' only apply to the
# ThresholdBeals method and do not affect these results.
# The result is cached because the computation is slow.

dd_file <- "results/darkdiversity/darkdiversity.RData"
if (!file.exists(dd_file)) {
  DarkDiv(x, r = x, method = "Hypergeometric", limit = "min",
          const = 0.01, removeAbsent = TRUE, wa = FALSE, weights = NULL) -> dd
  save(dd, file = dd_file)
}
load(dd_file)

### Plot specific species pools

data.frame(dd$Pool, check.names = FALSE) %>%
  rownames_to_column(var = "plotID") %>%
  gather(taxon, p, -plotID) %>%
  merge(header) %>%
  select(plotID:habitat) -> pools

### Observed species in the invaded habitats

merge(header, species) %>%
  select(habitat, taxon) %>%
  unique %>%
  mutate(observed = 1) -> observed

### Native species per habitat

merge(header, species) %>%
  select(habitat, taxon) %>%
  unique %>%
  merge(origin) %>%
  filter(origin == "native") %>%
  group_by(habitat) %>%
  tally %>%
  rename(threshold_n = n) -> threshold

### Habitat-specific species pools

pools %>%
  group_by(habitat, taxon) %>%
  summarise(p = sum(p))  %>%
  merge(observed, all.x = TRUE) %>%
  mutate(observed = ifelse(is.na(observed), 0, observed)) -> habitatpools

### Criterion 1: top-N (as many candidates as realized species, highest probability first)

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
  merge(origin) %>%
  filter(origin == "native") %>%
  select(habitat, taxon) %>%
  mutate(p = 15, observed = 1) %>%
  rbind(darkdiv) -> pools_topn

write.csv(pools_topn, "results/darkdiversity/species-pools-topn.csv",
          row.names = FALSE, fileEncoding = "latin1")

pools_topn %>%
  select(-p) %>%
  spread(habitat, observed, fill = "-") %>%
  write.csv("results/supplementary/Supporting-Information-S3.csv", fileEncoding = "latin1", row.names = FALSE)

### Candidates for the probability criteria: native species absent from the type
### with mean membership probability > 0.5

realized <- merge(header, species) %>%
  distinct(habitat, taxon) %>%
  inner_join(origin, by = "taxon") %>%
  filter(origin == "native") %>%
  transmute(habitat, taxon, p = 15, observed = 1)

unrealized <- habitatpools %>%
  filter(observed == 0, p / 15 > 0.5) %>%
  inner_join(origin, by = "taxon") %>%
  filter(origin == "native") %>%
  transmute(habitat, taxon, p, observed = 0)

### Criterion 2: candidates with mean membership probability > 0.95

pools_p95 <- bind_rows(realized, unrealized %>% filter(p / 15 > 0.95))

write.csv(pools_p95, "results/darkdiversity/species-pools-p95.csv",
          row.names = FALSE, fileEncoding = "latin1")

### Criterion 3: candidates with mean membership probability > 0.5

pools_p05 <- bind_rows(realized, unrealized)

write.csv(pools_p05, "results/darkdiversity/species-pools-p05.csv",
          row.names = FALSE, fileEncoding = "latin1")
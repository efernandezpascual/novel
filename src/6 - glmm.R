rm(list = ls())

library(tidyverse); library(coda); source("src/colores.R")

### Open datasets

read.csv("data/sivim-eivs.csv", fileEncoding = "latin1") -> traits

list(topn = read.csv("results/darkdiversity/species-pools-topn.csv", fileEncoding = "latin1"),
     p95  = read.csv("results/darkdiversity/species-pools-p95.csv", fileEncoding = "latin1"),
     p05  = read.csv("results/darkdiversity/species-pools-p05.csv", fileEncoding = "latin1")) -> pools

### Data preparation

# Tree names: aggregates take the species name; spaces and hyphens become underscores
prepare <- function(d) {
  d %>%
    merge(traits, all.x = TRUE) %>%
    select(any_of(c("taxon", "observed", "habitat")), Disturbance.Frequency:Soil.Disturbance) %>%
    na.omit() %>%
    mutate(animal = gsub(" aggr.", "", taxon),
           animal = gsub("[ -]", "_", animal))
}

# Overall model: one row per species (realized if present in any type;
# unrealized if candidate in some type and present in none)
one_row <- function(d) {
  d %>% group_by(taxon) %>% summarise(observed = max(observed), .groups = "drop")
}

### Correlation among indicator values (variable selection)

prepare(pools$topn) -> df1

df1 %>% select(Disturbance.Frequency:Soil.Disturbance) %>% unique %>% cor
df1 %>% select(Mowing.Frequency, Soil.Disturbance, Grazing.Pressure, EIVE1.M, EIVE1.N) %>% unique %>% cor
df1 %>% select(Disturbance.Frequency:Soil.Disturbance) %>% unique %>% FactoMineR::PCA()
df1 %>% select(Mowing.Frequency, Soil.Disturbance, Grazing.Pressure, EIVE1.M, EIVE1.N) %>% unique %>% FactoMineR::PCA()

### Make tree (union of species in all pools; cached)

tree_file <- "results/tree/phylogeny.tre"
if (!file.exists(tree_file)) {
  bind_rows(pools) %>%
    distinct(taxon) %>%
    prepare() %>%
    distinct(animal) %>%
    separate(animal, into = c("genus", "epithet"), sep = "_", remove = FALSE, extra = "drop") %>%
    transmute(species = gsub("_", " ", animal), genus) %>%
    arrange(species) -> ranks1
  read.csv("data/plant_genus_list.csv") -> gen.list
  ape::read.tree("data/plant_megatree.tre") -> megatree
  U.PhyloMaker::phylo.maker(ranks1, megatree, gen.list, nodes.type = 1, scenario = 3) -> tree
  ape::write.tree(tree$phylo, tree_file)
  write.csv(tree$sp.list, "results/tree/phylogeny_splist.csv", row.names = FALSE, fileEncoding = "latin1")
}

### Read tree

phangorn::nnls.tree(cophenetic(ape::read.tree(tree_file)),
                    ape::read.tree(tree_file), method = "ultrametric") ->
  nnls_orig

nnls_orig$node.label <- NULL

### Set number of iterations

nite = 500000
nbur = 50000
nthi = 50

### Priors

list(R = list(V = 1, fix = 1), 
     G = list(G1 = list(V = diag(1), nu = 2.02))) ->
  priors # Priors for one random variable

### Fit models (cached: a model is only fitted if its file does not exist)

fit <- function(d, file) {
  if (file.exists(file)) return(invisible(NULL))
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  set.seed(1)
  MCMCglmm::MCMCglmm(observed ~
                       scale(Soil.Disturbance) + scale(Grazing.Pressure) + scale(Mowing.Frequency) +
                       scale(EIVE1.M) + scale(EIVE1.N),
                     random = ~ animal,
                     family = "threshold", pedigree = nnls_orig, prior = priors, data = d,
                     nitt = nite, thin = nthi, burnin = nbur,
                     verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> m
  save(m, file = file)
}

for (cr in names(pools)) {
  d <- prepare(pools[[cr]]) %>% filter(animal %in% nnls_orig$tip.label)
  for (h in sort(unique(d$habitat))) {
    fit(d %>% filter(habitat == h), file.path("results/models", cr, paste0("m-", tolower(h), ".Rdata")))
  }
  fit(prepare(one_row(pools[[cr]])) %>% filter(animal %in% nnls_orig$tip.label),
      file.path("results/models", cr, "m-all.Rdata"))
}

### Load main models (top-N criterion)

load_model <- function(f) { e <- new.env(); load(f, envir = e); get(ls(e)[1], envir = e) }

types <- c(Overall = "all", Robinia = "robinia", Bidens = "bidens", Cortaderia = "cortaderia",
           Paspalum = "paspalum", Buddleja = "buddleja", Oenothera = "oenothera", Reynoutria = "reynoutria")

models <- map(types, ~ load_model(file.path("results/models/topn", paste0("m-", .x, ".Rdata"))))

### Model diagnostics

# walk(models, plot)

### Model summary

map(models, summary)

### Export model tables

iwalk(models, ~ write.csv(summary(.x)$solutions,
                          file.path("results/models/topn", paste0("summary-", types[.y], ".csv")),
                          fileEncoding = "latin1"))

### Phylo signal http://www.mpcm-evolution.com/practice/online-practical-material-chapter-11/chapter-11-1-simple-model-mcmcglmm

map_dfr(models, function(m) {
  lambda <- coda::mcmc(as.numeric(m$VCV[, "animal"] / (m$VCV[, "animal"] + m$VCV[, "units"])))
  hpd <- coda::HPDinterval(lambda)
  tibble(lambda = round(mean(lambda), 2),
         min = round(hpd[, 1], 2),
         max = round(hpd[, 2], 2))
}, .id = "Model")

### Random effects

map(models, ~ summary(.x)$Gcovariances)

### Figure for manuscript

imap_dfr(models, ~ summary(.x)$solutions %>%
           data.frame %>%
           rownames_to_column(var = "Effect") %>%
           mutate(Model = .y)) %>%
  group_by %>%
  mutate(Effect = fct_recode(Effect,
                             "Moisture" = "scale(EIVE1.M)",
                             "Nutrients" = "scale(EIVE1.N)",
                             "Grazing" = "scale(Grazing.Pressure)",
                             "Mowing" = "scale(Mowing.Frequency)",
                             "Soil disturbance" = "scale(Soil.Disturbance)"))%>% 
  mutate(Effect = fct_relevel(Effect,
                              "Moisture",
                              "Nutrients",
                              "Soil disturbance",
                              "Grazing",
                              "Mowing")) %>%
  mutate(
    Model = fct_relevel(as.factor(Model), c("Robinia",
                                            "Reynoutria",
                                            "Paspalum",
                                            "Oenothera",
                                            "Cortaderia",
                                            "Buddleja",
                                            "Bidens",
                                            "Overall"))) %>%
  filter(! Effect == "(Intercept)") %>%
  ggplot(aes(y = Model, x = post.mean,
             xmin = l.95..CI, xmax = u.95..CI,
             color = Model)) +
  facet_wrap(~ Effect, nrow = 1) +
  geom_point(size = 2) +
  labs(x = "Effect size on the probability of a species being realized") +
  geom_errorbarh(height = .3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c(rev(colores), "grey10")) +
  ggthemes::theme_tufte() +
  theme(text = element_text(family = "sans"),
        strip.background = element_blank(),
        legend.position = "none", 
        #legend.direction = "vertical",
        legend.title = element_text(size = 10),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"),
        legend.text = element_text(size = 10), 
        panel.background = element_rect(color = "black", fill = NULL),
        strip.text = element_text(size = 9.8, margin = margin(l = 0, b = 4)),
        #strip.text = element_blank(),
        plot.title = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.text.y = element_text(color = "black", face = "italic"),
        axis.text.x = element_text(size = 7, color = "black"),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")) ->
  fig4; fig4

ggsave(fig4, file = "results/figures/Fig4.png", bg = "white", 
       path = NULL, scale = 1, width = 180, height = 80, units = "mm", dpi = 600)

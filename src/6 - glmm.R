library(U.PhyloMaker); library(tidyverse)

rm(list = ls())

### Open datasets

read.csv("data/sivim-eivs.csv", fileEncoding = "latin1") %>%
  mutate(taxon = gsub("-", " ", taxon)) -> traits
read.csv("results/darkdiversity/species-pools-native.csv", fileEncoding = "latin1") -> pools

merge(pools, traits, all.x = TRUE) %>%
  select(taxon, observed, habitat, Disturbance.Frequency:Soil.Disturbance) %>%
  na.omit %>%
  mutate(animal = gsub(" ", "_", taxon)) -> df1

### Make tree

# df1 %>%
#   select(taxon) %>%
#   unique %>%
#   separate(taxon, into = c("genus", "epithet"), sep = " ", remove = FALSE) %>%
#   select(taxon, genus) %>%
#   unique %>%
#   arrange(taxon) ->
#   ranks1
# 
# read.csv("data/plant_genus_list.csv") -> gen.list
# read.tree("data/plant_megatree.tre") -> megatree
# 
# phylo.maker(ranks1,
#             megatree,
#             gen.list,
#             nodes.type = 1,
#             scenario = 3) -> tree
# 
# write.tree(tree$phylo, "results/tree/disturbance.tre")
# write.csv(tree$sp.list, "results/tree/tree_splist.csv", row.names = FALSE, fileEncoding = "latin1")

### Read tree

phangorn::nnls.tree(cophenetic(ape::read.tree("results/tree/disturbance.tre")),
                    ape::read.tree("results/tree/disturbance.tre"), method = "ultrametric") ->
  nnls_orig

nnls_orig$node.label <- NULL

### Set number of iterations

nite = 100
nbur = 10
nthi = 5

nite = 500000
nbur = 50000
nthi = 50

### Priors

list(R = list(V = 1, fix = 1), 
     G = list(G1 = list(V = diag(1), nu = 2.02))) ->
  priors # Priors for one random variable

### Model for Robinia

df1 %>%
  filter(animal %in% nnls_orig$tip.label) %>%
  filter(habitat == "Robinia") -> df2

MCMCglmm::MCMCglmm(observed ~ scale(Disturbance.Frequency) + scale(Disturbance.Severity) +
                     scale(Soil.Disturbance) + scale(Grazing.Pressure) + scale(Mowing.Frequency) +
                     scale(EIVE1.T) + scale(EIVE1.M) + scale(EIVE1.L) + scale(EIVE1.N),
                   random = ~ animal,
                   family = "threshold", pedigree = nnls_orig, prior = priors, data = df2,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> m2

save(m2, file = "results/models/m-robinia.Rdata")

## Model for Bidens

df1 %>%
  filter(animal %in% nnls_orig$tip.label) %>%
  filter(habitat == "Bidens") -> df2

MCMCglmm::MCMCglmm(observed ~ scale(Disturbance.Frequency) + scale(Disturbance.Severity) +
                     scale(Soil.Disturbance) + scale(Grazing.Pressure) + scale(Mowing.Frequency) +
                     scale(EIVE1.T) + scale(EIVE1.M) + scale(EIVE1.L) + scale(EIVE1.N),
                   random = ~ animal,
                   family = "threshold", pedigree = nnls_orig, prior = priors, data = df2,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> m3

save(m3, file = "results/models/m-bidens.Rdata")

### Model for Cortaderia

df1 %>%
  filter(animal %in% nnls_orig$tip.label) %>%
  filter(habitat == "Cortaderia") -> df2

MCMCglmm::MCMCglmm(observed ~ scale(Disturbance.Frequency) + scale(Disturbance.Severity) +
                     scale(Soil.Disturbance) + scale(Grazing.Pressure) + scale(Mowing.Frequency) +
                     scale(EIVE1.T) + scale(EIVE1.M) + scale(EIVE1.L) + scale(EIVE1.N),
                   random = ~ animal,
                   family = "threshold", pedigree = nnls_orig, prior = priors, data = df2,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> m4

save(m4, file = "results/models/m-cortaderia.Rdata")

### Model for Paspalum

df1 %>%
  filter(animal %in% nnls_orig$tip.label) %>%
  filter(habitat == "Paspalum") -> df2

MCMCglmm::MCMCglmm(observed ~ scale(Disturbance.Frequency) + scale(Disturbance.Severity) +
                     scale(Soil.Disturbance) + scale(Grazing.Pressure) + scale(Mowing.Frequency) +
                     scale(EIVE1.T) + scale(EIVE1.M) + scale(EIVE1.L) + scale(EIVE1.N),
                   random = ~ animal,
                   family = "threshold", pedigree = nnls_orig, prior = priors, data = df2,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> m5

save(m5, file = "results/models/m-paspalum.Rdata")

### Model for Buddleja

df1 %>%
  filter(animal %in% nnls_orig$tip.label) %>%
  filter(habitat == "Buddleja") -> df2

MCMCglmm::MCMCglmm(observed ~ scale(Disturbance.Frequency) + scale(Disturbance.Severity) +
                     scale(Soil.Disturbance) + scale(Grazing.Pressure) + scale(Mowing.Frequency) +
                     scale(EIVE1.T) + scale(EIVE1.M) + scale(EIVE1.L) + scale(EIVE1.N),
                   random = ~ animal,
                   family = "threshold", pedigree = nnls_orig, prior = priors, data = df2,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> m6

save(m6, file = "results/models/m-buddleja.Rdata")

# ### Model for Oenothera

df1 %>%
  filter(animal %in% nnls_orig$tip.label) %>%
  filter(habitat == "Oenothera") -> df2

MCMCglmm::MCMCglmm(observed ~ scale(Disturbance.Frequency) + scale(Disturbance.Severity) +
                     scale(Soil.Disturbance) + scale(Grazing.Pressure) + scale(Mowing.Frequency) +
                     scale(EIVE1.T) + scale(EIVE1.M) + scale(EIVE1.L) + scale(EIVE1.N),
                   random = ~ animal,
                   family = "threshold", pedigree = nnls_orig, prior = priors, data = df2,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> m7

save(m7, file = "results/models/m-oenothera.Rdata")

### Model for Reynoutria

df1 %>%
  filter(animal %in% nnls_orig$tip.label) %>%
  filter(habitat == "Reynoutria") -> df2

MCMCglmm::MCMCglmm(observed ~ scale(Disturbance.Frequency) + scale(Disturbance.Severity) +
                     scale(Soil.Disturbance) + scale(Grazing.Pressure) + scale(Mowing.Frequency) +
                     scale(EIVE1.T) + scale(EIVE1.M) + scale(EIVE1.L) + scale(EIVE1.N),
                   random = ~ animal,
                   family = "threshold", pedigree = nnls_orig, prior = priors, data = df2,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> m8

save(m8, file = "results/models/m-reynoutria.Rdata")

### Model for all species

read.csv("data/sivim-eivs.csv", fileEncoding = "latin1") %>%
  mutate(taxon = gsub("-", " ", taxon)) -> traits
read.csv("results/species-pools-native.csv", fileEncoding = "latin1")  %>%
  select(taxon, observed) %>% unique  -> pools

merge(pools, traits, all.x = TRUE) %>%
  select(taxon, observed, Disturbance.Frequency:Soil.Disturbance) %>%
  na.omit %>%
  mutate(animal = gsub(" ", "_", taxon)) %>%
  filter(animal %in% nnls_orig$tip.label) -> df1

MCMCglmm::MCMCglmm(observed ~ scale(Disturbance.Frequency) + scale(Disturbance.Severity) +
                     scale(Soil.Disturbance) + scale(Grazing.Pressure) + scale(Mowing.Frequency) +
                     scale(EIVE1.T) + scale(EIVE1.M) + scale(EIVE1.L) + scale(EIVE1.N),
                   random = ~ animal,
                   family = "threshold", pedigree = nnls_orig, prior = priors, data = df1,
                   nitt = nite, thin = nthi, burnin = nbur,
                   verbose = FALSE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE, pr = FALSE, pl = FALSE) -> m1

save(m1, file = "results/models/m-all.Rdata")

### Load models

load(file = "results/models/m-all.Rdata")
load(file = "results/models/m-bidens.Rdata")
load(file = "results/models/m-buddleja.Rdata")
load(file = "results/models/m-cortaderia.Rdata")
load(file = "results/models/m-oenothera.Rdata")
load(file = "results/models/m-paspalum.Rdata")
load(file = "results/models/m-reynoutria.Rdata")
load(file = "results/models/m-robinia.Rdata")

### Model diagnostics

# plot(m1)
# plot(m2)
# plot(m3)
# plot(m4)
# plot(m5)
# plot(m6)
# plot(m7)
# plot(m8)

### Model summary

summary(m1) # All
summary(m2) # Robinia
summary(m3) # Bidens
summary(m4) # Cortaderia
summary(m5) # Paspalum
summary(m6) # Buddleja
summary(m7) # Oenothera
summary(m8) # Reynoutria

### Export model tables

write.csv(summary(m1)$solutions, "results/models/summary-m1.csv", row.names = FALSE, fileEncoding = "latin1")
write.csv(summary(m2)$solutions, "results/models/summary-m2.csv", row.names = FALSE, fileEncoding = "latin1")
write.csv(summary(m3)$solutions, "results/models/summary-m3.csv", row.names = FALSE, fileEncoding = "latin1")
write.csv(summary(m4)$solutions, "results/models/summary-m4.csv", row.names = FALSE, fileEncoding = "latin1")
write.csv(summary(m5)$solutions, "results/models/summary-m5.csv", row.names = FALSE, fileEncoding = "latin1")
write.csv(summary(m6)$solutions, "results/models/summary-m6.csv", row.names = FALSE, fileEncoding = "latin1")
write.csv(summary(m7)$solutions, "results/models/summary-m7.csv", row.names = FALSE, fileEncoding = "latin1")
write.csv(summary(m8)$solutions, "results/models/summary-m8.csv", row.names = FALSE, fileEncoding = "latin1")

### Phylo signal http://www.mpcm-evolution.com/practice/online-practical-material-chapter-11/chapter-11-1-simple-model-mcmcglmm

lambda <- m1$VCV[,"animal"]/(m1$VCV[,"animal"] + m1$VCV[,"units"])
data.frame(
  lambda = mean(lambda) %>% round(2),
  min = coda::HPDinterval(lambda)[, 1] %>% round(2),
  max = coda::HPDinterval(lambda)[, 2] %>% round(2))

lambda <- m2$VCV[,"animal"]/(m2$VCV[,"animal"] + m2$VCV[,"units"])
data.frame(
  lambda = mean(lambda) %>% round(2),
  min = coda::HPDinterval(lambda)[, 1] %>% round(2),
  max = coda::HPDinterval(lambda)[, 2] %>% round(2))

lambda <- m3$VCV[,"animal"]/(m3$VCV[,"animal"] + m3$VCV[,"units"])
data.frame(
  lambda = mean(lambda) %>% round(2),
  min = coda::HPDinterval(lambda)[, 1] %>% round(2),
  max = coda::HPDinterval(lambda)[, 2] %>% round(2))

lambda <- m4$VCV[,"animal"]/(m4$VCV[,"animal"] + m4$VCV[,"units"])
data.frame(
  lambda = mean(lambda) %>% round(2),
  min = coda::HPDinterval(lambda)[, 1] %>% round(2),
  max = coda::HPDinterval(lambda)[, 2] %>% round(2))

lambda <- m5$VCV[,"animal"]/(m5$VCV[,"animal"] + m5$VCV[,"units"])
data.frame(
  lambda = mean(lambda) %>% round(2),
  min = coda::HPDinterval(lambda)[, 1] %>% round(2),
  max = coda::HPDinterval(lambda)[, 2] %>% round(2))

lambda <- m6$VCV[,"animal"]/(m6$VCV[,"animal"] + m6$VCV[,"units"])
data.frame(
  lambda = mean(lambda) %>% round(2),
  min = coda::HPDinterval(lambda)[, 1] %>% round(2),
  max = coda::HPDinterval(lambda)[, 2] %>% round(2))

lambda <- m7$VCV[,"animal"]/(m7$VCV[,"animal"] + m7$VCV[,"units"])
data.frame(
  lambda = mean(lambda) %>% round(2),
  min = coda::HPDinterval(lambda)[, 1] %>% round(2),
  max = coda::HPDinterval(lambda)[, 2] %>% round(2))

lambda <- m8$VCV[,"animal"]/(m8$VCV[,"animal"] + m8$VCV[,"units"])
data.frame(
  lambda = mean(lambda) %>% round(2),
  min = coda::HPDinterval(lambda)[, 1] %>% round(2),
  max = coda::HPDinterval(lambda)[, 2] %>% round(2))

### Random effects

summary(m1)$Gcovariances
summary(m2)$Gcovariances
summary(m3)$Gcovariances
summary(m4)$Gcovariances
summary(m5)$Gcovariances
summary(m6)$Gcovariances
summary(m7)$Gcovariances
summary(m8)$Gcovariances

### Figure for manuscript

rbind(
  summary(m1)$solutions %>%
    data.frame %>%
    rownames_to_column(var = "Effect") %>%
    mutate(Model = "Overall"),
  summary(m2)$solutions %>%
    data.frame %>%
    rownames_to_column(var = "Effect") %>%
    mutate(Model = "Robinia"),
  summary(m3)$solutions %>%
    data.frame %>%
    rownames_to_column(var = "Effect") %>%
    mutate(Model = "Bidens"),
  summary(m4)$solutions %>%
    data.frame %>%
    rownames_to_column(var = "Effect") %>%
    mutate(Model = "Cortaderia"),
  summary(m5)$solutions %>%
    data.frame %>%
    rownames_to_column(var = "Effect") %>%
    mutate(Model = "Paspalum"),
  summary(m6)$solutions %>%
    data.frame %>%
    rownames_to_column(var = "Effect") %>%
    mutate(Model = "Buddleja"),
  summary(m7)$solutions %>%
    data.frame %>%
    rownames_to_column(var = "Effect") %>%
    mutate(Model = "Oenothera"),
  summary(m8)$solutions %>%
    data.frame %>%
    rownames_to_column(var = "Effect") %>%
    mutate(Model = "Reynoutria")) %>%
  group_by %>%
  mutate(Effect = fct_recode(Effect,
                             "Disturbance frequency" = "scale(Disturbance.Frequency)",
                             "Disturbance severity" = "scale(Disturbance.Severity)",
                             "Light" = "scale(EIVE1.L)",
                             "Moisture" = "scale(EIVE1.M)",
                             "Nutrients" = "scale(EIVE1.N)",
                             "Temperature" = "scale(EIVE1.T)",
                             "Grazing" = "scale(Grazing.Pressure)",
                             "Mowing" = "scale(Mowing.Frequency)",
                             "Soil disturbance" = "scale(Soil.Disturbance)"))%>% 
  mutate(Effect = fct_relevel(Effect,
                              "Temperature",
                              "Moisture",
                              "Light",
                              "Nutrients",
                              "Disturbance frequency",
                              "Disturbance severity",
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
  # mutate(
  #   Model = fct_recode(Model, "Wetlands - High disturbance" = "m2.8",
  #                      "Wetlands - Low disturbance" = "m2.7",
  #                      "Cold stress - High disturbance" = "m2.6",
  #                      "Cold stress - Low disturbance" = "m2.5",
  #                      "Drought stress - High disturbance" = "m2.4",
  #                      "Drought stress - Low disturbance" =  "m2.3",
  #                      "Low stress - High disturbance" = "m2.2",
  #                      "Low stress - Low disturbance" = "m2.1")) %>%
  filter(! Effect == "(Intercept)") %>%
  #filter(pMCMC < 0.05) %>%
  ggplot(aes(y = Model, x = post.mean,
             xmin = l.95..CI, xmax = u.95..CI,
             color = Model)) +
  facet_wrap(~ Effect, nrow = 3) +
  geom_point(size = 2) +
  labs(x = "Effect size") +
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
  fig6; fig6

ggsave(fig6, file = "results/figures/Fig4.png", bg = "white", 
       path = NULL, scale = 1, width = 180, height = 120, units = "mm", dpi = 600)
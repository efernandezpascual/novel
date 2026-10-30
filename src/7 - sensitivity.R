library(tidyverse)

rm(list = ls())

# =========================================================
# Sensitivity of the realized vs unrealized comparison to the
# definition of the unrealized component (Supporting Information S5).
# Criteria (nested cuts of the same probability ranking):
#   topn : as many candidates as realized species, highest probability first
#   p95  : candidates with mean membership probability > 0.95
#   p05  : candidates with mean membership probability > 0.5
# Requires the outputs of 5 - darkdiv.R and 6 - glmm.R. No models are fitted here.
# =========================================================

### Data

read.csv("data/sivim-eivs.csv", fileEncoding = "latin1") -> traits

list(topn = read.csv("results/darkdiversity/species-pools-topn.csv", fileEncoding = "latin1"),
     p95  = read.csv("results/darkdiversity/species-pools-p95.csv", fileEncoding = "latin1"),
     p05  = read.csv("results/darkdiversity/species-pools-p05.csv", fileEncoding = "latin1")) -> pools

criteria <- c(topn = "Highest-ranked", p95 = "P > 0.95", p05 = "P > 0.5")

### 1. Pool sizes per criterion

imap_dfr(pools, ~ .x %>% count(habitat, observed) %>% mutate(criterion = .y)) %>%
  mutate(component = ifelse(observed == 1, "realized", "unrealized")) %>%
  select(-observed) %>%
  pivot_wider(names_from = c(criterion, component), values_from = n) %>%
  select(habitat, realized = topn_realized, ends_with("unrealized")) -> sizes

# Overall model: one row per species (realized if present in any type)
imap_dfr(pools, ~ .x %>% group_by(taxon) %>% summarise(observed = max(observed), .groups = "drop") %>%
           count(observed) %>% mutate(criterion = .y)) %>%
  mutate(component = ifelse(observed == 1, "realized", "unrealized")) %>%
  select(-observed) %>%
  pivot_wider(names_from = c(criterion, component), values_from = n) %>%
  transmute(habitat = "Overall", realized = topn_realized,
            topn_unrealized, p95_unrealized, p05_unrealized) %>%
  bind_rows(sizes, .) -> sizes

sizes
write.csv(sizes, "results/supplementary/S5-pool-sizes.csv", row.names = FALSE, fileEncoding = "latin1")

### 2. Model coefficients per criterion

load_model <- function(f) { e <- new.env(); load(f, envir = e); get(ls(e)[1], envir = e) }

types <- c(Overall = "all", Bidens = "bidens", Buddleja = "buddleja", Cortaderia = "cortaderia",
           Oenothera = "oenothera", Paspalum = "paspalum", Reynoutria = "reynoutria", Robinia = "robinia")

effects <- c("scale(EIVE1.M)" = "Moisture", "scale(EIVE1.N)" = "Nutrients",
             "scale(Soil.Disturbance)" = "Soil disturbance", "scale(Grazing.Pressure)" = "Grazing",
             "scale(Mowing.Frequency)" = "Mowing")

map_dfr(names(criteria), function(cr) {
  imap_dfr(types, function(ty, model) {
    s <- summary(load_model(file.path("results/models", cr, paste0("m-", ty, ".Rdata"))))$solutions
    tibble(criterion = cr, Model = model, term = rownames(s),
           mean = s[, "post.mean"], lower = s[, "l-95% CI"], upper = s[, "u-95% CI"])
  })
}) %>%
  filter(term %in% names(effects)) %>%
  mutate(Effect = factor(effects[term], levels = effects),
         criterion = factor(criteria[criterion], levels = criteria),
         Model = factor(Model, levels = names(types))) -> coefs

# Table: posterior mean [95% CI], * when the interval excludes zero
coefs %>%
  mutate(cell = sprintf("%.2f [%.2f, %.2f]%s", mean, lower, upper,
                        ifelse(lower > 0 | upper < 0, "*", ""))) %>%
  select(Model, criterion, Effect, cell) %>%
  pivot_wider(names_from = Effect, values_from = cell) %>%
  arrange(Model, criterion) -> coef_table

print(coef_table, n = Inf)
write.csv(coef_table, "results/supplementary/S5-coefficients.csv", row.names = FALSE, fileEncoding = "latin1")

# Figure: effect sizes per criterion
coefs %>%
  mutate(Model = fct_rev(Model)) %>%
  ggplot(aes(y = Model, x = mean, xmin = lower, xmax = upper, color = criterion)) +
  facet_wrap(~ Effect, nrow = 1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(position = position_dodge(width = 0.6), size = 0.2) +
  scale_color_manual(values = c("grey10", "#2c7fb8", "#d95f0e"), name = "Unrealized component") +
  labs(x = "Effect size on the probability of a species being realized") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(face = "italic"),
        legend.position = "bottom") -> fig_coef

ggsave(fig_coef, file = "results/supplementary/FigS5-coefficients.png", bg = "white",
       width = 180, height = 110, units = "mm", dpi = 600)

### 3. Candidates along the probability ranking

# All candidates with probability > 0.5 (the broadest pool), in 0.05-wide bins
pools$p05 %>%
  mutate(prob = p / 15) %>%
  inner_join(traits %>% select(taxon, EIVE1.M, EIVE1.N, Soil.Disturbance, Grazing.Pressure, Mowing.Frequency),
             by = "taxon") -> d

# Lowest probability included by each criterion, per type
pools$topn %>%
  filter(observed == 0) %>%
  group_by(habitat) %>%
  summarise(cut = min(p / 15), .groups = "drop") %>%
  mutate(criterion = "Highest-ranked") %>%
  bind_rows(tibble(habitat = unique(d$habitat), cut = 0.95, criterion = "P > 0.95")) -> cuts

cuts %>% filter(criterion == "Highest-ranked")   # probability at which the Highest-ranked criterion stops

d %>%
  filter(observed == 0) %>%
  mutate(mid = 0.5 + 0.05 * as.integer(cut(prob, breaks = seq(0.5, 1, 0.05), include.lowest = TRUE)) - 0.025) -> cand

# (a) Number of candidates per probability bin
ggplot(cand %>% count(habitat, mid), aes(mid, n)) +
  geom_col(fill = "grey60") +
  geom_vline(data = cuts, aes(xintercept = cut, linetype = criterion)) +
  facet_wrap(~ habitat) +
  scale_x_reverse() +
  scale_linetype_manual(values = c("Highest-ranked" = "solid", "P > 0.95" = "dashed"), 
                        breaks = c("Highest-ranked", "P > 0.95"),
                        name = NULL) +
  labs(x = "Mean membership probability", y = "Number of native candidates") +
  theme_bw() +
  theme(legend.position = "bottom") -> fig_counts

ggsave(fig_counts, file = "results/supplementary/FigS5-counts.png", bg = "white",
       width = 180, height = 150, units = "mm", dpi = 600)

# (b) Mean standardized indicator values per bin vs. realized species
vars <- c(EIVE1.M = "Moisture", EIVE1.N = "Nutrients", Soil.Disturbance = "Soil disturbance",
          Grazing.Pressure = "Grazing", Mowing.Frequency = "Mowing")

d %>%
  mutate(across(all_of(names(vars)), ~ as.numeric(scale(.x)))) %>%
  pivot_longer(all_of(names(vars)), names_to = "var") %>%
  mutate(Effect = factor(vars[var], levels = vars)) -> z

z %>%
  filter(observed == 1) %>%
  group_by(habitat, Effect) %>%
  summarise(realized = mean(value), .groups = "drop") -> z_realized

z %>%
  filter(observed == 0) %>%
  mutate(mid = 0.5 + 0.05 * as.integer(cut(prob, breaks = seq(0.5, 1, 0.05), include.lowest = TRUE)) - 0.025) %>%
  group_by(habitat, Effect, mid) %>%
  summarise(mean = mean(value), n = n(), .groups = "drop") -> z_bins

ggplot(z_bins, aes(mid, mean)) +
  geom_hline(data = z_realized, aes(yintercept = realized), colour = "#d95f0e") +
  geom_vline(data = cuts, aes(xintercept = cut, linetype = criterion)) +
  geom_line() +
  geom_point(aes(size = n)) +
  facet_grid(Effect ~ habitat, scales = "free_y", labeller = labeller(Effect = label_wrap_gen(10))) +
  scale_x_reverse() +
  scale_size_area(max_size = 2.5, name = "Candidates") +
  scale_linetype_manual(values = c("Highest-ranked" = "solid", "P > 0.95" = "dashed"), 
                        breaks = c("Highest-ranked", "P > 0.95"),
                        name = NULL) +
  labs(x = "Mean membership probability",
       y = "Mean standardized indicator value of candidates") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 6)) -> fig_traits

ggsave(fig_traits, file = "results/supplementary/FigS5-traits.png", bg = "white",
       width = 180, height = 150, units = "mm", dpi = 600)

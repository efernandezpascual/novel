library(tidyverse)
library(twinspanR)
library(vegan)
library(cluster)
library(indicspecies)
library(permute)   # <-- CLAVE: aquí está how()

## ---- 0) Datos (df1 debe ser parcelas x especies, SOLO numérico) ----
if (!all(vapply(as.data.frame(df1), is.numeric, logical(1)))) {
  stop("df1 contiene columnas no numéricas. Elimina plotID/Cluster u otras columnas antes.")
}

comm <- as.matrix(df1)
mode(comm) <- "numeric"
comm[is.na(comm)] <- 0
if (is.null(rownames(comm))) rownames(comm) <- paste0("plot_", seq_len(nrow(comm)))

message("df1: ", nrow(comm), " parcelas x ", ncol(comm), " especies")

## ---- 1) Parámetros (tuyos) ----
tw_params <- list(
  modif = TRUE,
  cut.levels = c(0, 15, 25),
  min.group.size = 10,
  diss = "multi.sorensen",
  mean.median = "mean",
  show.output.on.console = FALSE,
  quiet = TRUE
)

run_tw <- function(comm_mat, k, params = tw_params) {
  twinspanR::twinspan(
    comm_mat,
    modif = params$modif,
    cut.levels = params$cut.levels,
    min.group.size = params$min.group.size,
    clusters = k,
    diss = params$diss,
    mean.median = params$mean.median,
    show.output.on.console = params$show.output.on.console,
    quiet = params$quiet
  )
}

## ---- 2) Distancia para evaluación ----
dist_eval <- vegan::vegdist(comm, method = "bray")

## ---- 3) Métricas ----
silhouette_mean <- function(groups, dist_obj) {
  sil <- cluster::silhouette(as.integer(groups), dist_obj)
  mean(sil[, "sil_width"])
}

within_between_ratio <- function(groups, dist_obj) {
  D <- as.matrix(dist_obj)
  g <- groups
  same <- outer(g, g, FUN = "==")
  diag(same) <- FALSE
  within <- mean(D[same], na.rm = TRUE)
  between <- mean(D[!same], na.rm = TRUE)
  within / between
}

indval_summary <- function(comm_mat, groups, nperm = 999, seed = 1) {
  stopifnot(nrow(comm_mat) == length(groups))
  set.seed(seed)
  
  iv <- indicspecies::multipatt(
    comm_mat, groups,
    func = "IndVal.g",
    control = permute::how(nperm = nperm)  # <-- AQUÍ la corrección
  )
  
  s <- as.data.frame(iv$sign) %>% tibble::rownames_to_column("species")
  sig <- s %>% dplyr::filter(!is.na(p.value), p.value <= 0.05)
  
  tibble(
    n_indic_sig = nrow(sig),
    indval_median_sig = if (nrow(sig) > 0) median(sig$stat, na.rm = TRUE) else NA_real_,
    indval_max_sig = if (nrow(sig) > 0) max(sig$stat, na.rm = TRUE) else NA_real_
  )
}

## ---- 4) Intentar k (definición realista de viable) ----
try_k <- function(comm_mat, k) {
  out <- try({
    tw <- run_tw(comm_mat, k, tw_params)
    grp <- factor(cut(tw))
    if (length(grp) != nrow(comm_mat)) stop("cut(tw) no devuelve longitud nrow")
    if (any(is.na(grp))) stop("grupos contienen NA")
    if (nlevels(grp) < 2) stop("solo 1 grupo")
    list(ok = TRUE, tw = tw, grp = grp, msg = NA_character_)
  }, silent = TRUE)
  
  if (inherits(out, "try-error")) {
    list(ok = FALSE, tw = NULL, grp = NULL, msg = as.character(out))
  } else out
}

## ---- 5) Rango de k a explorar ----
k_values <- 2:min(20, nrow(comm) - 1)

## ---- 6) Evaluar todos los k que funcionen ----
results <- purrr::map(k_values, function(k) {
  message("Probando k = ", k)
  res <- try_k(comm, k)
  
  if (!res$ok) {
    message("  -> k falló: ", res$msg)
    return(NULL)
  }
  
  grp <- res$grp
  sil <- silhouette_mean(grp, dist_eval)
  wbr <- within_between_ratio(grp, dist_eval)
  ivs <- indval_summary(comm, grp, nperm = 999, seed = 1)
  
  tibble(
    k = k,
    silhouette = sil,
    within_between_ratio = wbr,
    n_indic_sig = ivs$n_indic_sig,
    indval_median_sig = ivs$indval_median_sig,
    indval_max_sig = ivs$indval_max_sig,
    n_groups_real = nlevels(grp),
    min_group_size = min(table(grp)),
    tw_obj = list(res$tw),
    groups_obj = list(grp)
  )
}) %>%
  purrr::compact() %>%
  bind_rows()

if (nrow(results) == 0) stop("Ningún k funcionó. Revisa df1 y/o parámetros de twinspan.")

## ---- 7) Tabla resumen + gráficos ----
summary_tbl <- results %>%
  select(k, silhouette, within_between_ratio,
         n_indic_sig, indval_median_sig, indval_max_sig,
         n_groups_real, min_group_size) %>%
  arrange(desc(silhouette))

print(summary_tbl)

ggplot(summary_tbl, aes(k, silhouette)) + geom_line() + geom_point() +
  labs(title = "Silhouette (Bray–Curtis)", x = "k", y = "Silhouette")

ggplot(summary_tbl, aes(k, within_between_ratio)) + geom_line() + geom_point() +
  labs(title = "Within/Between ratio (menor = mejor)", x = "k", y = "Ratio")

ggplot(summary_tbl, aes(k, n_indic_sig)) + geom_line() + geom_point() +
  labs(title = "Nº especies indicadoras (IndVal, p<=0.05)", x = "k", y = "N indicadoras")

## ---- 8) Elegir k_best (heurística) ----
candidate <- summary_tbl %>%
  mutate(
    z_sil = as.numeric(scale(silhouette)),
    z_ind = as.numeric(scale(n_indic_sig)),
    z_rat = as.numeric(scale(within_between_ratio)),
    penalty = pmax(0, 3 - min_group_size),
    score = z_sil + z_ind - z_rat - 0.25 * as.numeric(scale(penalty))
  ) %>%
  arrange(desc(score))

print(candidate)

k_best <- candidate$k[1]
message("k recomendado (heurístico): ", k_best)

## ---- 9) Tu salida final tipo kclusters ----
best_row <- which(results$k == k_best)
best_groups <- results$groups_obj[[best_row]]

kclusters <- tibble(
  plotID = rownames(comm),
  Cluster = as.character(best_groups)
) %>%
  left_join(as_tibble(comm, rownames = "plotID"), by = "plotID")

print(head(kclusters))

## ---- 10) (Opcional) Resumen de indicadoras para k_best ----
set.seed(1)
iv_best <- indicspecies::multipatt(
  comm, best_groups,
  func = "IndVal.g",
  control = permute::how(nperm = 999)
)
summary(iv_best)

# ##############################################################################
# ##   PROVA DEI DUE NULLI SPAZIALI SU TUTTI I DATASET                        ##
# ##                                                                          ##
# ##   Lancia prova_nulli.R per ogni checkpoint, in processi R separati e in  ##
# ##   parallelo, poi mette insieme tabelle e grafici:                        ##
# ##     falsi centri (calibrazione dei test), istogrammi dei p-value,        ##
# ##     Golden Standard (boxplot, ROC, AUC con IC) e GSEA (AUC di contesto,  ##
# ##     NES, pathway con FDR < 0.05).                                        ##
# ##############################################################################

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })

SCRIPT_DIR    <- "/home/ammendola/laboratorio/spatial/RGdist"   # cartella con prova_nulli.R
CONFRONTO_DIR <- "/home/ammendola/laboratorio/spatial/test_arricchimenti/prova_nulli_confronto"
RIFAI    <- TRUE     # FALSE = non rilancia, rifa' solo tabelle e grafici
N_CPU    <- 40
MAX_JOBS <- 6        # dataset in parallelo: abbassalo se la RAM non basta

CHECKPOINTS <- c(
  "/home/ammendola/laboratorio/spatial/topo/code/benchmark_v7/heme_1000/checkpoint_pre_grafici.RData",
  "/home/ammendola/laboratorio/spatial/brain_injury/code/stab_inj1/benchmark_v7/stab_inj1/checkpoint_pre_grafici.RData",
  "/home/ammendola/laboratorio/spatial/brain_injury/code/stab_inj2/benchmark_v7/stab_inj2/checkpoint_pre_grafici.RData",
  "/home/ammendola/laboratorio/spatial/glioblastoma_mice/code/benchmark_v7/gbm/checkpoint_pre_grafici.RData",
  "/home/ammendola/laboratorio/spatial/periodontitis/code/benchmark_v7/perio/checkpoint_pre_grafici.RData",
  "/home/ammendola/laboratorio/spatial/tumor/full_tissue/code/benchmark_v7/head_and_neck/checkpoint_pre_grafici.RData"
)

core <- file.path(SCRIPT_DIR, "prova_nulli.R")
stopifnot(file.exists(core))
dir.create(CONFRONTO_DIR, recursive = TRUE, showWarnings = FALSE)
rscript <- file.path(R.home("bin"), "Rscript")

# --- 1. Un processo per dataset ----------------------------------------------------------
presenti <- CHECKPOINTS[file.exists(CHECKPOINTS)]
for (ck in setdiff(CHECKPOINTS, presenti)) cat("MANCA:", ck, "\n")
n_jobs  <- max(1L, min(MAX_JOBS, length(presenti)))
threads <- max(1L, floor(N_CPU / n_jobs))
cat("Dataset:", length(presenti), "| in parallelo:", n_jobs, "| processori per dataset:", threads, "\n")
lancia <- function(ck) {
  t0  <- Sys.time()
  env <- c(paste0("PROVA_THREADS=", threads), paste0("OMP_NUM_THREADS=", threads),
           paste0("OPENBLAS_NUM_THREADS=", threads), paste0("MKL_NUM_THREADS=", threads))
  out <- suppressWarnings(system2(rscript, shQuote(c(core, ck)), stdout = TRUE, stderr = TRUE, env = env))
  writeLines(out, file.path(dirname(ck), "prova_nulli_log.txt"))
  st <- attr(out, "status")
  list(ck = ck, ok = is.null(st) || st == 0, out = out,
       min = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1))
}
if (RIFAI && length(presenti) > 0) {
  t_tot <- Sys.time()
  esiti <- parallel::mclapply(presenti, lancia, mc.cores = n_jobs, mc.preschedule = FALSE)
  for (e in esiti) {
    if (inherits(e, "try-error")) { cat("\nERRORE nel processo:", as.character(e), "\n"); next }
    cat("\n#################", basename(dirname(e$ck)), "#################\n")
    cat(e$out, sep = "\n")
    cat("--> ", basename(dirname(e$ck)), ": ", if (e$ok) "OK" else "ERRORE (vedi sopra)", " in ", e$min, " min\n", sep = "")
  }
  cat("\nTempo totale:", round(as.numeric(difftime(Sys.time(), t_tot, units = "mins")), 1), "min\n")
}

# --- 2. Tabelle -------------------------------------------------------------------------------
leggi_t <- function(f) bind_rows(lapply(CHECKPOINTS, function(ck) {
  p <- file.path(dirname(ck), f); if (file.exists(p)) suppressMessages(read_csv(p, show_col_types = FALSE)) else NULL
}))
riass <- leggi_t("prova_nulli_riassunto.csv")
geni  <- leggi_t("prova_nulli_tab_geni.csv")
gsea  <- leggi_t("prova_nulli_tab_gsea.csv")
falsi <- leggi_t("prova_nulli_falsi_centri.csv")
pv    <- leggi_t("prova_nulli_pvalori_long.csv")
if (nrow(geni) == 0) stop("Nessuna tabella trovata: controlla i log prova_nulli_log.txt")
larga  <- function(d, col) as.data.frame(d %>% select(Method, dataset, all_of(col)) %>%
                                           pivot_wider(names_from = dataset, values_from = all_of(col)))
mostra <- function(titolo, d) { cat("\n>>>", titolo, "<<<\n"); print(d, row.names = FALSE) }

mostra("Riassunto per dataset", as.data.frame(riass))
fs <- falsi %>% group_by(dataset, Method) %>%
  summarise(frac_p05 = round(mean(frac_p05), 3), n_bh = round(mean(n_bh), 1), .groups = "drop")
mostra("Falsi centri: frazione di geni con p < 0.05 (calibrato = 0.05)", larga(fs, "frac_p05"))
mostra("Falsi centri: geni BH-significativi, media (calibrato = circa 0)", larga(fs, "n_bh"))
mostra("% geni significativi", larga(geni, "pct_sig"))
mostra("Golden Standard recuperato (%)", larga(geni, "GS_recall"))
mostra("Golden Standard: fold sul caso", larga(geni, "GS_fold"))
mostra("AUC sul Golden Standard", larga(geni, "AUC_GS"))
mostra("GSEA: AUC contesto vs altri sul NES", larga(gsea, "auc_ctx_nes"))
mostra("GSEA: correlazione del NES con il controllo sull'espressione", larga(gsea, "cor_con_espressione"))
mostra("GSEA: pathway con FDR < 0.05 (solo diagnostica)", larga(gsea, "n_p05"))
mostra("GSEA: pathway di contesto con FDR < 0.05 (solo diagnostica)", larga(gsea, "n_p05_ctx"))
for (x in list(list(riass, "riassunto"), list(geni, "tab_geni"), list(gsea, "tab_gsea"), list(falsi, "falsi_centri")))
  write_csv(x[[1]], file.path(CONFRONTO_DIR, paste0("nulli_", x[[2]], ".csv")))

# --- 3. Grafici della calibrazione ------------------------------------------------------------
COL_T <- c(GLM_linear = "#7570b3", GLM_linear_corr = "#1b9e77", GLM_linear_spettrale = "#e6ab02")
tema0 <- theme_minimal(base_size = 9) + theme(strip.text = element_text(face = "bold"))
if (nrow(falsi) > 0) {
  pf <- falsi %>% mutate(Method = factor(Method, levels = rev(names(COL_T)))) %>%
    ggplot(aes(x = frac_p05, y = Method, colour = Method)) +
    geom_vline(xintercept = 0.05, linetype = "dashed", colour = "grey40") +
    geom_point(size = 1.8, alpha = 0.8, position = position_jitter(height = 0.15, width = 0, seed = 1)) +
    facet_wrap(~ dataset, ncol = 3, scales = "free_x") +
    scale_colour_manual(values = COL_T, guide = "none") + tema0 +
    labs(title = "Falsi centri: frazione di geni con p < 0.05", subtitle = "un punto per falso centro | linea = 0.05 (test calibrato)",
         x = "Frazione di geni con p < 0.05", y = NULL)
  print(pf); ggsave(file.path(CONFRONTO_DIR, "nulli_00a_falsi_centri_frac.png"), pf, width = 11, height = 5.5, dpi = 150)
  pb <- falsi %>% mutate(Method = factor(Method, levels = rev(names(COL_T)))) %>%
    ggplot(aes(x = n_bh + 1, y = Method, colour = Method)) +
    geom_point(size = 1.8, alpha = 0.8, position = position_jitter(height = 0.15, width = 0, seed = 1)) + scale_x_log10() +
    facet_wrap(~ dataset, ncol = 3) + scale_colour_manual(values = COL_T, guide = "none") + tema0 +
    labs(title = "Falsi centri: geni BH-significativi", subtitle = "scala log (n + 1) | 1 = nessun gene, come atteso con un test calibrato",
         x = "Geni con FDR < 0.05 (+1)", y = NULL)
  print(pb); ggsave(file.path(CONFRONTO_DIR, "nulli_00b_falsi_centri_bh.png"), pb, width = 11, height = 5.5, dpi = 150)
}
if (nrow(pv) > 0) {
  ph <- pv %>% pivot_longer(c(p_QL, p_corr, p_spet), names_to = "Test", values_to = "p") %>% filter(!is.na(p)) %>%
    mutate(Test = recode(Test, p_QL = "GLM parametrico (QL)", p_corr = "nullo corr", p_spet = "nullo spettrale")) %>%
    ggplot(aes(x = p)) + geom_histogram(breaks = seq(0, 1, 0.05), fill = "grey55", colour = "white") +
    facet_grid(dataset ~ Test, scales = "free_y") + tema0 + theme(strip.text.y = element_text(angle = 0)) +
    labs(title = "P-value al centro vero", subtitle = "barra a sinistra = geni con segnale | piatto altrove = test calibrato sui nulli",
         x = "p", y = "geni")
  print(ph); ggsave(file.path(CONFRONTO_DIR, "nulli_00c_istogrammi_p.png"), ph, width = 10, height = 1.6 * n_distinct(pv$dataset) + 1, dpi = 150)
}

# --- 4. Golden Standard e GSEA (come la pipeline) -----------------------------------------------
ALPHA_G <- 0.05; B_GS <- 1000; B_CTX <- 1000; SEED_G <- 20261029
LIV <- c("GLM_linear", "GLM_linear_corr", "GLM_linear_spettrale", "SPATA2", "SpatialDE_orig", "Expression_control")
COL <- c(COL_T, SPATA2 = "#d95f02", SpatialDE_orig = "#e7298a", Expression_control = "grey50")
leggi_g <- function(f) bind_rows(lapply(CHECKPOINTS, function(ck) {
  p <- file.path(dirname(ck), f)
  if (!file.exists(p)) return(NULL)
  d <- suppressMessages(read_csv(p, show_col_types = FALSE))
  if (!"dataset" %in% names(d)) d$dataset <- basename(dirname(ck))
  d
}))
rk_all <- leggi_g("prova_nulli_ranks.csv")
gs_all <- leggi_g("prova_nulli_gsea.csv")
if (nrow(rk_all) == 0) stop("Manca prova_nulli_ranks.csv: rilancia prova_nulli_tutti.R con RIFAI <- TRUE")
DS  <- unique(rk_all$dataset)
LIV <- intersect(LIV, unique(c(rk_all$Method, gs_all$Method)))
rk_all <- rk_all %>% filter(tested, !is.na(rank)) %>%
  mutate(Method = factor(Method, levels = LIV), dataset = factor(dataset, levels = DS))
gs_all <- gs_all %>% mutate(Method = factor(Method, levels = LIV), dataset = factor(dataset, levels = DS))

salva_g <- function(p, nome, w, h) {
  print(p); ggsave(file.path(CONFRONTO_DIR, paste0(nome, ".png")), p, width = w, height = h, dpi = 150)
}
asse_y <- scale_y_discrete(limits = rev)
tema   <- theme_minimal(base_size = 9) + theme(strip.text = element_text(face = "bold"))
largh  <- 4 * min(3, length(DS)); alt <- function(h1) h1 * ceiling(length(DS) / 3)

# AUC di Mann-Whitney: rank basso = gene in cima
auc_mw <- function(rank, gs) {
  n1 <- sum(gs); n0 <- length(gs) - n1
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(rank(-rank)[gs]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# --- 13. Boxplot dei rank del Golden Standard -----------------------------------------------
gsg <- rk_all %>% filter(is_GS)
if (nrow(gsg) > 0) {
  med <- gsg %>% group_by(dataset, Method) %>% summarise(m = median(rank), .groups = "drop")
  p13 <- ggplot(gsg, aes(x = rank, y = Method, fill = Method)) +
    geom_boxplot(outlier.size = 0.5, width = 0.6) +
    geom_text(data = med, aes(x = m, y = Method, label = round(m)), vjust = -1.3, size = 2.4, inherit.aes = FALSE) +
    facet_wrap(~ dataset, scales = "free_x", ncol = 3) +
    scale_fill_manual(values = COL, guide = "none") + asse_y + tema +
    labs(title = "Rank dei geni del Golden Standard", subtitle = "piu' a sinistra = meglio | numero = mediana",
         x = "Rank", y = NULL)
  salva_g(p13, "nulli_13_boxplot_gs", largh, alt(3.5))

  # --- 14. ROC --------------------------------------------------------------------------------
  roc <- rk_all %>% group_by(dataset, Method) %>% arrange(rank, .by_group = TRUE) %>%
    mutate(TPR = cumsum(is_GS) / sum(is_GS), FPR = cumsum(!is_GS) / sum(!is_GS)) %>%
    filter(sum(is_GS) > 0) %>% select(dataset, Method, FPR, TPR) %>% ungroup()
  roc <- bind_rows(roc, roc %>% distinct(dataset, Method) %>% mutate(FPR = 0, TPR = 0)) %>%
    arrange(dataset, Method, FPR, TPR)
  auc_roc <- roc %>% group_by(dataset, Method) %>%
    summarise(AUC = sum(diff(FPR) * (head(TPR, -1) + tail(TPR, -1)) / 2), .groups = "drop") %>%
    group_by(dataset) %>% mutate(y = 0.04 + 0.075 * (rev(seq_len(n())) - 1)) %>% ungroup()
  p14 <- ggplot(roc, aes(x = FPR, y = TPR, colour = Method)) +
    geom_abline(linetype = "dashed", colour = "grey60") + geom_line(linewidth = 0.6) +
    geom_text(data = auc_roc, aes(x = 1, y = y, label = sprintf("%.3f", AUC)), hjust = 1, size = 2.6,
              show.legend = FALSE) +
    facet_wrap(~ dataset, ncol = 3) + scale_colour_manual(values = COL) + coord_equal() + tema +
    theme(legend.position = "bottom") +
    labs(title = "ROC sul Golden Standard", subtitle = "numeri = AUC", x = "False positive rate",
         y = "True positive rate", colour = NULL)
  salva_g(p14, "nulli_14_roc", largh, alt(4) + 0.8)

  # --- 14b. AUC con IC bootstrap (geni ricampionati, come la pipeline) ------------------------
  combos <- rk_all %>% distinct(dataset, Method)
  boot1 <- function(i) {
    d <- rk_all %>% filter(dataset == combos$dataset[i], Method == combos$Method[i])
    if (sum(d$is_GS) == 0) return(NULL)
    set.seed(SEED_G)
    bv <- replicate(B_GS, { idx <- sample.int(nrow(d), nrow(d), replace = TRUE); auc_mw(d$rank[idx], d$is_GS[idx]) })
    data.frame(dataset = combos$dataset[i], Method = combos$Method[i], AUC = auc_mw(d$rank, d$is_GS),
               lo = quantile(bv, 0.025, na.rm = TRUE), hi = quantile(bv, 0.975, na.rm = TRUE))
  }
  cat("Bootstrap AUC Golden Standard:", nrow(combos), "combinazioni x", B_GS, "repliche ...\n")
  auc_ci <- bind_rows(parallel::mclapply(seq_len(nrow(combos)), boot1, mc.cores = min(N_CPU, nrow(combos))))
  write_csv(auc_ci, file.path(CONFRONTO_DIR, "nulli_auc_gs_bootstrap.csv"))
  p14b <- ggplot(auc_ci, aes(x = AUC, y = Method, xmin = lo, xmax = hi, colour = Method)) +
    geom_pointrange(size = 0.3) + facet_wrap(~ dataset, ncol = 3) +
    scale_colour_manual(values = COL, guide = "none") + asse_y + tema +
    labs(title = "AUC sul Golden Standard con intervallo bootstrap (95%)",
         subtitle = paste0(B_GS, " repliche | intervalli marginali: la non sovrapposizione non e' un test"),
         x = "AUC", y = NULL)
  salva_g(p14b, "nulli_14b_auc_gs_bootstrap", largh, alt(3))
} else cat("Nessun gene del Golden Standard: 13, 14 e 14b saltati\n")

# --- 16. Arricchimento --------------------------------------------------------------------------
if (nrow(gs_all) > 0) {
  auc_r <- function(r, i, n1, n0) (sum(r[i]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  ctx <- bind_rows(lapply(split(gs_all, list(gs_all$dataset, gs_all$Method), drop = TRUE), function(d) {
    d <- d[!is.na(d$NES), ]; y <- d$is_context; n1 <- sum(y); n0 <- sum(!y)
    if (n1 < 5 || n0 < 5) return(NULL)
    set.seed(SEED_G)
    bt <- replicate(B_CTX, auc_r(rank(d$NES[c(sample(which(y), replace = TRUE), sample(which(!y), replace = TRUE))]),
                                 seq_len(n1), n1, n0))
    data.frame(dataset = d$dataset[1], Method = d$Method[1], auc = auc_r(rank(d$NES), which(y), n1, n0),
               lo = quantile(bt, 0.025), hi = quantile(bt, 0.975), n_ctx = n1)
  }))
  write_csv(ctx, file.path(CONFRONTO_DIR, "nulli_auc_contesto_bootstrap.csv"))
  p04 <- ggplot(ctx, aes(x = auc, y = Method, colour = Method)) +
    geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.3, orientation = "y") + geom_point(size = 2) +
    facet_wrap(~ dataset, ncol = 3) + scale_colour_manual(values = COL, guide = "none") + asse_y + tema +
    labs(title = "GSEA: pathway di contesto contro gli altri, sul NES",
         subtitle = "AUC (0.5 = nessuna preferenza) | IC 95% bootstrap sui pathway", x = "AUC", y = NULL)
  salva_g(p04, "nulli_04_auc_contesto", largh, alt(3))

  p04c <- gs_all %>% mutate(Pathway = ifelse(is_context, "Contesto", "Altri")) %>%
    ggplot(aes(x = NES, y = Method, fill = Pathway)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_boxplot(outlier.size = 0.3, linewidth = 0.25) +
    facet_wrap(~ dataset, ncol = 3, scales = "free_x") +
    scale_fill_manual(values = c(Contesto = "#d95f02", Altri = "grey75")) + asse_y + tema +
    theme(legend.position = "bottom") +
    labs(title = "GSEA: NES dei pathway di contesto e degli altri", subtitle = "nessun filtro di significativita'",
         x = "NES", y = NULL, fill = NULL)
  salva_g(p04c, "nulli_04c_nes_boxplot", largh, alt(3.2) + 0.5)

  p06 <- gs_all %>% group_by(dataset, Method) %>%
    summarise(Tutti = sum(p_adjust < ALPHA_G, na.rm = TRUE),
              Contesto = sum(p_adjust < ALPHA_G & is_context, na.rm = TRUE), .groups = "drop") %>%
    pivot_longer(c(Tutti, Contesto), names_to = "Pathway", values_to = "n") %>%
    mutate(Pathway = factor(Pathway, levels = c("Tutti", "Contesto"))) %>%
    ggplot(aes(x = n, y = Method, fill = Method)) + geom_col() +
    geom_text(aes(label = n), hjust = -0.15, size = 2.4) +
    facet_grid(dataset ~ Pathway, scales = "free_x") +
    scale_fill_manual(values = COL, guide = "none") + asse_y +
    scale_x_continuous(expand = expansion(mult = c(0, 0.2))) + tema +
    theme(strip.text.y = element_text(angle = 0)) +
    labs(title = paste0("GSEA: pathway con FDR < ", ALPHA_G, " (solo diagnostica)"), x = "Numero di pathway", y = NULL)
  salva_g(p06, "nulli_06_pathway_fdr", 9, 1.6 * length(DS) + 1)
}
cat("Grafici di tutti i dataset in", CONFRONTO_DIR, "\n")

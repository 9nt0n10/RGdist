# ##############################################################################
# ##   GRAFICI DELLA PROVA SUL P EMPIRICO, TUTTI I DATASET INSIEME            ##
# ##                                                                          ##
# ##   Gli stessi grafici della pipeline_v7 (13 boxplot, 14 ROC, 14b AUC con  ##
# ##   IC bootstrap, 16 arricchimento), un pannello per dataset.             ##
# ##   Legge i CSV scritti da prova_p_empirico.R nella cartella di ogni       ##
# ##   checkpoint: prova_p_empirico_ranks.csv e prova_p_empirico_gsea.csv.    ##
# ##   Lo lancia prova_p_empirico_tutti.R alla fine; da solo serve che        ##
# ##   CHECKPOINTS e CONFRONTO_DIR siano definiti (come in quello script).    ##
# ##############################################################################

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })
stopifnot(exists("CHECKPOINTS"), exists("CONFRONTO_DIR"))
if (!exists("N_CPU")) N_CPU <- 4

ALPHA_G <- 0.05
B_GS    <- 1000        # repliche bootstrap per l'AUC sul Golden Standard (come la 14b)
B_CTX   <- 1000        # repliche bootstrap per l'AUC di contesto (come la 16)
SEED_G  <- 20261028
LIV <- c("GLM_linear", "GLM_linear_emp", "GLM_linear_emp_abs", "SPATA2", "SpatialDE_orig", "Expression_control")
COL <- c(GLM_linear = "#7570b3", GLM_linear_emp = "#1b9e77", GLM_linear_emp_abs = "#66a61e",
         SPATA2 = "#d95f02", SpatialDE_orig = "#e7298a", Expression_control = "grey50")

leggi_g <- function(f) bind_rows(lapply(CHECKPOINTS, function(ck) {
  p <- file.path(dirname(ck), f)
  if (!file.exists(p)) return(NULL)
  d <- suppressMessages(read_csv(p, show_col_types = FALSE))
  if (!"dataset" %in% names(d)) d$dataset <- basename(dirname(ck))
  d
}))
rk_all <- leggi_g("prova_p_empirico_ranks.csv")
gs_all <- leggi_g("prova_p_empirico_gsea.csv")
if (nrow(rk_all) == 0) stop("Manca prova_p_empirico_ranks.csv: rilancia prova_p_empirico_tutti.R con RIFAI <- TRUE")
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
  salva_g(p13, "13_boxplot_gs_tutti", largh, alt(3.5))

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
  salva_g(p14, "14_roc_tutti", largh, alt(4) + 0.8)

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
  write_csv(auc_ci, file.path(CONFRONTO_DIR, "auc_gs_bootstrap_tutti.csv"))
  p14b <- ggplot(auc_ci, aes(x = AUC, y = Method, xmin = lo, xmax = hi, colour = Method)) +
    geom_pointrange(size = 0.3) + facet_wrap(~ dataset, ncol = 3) +
    scale_colour_manual(values = COL, guide = "none") + asse_y + tema +
    labs(title = "AUC sul Golden Standard con intervallo bootstrap (95%)",
         subtitle = paste0(B_GS, " repliche | intervalli marginali: la non sovrapposizione non e' un test"),
         x = "AUC", y = NULL)
  salva_g(p14b, "14b_auc_gs_bootstrap_tutti", largh, alt(3))
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
  write_csv(ctx, file.path(CONFRONTO_DIR, "auc_contesto_bootstrap_tutti.csv"))
  p04 <- ggplot(ctx, aes(x = auc, y = Method, colour = Method)) +
    geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.3, orientation = "y") + geom_point(size = 2) +
    facet_wrap(~ dataset, ncol = 3) + scale_colour_manual(values = COL, guide = "none") + asse_y + tema +
    labs(title = "GSEA: pathway di contesto contro gli altri, sul NES",
         subtitle = "AUC (0.5 = nessuna preferenza) | IC 95% bootstrap sui pathway", x = "AUC", y = NULL)
  salva_g(p04, "04_auc_contesto_tutti", largh, alt(3))

  p04c <- gs_all %>% mutate(Pathway = ifelse(is_context, "Contesto", "Altri")) %>%
    ggplot(aes(x = NES, y = Method, fill = Pathway)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_boxplot(outlier.size = 0.3, linewidth = 0.25) +
    facet_wrap(~ dataset, ncol = 3, scales = "free_x") +
    scale_fill_manual(values = c(Contesto = "#d95f02", Altri = "grey75")) + asse_y + tema +
    theme(legend.position = "bottom") +
    labs(title = "GSEA: NES dei pathway di contesto e degli altri", subtitle = "nessun filtro di significativita'",
         x = "NES", y = NULL, fill = NULL)
  salva_g(p04c, "04c_nes_boxplot_tutti", largh, alt(3.2) + 0.5)

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
  salva_g(p06, "06_pathway_fdr_tutti", 9, 1.6 * length(DS) + 1)
}
cat("Grafici di tutti i dataset in", CONFRONTO_DIR, "\n")

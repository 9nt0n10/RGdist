# ##############################################################################
# ##   CONFRONTO FRA DUE RUN DELLA PIPELINE: rm_zero_infl TRUE contro FALSE   ##
# ##                                                                          ##
# ##   Script autonomo: legge i CSV che pipeline_v7 scrive in OUTPUT_DIR      ##
# ##     risultati_tutti_metodi_long.csv  (rank, fdr, Golden Standard)        ##
# ##     enrichment_gsea_pathways.csv     (NES e FDR di ogni pathway)         ##
# ##   e mette a confronto i due run, metodo per metodo:                      ##
# ##     - geni nella lista, geni significativi                               ##
# ##     - Golden Standard: geni nella lista, recuperati (FDR < ALPHA), AUC   ##
# ##       sulla propria lista e con i GS non testati in fondo, rank relativo ##
# ##     - GSEA: pathway con FDR < ALPHA, di contesto, frazione di contesto,  ##
# ##       AUC contesto/altri sul NES (IC bootstrap), NES mediani,            ##
# ##       correlazione con il controllo sull'espressione                     ##
# ##     - SPATA2 nel dettaglio: geni del GS persi, rank prima e dopo,        ##
# ##       sovrapposizione dei significativi, NES dei pathway nei due run     ##
# ##   Fra i due run cambia solo SPATA2: per gli altri metodi le differenze   ##
# ##   devono essere zero (controllo di coerenza).                            ##
# ##   Tutto in OUT_DIR: tabelle CSV e grafici PNG.                           ##
# ##############################################################################

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })

# --- da impostare --------------------------------------------------------------------------
# Per ogni dataset le due cartelle di output (OUTPUT_DIR della pipeline) dei due run.
RUNS <- list(
  perio = c(zero_infl_TRUE  = "/home/ammendola/laboratorio/spatial/periodontitis/code/benchmark_v7/perio_zi_true",
            zero_infl_FALSE = "/home/ammendola/laboratorio/spatial/periodontitis/code/benchmark_v7/perio")
)
OUT_DIR <- "/home/ammendola/laboratorio/spatial/test_arricchimenti/confronto_zero_infl"
ALPHA   <- 0.05
B_BOOT  <- 1000
SEED    <- 20261030
# ---------------------------------------------------------------------------------------------

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
RUN_LIV <- c("zero_infl_FALSE", "zero_infl_TRUE")
COL_RUN <- c(zero_infl_FALSE = "#1b9e77", zero_infl_TRUE = "#d95f02")

leggi <- function(dir, f) {
  p <- file.path(dir, f)
  if (!file.exists(p)) { cat("MANCA:", p, "\n"); return(NULL) }
  suppressMessages(read_csv(p, show_col_types = FALSE))
}
long <- bind_rows(lapply(names(RUNS), function(ds) bind_rows(lapply(names(RUNS[[ds]]), function(r) {
  d <- leggi(RUNS[[ds]][[r]], "risultati_tutti_metodi_long.csv")
  if (is.null(d)) return(NULL)
  d %>% mutate(dataset = ds, run = r, Method = as.character(Method),
               is_GS = is_GS %in% TRUE, sig = !is.na(fdr) & fdr < ALPHA)
}))))
gsea <- bind_rows(lapply(names(RUNS), function(ds) bind_rows(lapply(names(RUNS[[ds]]), function(r) {
  d <- leggi(RUNS[[ds]][[r]], "enrichment_gsea_pathways.csv")
  if (is.null(d)) return(NULL)
  d %>% mutate(dataset = ds, run = r, Method = as.character(Method), is_context = is_context %in% TRUE)
}))))
stopifnot(nrow(long) > 0)
long <- long %>% filter(!is.na(rank)) %>% mutate(run = factor(run, levels = RUN_LIV))
if (nrow(gsea) > 0) gsea <- gsea %>% mutate(run = factor(run, levels = RUN_LIV))

auc_mw <- function(rank, gs) {
  n1 <- sum(gs); n0 <- sum(!gs)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(rank(-rank)[gs]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
mostra <- function(titolo, d) { cat("\n>>>", titolo, "<<<\n"); print(as.data.frame(d), row.names = FALSE, digits = 3) }
salva  <- function(p, nome, w, h) { print(p); ggsave(file.path(OUT_DIR, paste0(nome, ".png")), p, width = w, height = h, dpi = 150) }
tema   <- theme_minimal(base_size = 9) + theme(strip.text = element_text(face = "bold"), legend.position = "bottom")
# posizione verticale di un metodo, spostata in alto o in basso secondo il run
y_run <- function(Method, run, lv) as.numeric(factor(Method, levels = lv)) + c(-0.18, 0.18)[as.integer(run)]
asse_metodi <- function(lv) scale_y_continuous(breaks = seq_along(lv), labels = lv, expand = expansion(add = 0.6))
ordine_metodi <- function(d, col) d %>% group_by(Method) %>% summarise(m = mean(.data[[col]], na.rm = TRUE)) %>%
  arrange(m) %>% pull(Method)

# Golden Standard di riferimento per ogni dataset e run: i geni GS presenti in almeno una lista
gs_rif <- long %>% filter(is_GS) %>% distinct(dataset, run, gene_id, gene)

# ===============================================================================================
# 1. Geni, Golden Standard, AUC
# ===============================================================================================
tab <- long %>% group_by(dataset, run, Method) %>% group_modify(function(d, k) {
  gs_tot <- gs_rif$gene_id[gs_rif$dataset == k$dataset & gs_rif$run == k$run]
  mn     <- setdiff(gs_tot, d$gene_id)
  rk_pen <- c(d$rank, max(d$rank) + seq_along(mn)); gs_pen <- c(d$is_GS, rep(TRUE, length(mn)))
  rr_gs  <- c(d$rank[d$is_GS] / nrow(d), rep(1, length(mn)))
  tibble(geni_lista = nrow(d), n_sig = sum(d$sig), pct_sig = round(100 * mean(d$sig), 1),
         GS_totali = length(gs_tot), GS_in_lista = sum(d$is_GS), GS_mancanti = length(mn),
         GS_recuperati = sum(d$sig & d$is_GS), GS_recall = round(100 * sum(d$sig & d$is_GS) / max(length(gs_tot), 1), 1),
         GS_fold = round((sum(d$sig & d$is_GS) / max(sum(d$sig), 1)) / (sum(d$is_GS) / nrow(d)), 2),
         AUC_lista = round(auc_mw(d$rank, d$is_GS), 3),
         AUC_GS_mancanti_in_fondo = round(auc_mw(rk_pen, gs_pen), 3),
         rank_rel_GS_mediana = round(median(rr_gs), 3))
}) %>% ungroup()
write_csv(tab, file.path(OUT_DIR, "tab_geni_gs_auc.csv"))

larga_run <- function(col) tab %>% select(dataset, Method, run, all_of(col)) %>%
  pivot_wider(names_from = run, values_from = all_of(col)) %>%
  mutate(differenza = .data[["zero_infl_TRUE"]] - .data[["zero_infl_FALSE"]]) %>% arrange(dataset, desc(abs(differenza)))
for (col in c("geni_lista", "n_sig", "GS_in_lista", "GS_recuperati", "GS_recall", "AUC_lista",
              "AUC_GS_mancanti_in_fondo", "rank_rel_GS_mediana"))
  mostra(paste0(col, " (TRUE - FALSE)"), larga_run(col))

ctrl <- larga_run("AUC_lista") %>% filter(Method != "SPATA2", abs(differenza) > 1e-6)
cat("\nControllo di coerenza: metodi diversi da SPATA2 che cambiano fra i due run:",
    if (nrow(ctrl) == 0) "nessuno (ok)" else paste(unique(ctrl$Method), collapse = ", "), "\n")

# --- grafici --------------------------------------------------------------------------------------
lv_auc <- ordine_metodi(tab, "AUC_GS_mancanti_in_fondo")
p1 <- tab %>% select(dataset, run, Method, `sulla propria lista` = AUC_lista,
                     `GS non testati in fondo` = AUC_GS_mancanti_in_fondo) %>%
  pivot_longer(-c(dataset, run, Method), names_to = "AUC", values_to = "v") %>%
  mutate(yy = y_run(Method, run, lv_auc)) %>%
  ggplot(aes(x = v, y = yy, colour = run, shape = AUC)) +
  geom_point(size = 2) + asse_metodi(lv_auc) +
  facet_wrap(~ dataset) + scale_colour_manual(values = COL_RUN) + tema +
  labs(title = "AUC sul Golden Standard", subtitle = "con i geni del GS non testati in fondo = confronto equo fra liste diverse",
       x = "AUC", y = NULL, colour = NULL, shape = NULL)
salva(p1, "01_auc_gs", 9, 0.3 * length(lv_auc) + 2.5)

rr <- long %>% group_by(dataset, run, Method) %>% group_modify(function(d, k) {
  gs_tot <- gs_rif$gene_id[gs_rif$dataset == k$dataset & gs_rif$run == k$run]
  mn <- setdiff(gs_tot, d$gene_id)
  bind_rows(tibble(rank_rel = d$rank[d$is_GS] / nrow(d), stato = "testato"),
            tibble(rank_rel = rep(1, length(mn)), stato = rep("non testato (in fondo)", length(mn))))
}) %>% ungroup() %>% mutate(Method = factor(Method, levels = rev(ordine_metodi(tab, "rank_rel_GS_mediana"))))
p2 <- ggplot(rr, aes(x = rank_rel, y = Method, fill = run)) +
  geom_boxplot(outlier.size = 0.4, linewidth = 0.3, position = position_dodge(width = 0.75)) +
  facet_wrap(~ dataset) + scale_fill_manual(values = COL_RUN) + scale_x_continuous(labels = scales::percent) + tema +
  labs(title = "Rank relativo dei geni del Golden Standard", subtitle = "posizione / lunghezza della lista | 0% = in cima | GS non testati al 100%",
       x = "Rank relativo", y = NULL, fill = NULL)
salva(p2, "02_boxplot_rank_relativo_gs", 9, 0.3 * n_distinct(rr$Method) + 2.5)

lv_gs <- ordine_metodi(tab, "GS_recuperati")
p3 <- tab %>% select(dataset, run, Method, `GS nella lista` = GS_in_lista, `GS recuperati (FDR < 0.05)` = GS_recuperati) %>%
  pivot_longer(-c(dataset, run, Method), names_to = "Misura", values_to = "n") %>%
  mutate(yy = y_run(Method, run, lv_gs)) %>%
  ggplot(aes(x = n, y = yy, fill = run)) + geom_col(orientation = "y", width = 0.34) +
  geom_text(aes(label = n), hjust = -0.15, size = 2.3) + asse_metodi(lv_gs) +
  facet_grid(dataset ~ Misura, scales = "free_x") + scale_fill_manual(values = COL_RUN) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) + tema +
  labs(title = "Geni del Golden Standard", subtitle = "totale del GS = geni presenti in almeno una lista", x = "Geni", y = NULL, fill = NULL)
salva(p3, "03_gs_recuperati", 10, 0.3 * length(lv_auc) + 2.5)

lv_sig <- ordine_metodi(tab, "n_sig")
p4 <- tab %>% select(dataset, run, Method, `Geni nella lista` = geni_lista, `Geni significativi` = n_sig) %>%
  pivot_longer(-c(dataset, run, Method), names_to = "Misura", values_to = "n") %>%
  mutate(yy = y_run(Method, run, lv_sig)) %>%
  ggplot(aes(x = n, y = yy, fill = run)) + geom_col(orientation = "y", width = 0.34) + asse_metodi(lv_sig) +
  facet_grid(dataset ~ Misura, scales = "free_x") + scale_fill_manual(values = COL_RUN) + tema +
  labs(title = "Geni nella lista e geni significativi", x = "Geni", y = NULL, fill = NULL)
salva(p4, "04_geni_lista_significativi", 10, 0.3 * length(lv_auc) + 2.5)

# ===============================================================================================
# 2. Arricchimento
# ===============================================================================================
if (nrow(gsea) > 0) {
  nes_expr <- gsea %>% filter(Method == "Expression_control") %>% select(dataset, run, pathway, NES_expr = NES)
  auc_r <- function(nes, ctx) {
    ok <- !is.na(nes); nes <- nes[ok]; ctx <- ctx[ok]; n1 <- sum(ctx); n0 <- sum(!ctx)
    if (n1 < 5 || n0 < 5) return(NA_real_)
    (sum(rank(nes)[ctx]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  }
  tab_e <- gsea %>% left_join(nes_expr, by = c("dataset", "run", "pathway")) %>%
    group_by(dataset, run, Method) %>% group_modify(function(d, k) {
      y <- d$is_context & !is.na(d$NES); x <- !d$is_context & !is.na(d$NES)
      set.seed(SEED)
      bt <- replicate(B_BOOT, auc_r(c(sample(d$NES[y], replace = TRUE), sample(d$NES[x], replace = TRUE)),
                                    c(rep(TRUE, sum(y)), rep(FALSE, sum(x)))))
      tibble(pathway_testati = sum(!is.na(d$NES)), pathway_contesto = sum(y),
             p05_tutti = sum(d$p_adjust < ALPHA, na.rm = TRUE),
             p05_contesto = sum(d$p_adjust < ALPHA & d$is_context, na.rm = TRUE),
             frazione_contesto = round(sum(d$p_adjust < ALPHA & d$is_context, na.rm = TRUE) /
                                         max(sum(d$p_adjust < ALPHA, na.rm = TRUE), 1), 3),
             frazione_contesto_attesa = round(mean(d$is_context[!is.na(d$NES)]), 3),
             AUC_contesto = round(auc_r(d$NES, d$is_context), 3),
             AUC_contesto_lo = round(quantile(bt, 0.025, na.rm = TRUE), 3),
             AUC_contesto_hi = round(quantile(bt, 0.975, na.rm = TRUE), 3),
             NES_contesto_mediana = round(median(d$NES[y]), 2), NES_altri_mediana = round(median(d$NES[x]), 2),
             cor_con_espressione = round(suppressWarnings(cor(d$NES, d$NES_expr, method = "spearman",
                                                              use = "pairwise.complete.obs")), 2))
    }) %>% ungroup()
  write_csv(tab_e, file.path(OUT_DIR, "tab_arricchimento.csv"))
  larga_e <- function(col) tab_e %>% select(dataset, Method, run, all_of(col)) %>%
    pivot_wider(names_from = run, values_from = all_of(col)) %>%
    mutate(differenza = .data[["zero_infl_TRUE"]] - .data[["zero_infl_FALSE"]]) %>% arrange(dataset, desc(abs(differenza)))
  for (col in c("p05_tutti", "p05_contesto", "frazione_contesto", "AUC_contesto", "NES_contesto_mediana", "cor_con_espressione"))
    mostra(paste0("GSEA: ", col, " (TRUE - FALSE)"), larga_e(col))
  cat("\nfrazione_contesto_attesa = quota di pathway di contesto fra tutti quelli testati (quella che si avrebbe a caso)\n")

  lv_e <- ordine_metodi(tab_e, "AUC_contesto")
  p5 <- tab_e %>% mutate(yy = y_run(Method, run, lv_e)) %>%
    ggplot(aes(x = AUC_contesto, y = yy, colour = run)) +
    geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(xmin = AUC_contesto_lo, xmax = AUC_contesto_hi), width = 0.15, orientation = "y") + geom_point(size = 1.8) +
    asse_metodi(lv_e) +
    facet_wrap(~ dataset) + scale_colour_manual(values = COL_RUN) + tema +
    labs(title = "GSEA: pathway di contesto contro gli altri, sul NES", subtitle = "0.5 = nessuna preferenza | IC 95% bootstrap sui pathway",
         x = "AUC", y = NULL, colour = NULL)
  salva(p5, "05_auc_contesto", 9, 0.3 * length(lv_e) + 2.5)

  p6 <- tab_e %>% select(dataset, run, Method, `Pathway con FDR < 0.05` = p05_tutti, `Di contesto` = p05_contesto,
                         `Frazione di contesto` = frazione_contesto) %>%
    pivot_longer(-c(dataset, run, Method), names_to = "Misura", values_to = "v") %>%
    mutate(yy = y_run(Method, run, lv_e),
           Misura = factor(Misura, levels = c("Pathway con FDR < 0.05", "Di contesto", "Frazione di contesto"))) %>%
    ggplot(aes(x = v, y = yy, fill = run)) + geom_col(orientation = "y", width = 0.34) + asse_metodi(lv_e) +
    facet_grid(dataset ~ Misura, scales = "free_x") + scale_fill_manual(values = COL_RUN) + tema +
    labs(title = "GSEA: pathway con FDR < 0.05 e frazione di contesto (solo descrittivo)", x = NULL, y = NULL, fill = NULL)
  salva(p6, "06_pathway_fdr_frazione", 12, 0.3 * length(lv_e) + 2.5)

  p7 <- gsea %>% mutate(Pathway = ifelse(is_context, "Contesto", "Altri"), Method = factor(Method, levels = lv_e)) %>%
    ggplot(aes(x = NES, y = Method, fill = interaction(Pathway, run, sep = " | "))) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_boxplot(outlier.size = 0.3, linewidth = 0.25) + facet_wrap(~ dataset) +
    scale_fill_manual(values = c("Contesto | zero_infl_FALSE" = "#1b9e77", "Altri | zero_infl_FALSE" = "#b3e2cd",
                                 "Contesto | zero_infl_TRUE" = "#d95f02", "Altri | zero_infl_TRUE" = "#fdcdac")) + tema +
    labs(title = "GSEA: NES dei pathway di contesto e degli altri", x = "NES", y = NULL, fill = NULL)
  salva(p7, "07_nes_boxplot", 10, 0.35 * length(lv_e) + 2.5)
}

# ===============================================================================================
# 3. SPATA2 nel dettaglio
# ===============================================================================================
for (ds in names(RUNS)) {
  sF <- long %>% filter(dataset == ds, run == "zero_infl_FALSE", Method == "SPATA2")
  sT <- long %>% filter(dataset == ds, run == "zero_infl_TRUE", Method == "SPATA2")
  if (nrow(sF) == 0 || nrow(sT) == 0) { cat("\n", ds, ": SPATA2 manca in uno dei due run\n"); next }
  cat("\n================ SPATA2,", ds, "================\n")
  cat("Geni nella lista: FALSE", nrow(sF), "| TRUE", nrow(sT), "| tolti dal filtro", length(setdiff(sF$gene_id, sT$gene_id)), "\n")
  sigF <- sF$gene_id[sF$sig]; sigT <- sT$gene_id[sT$sig]
  cat("Significativi: FALSE", length(sigF), "| TRUE", length(sigT), "| in comune", length(intersect(sigF, sigT)),
      "| Jaccard", round(length(intersect(sigF, sigT)) / max(length(union(sigF, sigT)), 1), 3), "\n")
  com <- intersect(sF$gene_id, sT$gene_id)
  cat("Spearman dei rank sui geni comuni:", round(cor(sF$rank[match(com, sF$gene_id)], sT$rank[match(com, sT$gene_id)],
                                                      method = "spearman"), 3), "\n")

  persi <- sF %>% filter(is_GS, !gene_id %in% sT$gene_id) %>% mutate(rank_rel_FALSE = round(rank / nrow(sF), 3)) %>%
    arrange(rank) %>% select(gene, gene_id, rank, rank_rel_FALSE, fdr, sig)
  tenuti <- sF %>% filter(is_GS, gene_id %in% sT$gene_id)
  cat("Geni del GS tolti dal filtro:", nrow(persi), "| rank relativo mediano nel run FALSE:",
      round(median(persi$rank_rel_FALSE), 3), "(GS tenuti:", round(median(tenuti$rank / nrow(sF)), 3), ")\n")
  cat("Di questi, significativi nel run FALSE:", sum(persi$sig), "\n")
  write_csv(persi, file.path(OUT_DIR, paste0("spata2_gs_persi_", ds, ".csv")))
  mostra(paste("Geni del GS tolti dal filtro (primi 40 per rank nel run FALSE),", ds), head(persi, 40))

  pr <- tibble(gene_id = com, rr_F = sF$rank[match(com, sF$gene_id)] / nrow(sF),
               rr_T = sT$rank[match(com, sT$gene_id)] / nrow(sT), GS = sF$is_GS[match(com, sF$gene_id)]) %>% arrange(GS)
  p8 <- ggplot(pr, aes(x = rr_F, y = rr_T, colour = GS)) + geom_point(size = 0.4, alpha = 0.5) +
    geom_abline(linetype = "dashed", colour = "grey40") + coord_equal() +
    scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "#d95f02")) +
    scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) + tema +
    labs(title = paste0(ds, ": SPATA2, rank relativo con e senza zero inflation"),
         subtitle = "geni testati in entrambi i run | arancio = Golden Standard",
         x = "rank relativo, zero_infl FALSE", y = "rank relativo, zero_infl TRUE", colour = "GS")
  salva(p8, paste0("08_spata2_rank_", ds), 6.5, 6.5)

  p9 <- bind_rows(tibble(rr = persi$rank_rel_FALSE, gruppo = "GS tolti dal filtro"),
                  tibble(rr = tenuti$rank / nrow(sF), gruppo = "GS tenuti")) %>%
    ggplot(aes(x = rr, fill = gruppo)) + geom_histogram(breaks = seq(0, 1, 0.05), colour = "white", position = "identity", alpha = 0.6) +
    scale_fill_manual(values = c("GS tolti dal filtro" = "#d95f02", "GS tenuti" = "#1b9e77")) +
    scale_x_continuous(labels = scales::percent) + tema +
    labs(title = paste0(ds, ": dove stavano, senza filtro, i geni del GS che il filtro toglie"),
         subtitle = "rank relativo nel run zero_infl FALSE", x = "Rank relativo", y = "geni", fill = NULL)
  salva(p9, paste0("09_spata2_gs_persi_", ds), 7, 4.5)

  if (nrow(gsea) > 0) {
    g2 <- gsea %>% filter(dataset == ds, Method == "SPATA2") %>% select(run, pathway, NES, p_adjust, is_context) %>%
      pivot_wider(names_from = run, values_from = c(NES, p_adjust))
    if (all(c("NES_zero_infl_FALSE", "NES_zero_infl_TRUE") %in% names(g2))) {
      p10 <- ggplot(g2 %>% arrange(is_context), aes(x = NES_zero_infl_FALSE, y = NES_zero_infl_TRUE, colour = is_context)) +
        geom_abline(linetype = "dashed", colour = "grey40") + geom_point(size = 0.7, alpha = 0.6) + coord_equal() +
        scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "#d95f02"), labels = c("altri", "contesto")) + tema +
        labs(title = paste0(ds, ": SPATA2, NES dei pathway con e senza zero inflation"),
             x = "NES, zero_infl FALSE", y = "NES, zero_infl TRUE", colour = NULL)
      salva(p10, paste0("10_spata2_nes_", ds), 6.5, 6.5)
      cambi <- g2 %>% filter(is_context) %>%
        mutate(sig_F = p_adjust_zero_infl_FALSE < ALPHA, sig_T = p_adjust_zero_infl_TRUE < ALPHA) %>%
        filter(sig_F != sig_T) %>% arrange(desc(abs(NES_zero_infl_TRUE - NES_zero_infl_FALSE)))
      mostra(paste("SPATA2: pathway di contesto che cambiano significativita',", ds), cambi)
    }
  }
}
cat("\nTabelle e grafici in", OUT_DIR, "\n")

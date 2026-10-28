# ##############################################################################
# ##   DIAGNOSTICA NUMERICA FRA DATASET (lettura dei risultati salvati)       ##
# ##                                                                          ##
# ##   Non rifitta niente: legge i CSV scritti da ogni run in                 ##
# ##   BASE_RESULTS/<SAMPLE_ID>/ e mette a confronto i numeri.                ##
# ##     A. geni testati e significativi (fdr < ALPHA) per metodo e dataset,  ##
# ##        e quanti di essi sono del Golden Standard                         ##
# ##     B. GSEA: NES contesto vs altri, AUC sul NES, correlazione del NES    ##
# ##        con quello del controllo sull'espressione                         ##
# ##     C. nullo spaziale (10f), AUC sul Golden Standard, calibrazione (10e) ##
# ##   Gli oggetti restano in memoria (diag_*) per guardarli meglio.          ##
# ##############################################################################

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })

BASE_RESULTS <- "/home/ammendola/laboratorio/spatial/topo/code/benchmark_v7"
ALPHA        <- 0.05
DATASETS     <- NULL      # NULL = tutte le cartelle con risultati; oppure c("heme_1000", ...)

if (is.null(DATASETS)) {
  dd <- list.dirs(BASE_RESULTS, recursive = FALSE, full.names = FALSE)
  DATASETS <- dd[file.exists(file.path(BASE_RESULTS, dd, "risultati_tutti_metodi_long.csv"))]
}
cat("Dataset trovati:", paste(DATASETS, collapse = ", "), "\n")

leggi <- function(ds, f) {
  p <- file.path(BASE_RESULTS, ds, f)
  if (!file.exists(p)) return(NULL)
  d <- suppressMessages(read_csv(p, show_col_types = FALSE)); d$dataset <- ds; d
}
leggi_tutti <- function(f) bind_rows(lapply(DATASETS, leggi, f = f))

# --- A. Geni ---------------------------------------------------------------------
long <- leggi_tutti("risultati_tutti_metodi_long.csv")
if (!"is_GS" %in% names(long)) long$is_GS <- NA
if (!"tested" %in% names(long)) long$tested <- TRUE
long$tested <- ifelse(is.na(long$tested), TRUE, long$tested)

diag_gene <- long %>% filter(tested) %>%
  group_by(dataset, Method) %>%
  summarise(n_tested   = n(),
            n_sig      = sum(!is.na(fdr) & fdr < ALPHA),
            pct_sig    = round(100 * n_sig / n_tested, 1),
            n_GS       = sum(is_GS %in% TRUE),
            GS_in_sig  = sum(is_GS %in% TRUE & !is.na(fdr) & fdr < ALPHA),
            .groups = "drop") %>%
  mutate(GS_base_rate = round(100 * n_GS / n_tested, 1),
         GS_pct_in_sig = round(100 * GS_in_sig / pmax(n_sig, 1), 1),
         GS_fold = round(GS_pct_in_sig / GS_base_rate, 2))

cat("\n>>> A1. % di geni testati con FDR <", ALPHA, "(righe = metodi, colonne = dataset) <<<\n")
print(as.data.frame(diag_gene %>% select(Method, dataset, pct_sig) %>%
                      pivot_wider(names_from = dataset, values_from = pct_sig)), row.names = FALSE)
cat("\n>>> A2. Numero di geni significativi <<<\n")
print(as.data.frame(diag_gene %>% select(Method, dataset, n_sig) %>%
                      pivot_wider(names_from = dataset, values_from = n_sig)), row.names = FALSE)
cat("\n>>> A3. Geni testati <<<\n")
print(as.data.frame(diag_gene %>% select(Method, dataset, n_tested) %>%
                      pivot_wider(names_from = dataset, values_from = n_tested)), row.names = FALSE)
cat("\n>>> A4. Golden Standard fra i significativi: % e fold rispetto al caso (1 = come a caso) <<<\n")
print(as.data.frame(diag_gene %>% select(Method, dataset, GS_fold) %>%
                      pivot_wider(names_from = dataset, values_from = GS_fold)), row.names = FALSE)

p_sig <- ggplot(diag_gene, aes(x = pct_sig, y = reorder(Method, pct_sig), colour = dataset)) +
  geom_point(size = 2) + theme_minimal(base_size = 9) +
  labs(title = paste0("% geni testati con FDR < ", ALPHA), x = "% geni", y = NULL, colour = NULL)
print(p_sig)

# --- B. GSEA ----------------------------------------------------------------------
gsea <- leggi_tutti("enrichment_gsea_pathways.csv")
if (!is.null(gsea) && nrow(gsea) > 0) {
  auc_ctx <- function(nes, ctx) {
    ok <- !is.na(nes); nes <- nes[ok]; ctx <- ctx[ok]
    n1 <- sum(ctx); n0 <- sum(!ctx)
    if (n1 < 5 || n0 < 5) return(NA_real_)
    (sum(rank(nes)[ctx]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  }
  nes_expr <- gsea %>% filter(Method == "Expression_control") %>%
    select(dataset, pathway, NES_expr = NES)

  diag_gsea <- gsea %>% left_join(nes_expr, by = c("dataset", "pathway")) %>%
    group_by(dataset, Method) %>%
    summarise(n_pathway     = n(),
              n_ctx         = sum(is_context),
              nes_ctx_med   = round(median(NES[is_context], na.rm = TRUE), 2),
              nes_other_med = round(median(NES[!is_context], na.rm = TRUE), 2),
              auc_ctx_nes   = round(auc_ctx(NES, is_context), 3),
              cor_con_espressione = round(suppressWarnings(cor(NES, NES_expr, method = "spearman",
                                                              use = "pairwise.complete.obs")), 2),
              n_p05_tutti   = sum(p_adjust < ALPHA, na.rm = TRUE),
              n_p05_ctx     = sum(p_adjust < ALPHA & is_context, na.rm = TRUE),
              .groups = "drop")

  cat("\n>>> B1. AUC contesto vs altri sul NES (0.5 = nessuna preferenza) <<<\n")
  print(as.data.frame(diag_gsea %>% select(Method, dataset, auc_ctx_nes) %>%
                        pivot_wider(names_from = dataset, values_from = auc_ctx_nes)), row.names = FALSE)
  cat("\n>>> B2. Correlazione di Spearman del NES con quello del controllo sull'espressione <<<\n")
  cat("    (vicino a 1 = i pathway del metodo seguono soprattutto il livello di espressione)\n")
  print(as.data.frame(diag_gsea %>% select(Method, dataset, cor_con_espressione) %>%
                        pivot_wider(names_from = dataset, values_from = cor_con_espressione)), row.names = FALSE)
  cat("\n>>> B3. Pathway con FDR <", ALPHA, "(solo diagnostica): tutti / di contesto <<<\n")
  print(as.data.frame(diag_gsea %>% mutate(n = paste0(n_p05_tutti, " / ", n_p05_ctx)) %>%
                        select(Method, dataset, n) %>%
                        pivot_wider(names_from = dataset, values_from = n)), row.names = FALSE)
  cat("\n>>> B4. NES mediano: contesto / altri <<<\n")
  print(as.data.frame(diag_gsea %>% mutate(n = paste0(nes_ctx_med, " / ", nes_other_med)) %>%
                        select(Method, dataset, n) %>%
                        pivot_wider(names_from = dataset, values_from = n)), row.names = FALSE)
} else cat("\nNessun enrichment_gsea_pathways.csv trovato.\n")

# --- C. Nullo spaziale, AUC, calibrazione --------------------------------------------
diag_null <- leggi_tutti("nullo_spaziale_centri.csv")
if (!is.null(diag_null) && nrow(diag_null) > 0) {
  cat("\n>>> C1. Nullo spaziale (10f): geni significativi al barcode vs altri centri <<<\n")
  print(as.data.frame(diag_null %>% select(dataset, Forma, n_sig_barcode, mediana_tutti, p95_tutti,
                                           p_empirico_tutti, mediana_lontani, p_empirico_lontani)),
        row.names = FALSE, digits = 3)
}
diag_auc <- leggi_tutti("roc_auc.csv")
if (!is.null(diag_auc) && nrow(diag_auc) > 0) {
  cat("\n>>> C2. AUC sul Golden Standard <<<\n")
  print(as.data.frame(diag_auc %>% select(Method, dataset, AUC) %>%
                        pivot_wider(names_from = dataset, values_from = AUC)), row.names = FALSE)
}
diag_calib <- leggi_tutti("calibrazione_profondita_score_test.csv")
if (!is.null(diag_calib) && nrow(diag_calib) > 0) {
  cat("\n>>> C3. Calibrazione sotto ipotesi nulla (10e): frac_p05_null vicino a 0.05 = calibrato <<<\n")
  print(as.data.frame(diag_calib), row.names = FALSE, digits = 3)
}

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
  mutate(GS_base_rate = round(100 * n_GS / n_tested, 2),
         GS_pct_in_sig = round(100 * GS_in_sig / pmax(n_sig, 1), 2),
         GS_fold   = round((GS_in_sig / pmax(n_sig, 1)) / (n_GS / n_tested), 2),
         GS_recall = round(100 * GS_in_sig / pmax(n_GS, 1), 1))

cat("\n>>> A1. % di geni testati con FDR <", ALPHA, "(righe = metodi, colonne = dataset) <<<\n")
print(as.data.frame(diag_gene %>% select(Method, dataset, pct_sig) %>%
                      pivot_wider(names_from = dataset, values_from = pct_sig)), row.names = FALSE)
cat("\n>>> A2. Numero di geni significativi <<<\n")
print(as.data.frame(diag_gene %>% select(Method, dataset, n_sig) %>%
                      pivot_wider(names_from = dataset, values_from = n_sig)), row.names = FALSE)
cat("\n>>> A3. Geni testati <<<\n")
print(as.data.frame(diag_gene %>% select(Method, dataset, n_tested) %>%
                      pivot_wider(names_from = dataset, values_from = n_tested)), row.names = FALSE)
cat("\n>>> A4. Golden Standard fra i significativi: fold rispetto al caso (1 = come a caso) <<<\n")
print(as.data.frame(diag_gene %>% select(Method, dataset, GS_fold) %>%
                      pivot_wider(names_from = dataset, values_from = GS_fold)), row.names = FALSE)
cat("\n>>> A5. Golden Standard recuperato (% dei geni GS con FDR <", ALPHA, ") <<<\n")
print(as.data.frame(diag_gene %>% select(Method, dataset, GS_recall) %>%
                      pivot_wider(names_from = dataset, values_from = GS_recall)), row.names = FALSE)

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

# --- D. SPATA2 contro gli altri metodi -----------------------------------------------
# sig_in_comune: geni significativi (fdr < ALPHA) sia per SPATA2 sia per l'altro metodo.
# overlap_top_N: i primi N geni dell'altro metodo (N = n. significativi di SPATA2) che
#   sono anche significativi per SPATA2: confronto alla pari, senza soglia dell'altro.
# spearman_rank: correlazione fra i rank dei due metodi sui geni testati da entrambi.
REF_M <- "SPATA2"
diag_spata <- bind_rows(lapply(DATASETS, function(ds) {
  d   <- long %>% filter(dataset == ds, tested, !is.na(rank))
  ref <- d %>% filter(Method == REF_M)
  if (nrow(ref) == 0) return(NULL)
  sig_ref <- ref$gene_id[!is.na(ref$fdr) & ref$fdr < ALPHA]
  n_top   <- length(sig_ref)
  bind_rows(lapply(setdiff(unique(d$Method), REF_M), function(m) {
    o      <- d %>% filter(Method == m)
    comuni <- intersect(o$gene_id, ref$gene_id)
    if (length(comuni) < 100) return(NULL)
    sig_o  <- o$gene_id[!is.na(o$fdr) & o$fdr < ALPHA]
    top_o  <- o$gene_id[order(o$rank)][seq_len(min(n_top, nrow(o)))]
    inter  <- length(intersect(sig_ref, sig_o))
    data.frame(dataset = ds, Method = m, n_sig_SPATA2 = n_top, n_sig_metodo = length(sig_o),
               sig_in_comune = inter,
               pct_SPATA2_anche_nel_metodo = round(100 * inter / max(n_top, 1), 1),
               jaccard_sig = round(inter / max(length(union(sig_ref, sig_o)), 1), 3),
               overlap_top_N = round(100 * length(intersect(sig_ref, top_o)) / max(n_top, 1), 1),
               spearman_rank = round(cor(ref$rank[match(comuni, ref$gene_id)],
                                         o$rank[match(comuni, o$gene_id)], method = "spearman"), 3),
               stringsAsFactors = FALSE)
  }))
}))
cat("\n>>> D1. SPATA2 contro gli altri metodi: sovrapposizione dei significativi e dei rank <<<\n")
print(as.data.frame(diag_spata %>% arrange(dataset, desc(spearman_rank))), row.names = FALSE)

# --- D2. Risoluzione dei p-value di SPATA2 -------------------------------------------
for (ds in DATASETS) {
  f <- file.path(BASE_RESULTS, ds, "spata2", "spata2_results.csv")
  if (!file.exists(f)) next
  sp <- suppressMessages(read_csv(f, show_col_types = FALSE))
  cat("\n>>> D2. SPATA2 p-value,", ds, "<<<\nColonne:", paste(names(sp), collapse = ", "), "\n")
  pc <- intersect(c("pval", "p_value", "pvalue", "p.value", "p"), names(sp))[1]
  if (!is.na(pc)) {
    v <- sp[[pc]]
    cat("p minimo:", signif(min(v, na.rm = TRUE), 3), "| geni al minimo:", sum(v == min(v, na.rm = TRUE), na.rm = TRUE),
        "| p NA:", sum(is.na(v)), "| p <= 1e-3:", sum(v <= 1e-3, na.rm = TRUE), "\n")
  } else cat("Nessuna colonna p-value riconosciuta fra le colonne sopra.\n")
  if ("fdr" %in% names(sp)) cat("fdr NA:", sum(is.na(sp$fdr)), "| fdr < ", ALPHA, ":", sum(sp$fdr < ALPHA, na.rm = TRUE), "\n")
}

# --- D3. Geni significativi per SPATA2 che non entrano nel confronto -------------------
# spata2_results.csv e' l'output completo di SPATA2; la pipeline tiene solo i geni del
# proprio universo (sezione 11b.4). Qui si vede cosa resta fuori fra i significativi.
for (ds in DATASETS) {
  f <- file.path(BASE_RESULTS, ds, "spata2", "spata2_results.csv")
  if (!file.exists(f)) next
  sp  <- suppressMessages(read_csv(f, show_col_types = FALSE))
  uni <- unique(long$gene_id[long$dataset == ds & long$Method != "SPATA2"])
  sp$sig <- !is.na(sp$fdr) & sp$fdr < ALPHA
  sp$in_universe <- !is.na(sp$gene_id) & sp$gene_id %in% uni
  sp$categoria <- ifelse(sp$in_universe, "nell'universo della pipeline",
                  ifelse(is.na(sp$gene_id), "non mappato a Ensembl",
                  ifelse(grepl("^mt-", sp$gene, ignore.case = TRUE), "mitocondriale",
                         "fuori universo (poco espresso o non annotato)")))
  cat("\n>>> D3. SPATA2, geni significativi dentro e fuori l'universo,", ds, "<<<\n")
  cat("Geni nell'output di SPATA2:", nrow(sp), "| universo della pipeline:", length(uni),
      "| significativi sul suo output:", sum(sp$sig),
      sprintf("(%.1f%%)", 100 * mean(sp$sig)), "\n")
  cat("Geni dell'universo assenti dall'output di SPATA2:", length(setdiff(uni, sp$gene_id)),
      "| gene_id duplicati nell'output:", sum(duplicated(sp$gene_id[!is.na(sp$gene_id)])), "\n")
  print(as.data.frame(sp %>% group_by(categoria) %>%
                        summarise(n_geni = n(), n_significativi = sum(sig),
                                  pct_significativi = round(100 * mean(sig), 1), .groups = "drop")),
        row.names = FALSE)
  fuori <- sp %>% filter(sig, !in_universe) %>% arrange(desc(tot_var)) %>%
    select(gene, gene_id, tot_var, p_value, fdr, categoria)
  write_csv(fuori, file.path(BASE_RESULTS, ds, "spata2", "spata2_sig_fuori_universo.csv"))
  cat("\nPrimi 25 significativi fuori dall'universo (per tot_var); elenco completo in spata2/spata2_sig_fuori_universo.csv:\n")
  print(as.data.frame(head(fuori, 25)), row.names = FALSE, digits = 3)
}

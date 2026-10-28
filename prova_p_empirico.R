# ##############################################################################
# ##   PROVA: P-VALUE EMPIRICO DAL NULLO SPAZIALE PER GLM_linear              ##
# ##                                                                          ##
# ##   Non rifa' la pipeline: carica il checkpoint della sezione 12b          ##
# ##   (checkpoint_pre_grafici.RData), che contiene gia' Z1_scan (score test  ##
# ##   lineare di ogni gene con ogni spot come centro, sezione 9c).           ##
# ##                                                                          ##
# ##   P empirico di un gene: frazione dei centri LONTANI dal barcode (oltre  ##
# ##   la mediana della distanza, come nella 10f) con Z almeno altrettanto    ##
# ##   estremo NELLA STESSA DIREZIONE del barcode (Z >= Z_obs se Z_obs > 0,  ##
# ##   Z <= Z_obs se Z_obs < 0). (k + 1) / (n + 1), poi BH sui geni.          ##
# ##   Ranking: fdr empirico, poi p empirico, poi z_std = (Z_obs - media dei  ##
# ##   centri lontani) / sd dei centri lontani, in valore assoluto.          ##
# ##   GLM_linear_emp_abs: versione precedente (|Z|, spareggio log p).        ##
# ##                                                                          ##
# ##   Confronto fra GLM_linear (p parametrico), GLM_linear_emp (p empirico), ##
# ##   SPATA2 e SpatialDE_orig: geni significativi, Golden Standard, GSEA     ##
# ##   come nella sezione 16 (solo geni testati, N..1, exponent = 0).         ##
# ##############################################################################

suppressPackageStartupMessages({
  library(SpatialExperiment); library(SingleCellExperiment); library(Matrix)
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
  library(edgeR); library(msigdbr); library(clusterProfiler)
})

# Checkpoint: primo argomento da riga di comando (lo usa prova_p_empirico_tutti.R),
# altrimenti CHECKPOINT_PATH gia' definito nella sessione, altrimenti heme_1000.
# Secondo argomento opzionale: lo script dei grafici, lanciato alla fine.
.prova_args <- commandArgs(trailingOnly = TRUE)
.prova_ck <- if (length(.prova_args) >= 1) .prova_args[1] else if (exists("CHECKPOINT_PATH")) CHECKPOINT_PATH else
  "/home/ammendola/laboratorio/spatial/topo/code/benchmark_v7/heme_1000/checkpoint_pre_grafici.RData"
load(.prova_ck)
CHECKPOINT_PATH <- .prova_ck
cat("\n========== PROVA P EMPIRICO:", SAMPLE_ID, "==========\n")
# Processori per questo run: PROVA_THREADS (lo imposta prova_p_empirico_tutti.R),
# altrimenti N_CORES del checkpoint. BLAS a 1 thread: i conti pesanti qui sono
# confronti fra matrici e GSEA, che si parallelizzano per metodo (sezione 5).
PROVA_THREADS <- suppressWarnings(as.integer(Sys.getenv("PROVA_THREADS", "")))
if (is.na(PROVA_THREADS) || PROVA_THREADS < 1) PROVA_THREADS <- if (exists("N_CORES")) N_CORES else 4L
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) RhpcBLASctl::blas_set_num_threads(1)
BiocParallel::register(BiocParallel::SerialParam())     # fgsea seriale dentro ogni processo
cat("Processori per questo run:", PROVA_THREADS, "\n")

METODI_PROVA <- c("GLM_linear", "GLM_linear_emp", "GLM_linear_emp_abs",
                  intersect(c("SPATA2", "SpatialDE_orig"), unique(as.character(all_results$Method))))
PROVA_SEED   <- 20261028
stopifnot(exists("Z1_scan"), exists("glm_lin"), exists("all_results"))

# --- 1. P empirico ---------------------------------------------------------------
d_bc    <- sqrt((spots$pixel_col - spots$pixel_col[center_idx])^2 +
                (spots$pixel_row - spots$pixel_row[center_idx])^2)
lontani <- which(d_bc > median(d_bc))
z_obs  <- Z1_scan[, center_idx]
Zl     <- Z1_scan[, lontani, drop = FALSE]
n_null <- rowSums(!is.na(Zl))
p_k    <- function(k) ifelse(is.na(z_obs) | n_null == 0, NA_real_, (k + 1) / (n_null + 1))

# con segno: solo i centri che danno lo stesso verso del barcode contano
k_sgn  <- ifelse(z_obs >= 0, rowSums(Zl >= z_obs, na.rm = TRUE), rowSums(Zl <= z_obs, na.rm = TRUE))
# versione precedente, in valore assoluto
k_abs  <- rowSums(abs(Zl) >= abs(z_obs), na.rm = TRUE)
# z standardizzato sui centri lontani dello stesso gene
mu_l   <- rowMeans(Zl, na.rm = TRUE)
sd_l   <- sqrt(rowSums((Zl - mu_l)^2, na.rm = TRUE) / pmax(n_null - 1, 1))
z_std  <- (z_obs - mu_l) / sd_l
rm(Zl); invisible(gc())

emp <- data.frame(gene_id = rownames(Z1_scan), z_barcode = z_obs, n_centri_null = n_null,
                  k_sgn = k_sgn, p_emp = p_k(k_sgn), k_abs = k_abs, p_emp_abs = p_k(k_abs),
                  media_lontani = mu_l, sd_lontani = sd_l, z_std = z_std, stringsAsFactors = FALSE)
emp$fdr_emp     <- p.adjust(emp$p_emp, "BH")
emp$fdr_emp_abs <- p.adjust(emp$p_emp_abs, "BH")
emp <- emp %>% left_join(glm_lin %>% select(gene_id, p_value, fdr, log_p), by = "gene_id")
rk <- function(...) { r <- order(order(...)); r }
emp$rank     <- rk(emp$fdr_emp, emp$p_emp, -abs(emp$z_std))
emp$rank_abs <- rk(emp$fdr_emp_abs, emp$p_emp_abs, emp$log_p)
emp$rank[is.na(emp$fdr_emp)]         <- NA_integer_
emp$rank_abs[is.na(emp$fdr_emp_abs)] <- NA_integer_
write_csv(emp, file.path(OUTPUT_DIR, "prova_p_empirico_glm_linear.csv"))

cat("Centri lontani usati:", length(lontani), "| p empirico minimo possibile:",
    signif(1 / (median(n_null) + 1), 3), "\n")
cat("Geni con p empirico:", sum(!is.na(emp$p_emp)), "| al minimo: con segno", sum(emp$k_sgn == 0, na.rm = TRUE),
    "| in valore assoluto", sum(emp$k_abs == 0, na.rm = TRUE), "\n")
cat("GLM_linear significativi: parametrico", sum(emp$fdr < ALPHA_FDR, na.rm = TRUE),
    "| empirico con segno", sum(emp$fdr_emp < ALPHA_FDR, na.rm = TRUE),
    "| empirico |Z|", sum(emp$fdr_emp_abs < ALPHA_FDR, na.rm = TRUE), "\n")
cat("Spearman fra z_std e log p del GLM (quanto lo spareggio nuovo segue ancora la potenza):",
    round(cor(abs(emp$z_std), -emp$log_p, method = "spearman", use = "complete.obs"), 3), "\n")

gs_ids <- unique(all_results$gene_id[all_results$is_GS %in% TRUE])
prova <- bind_rows(
  all_results %>% filter(Method %in% setdiff(METODI_PROVA, "GLM_linear_emp")) %>%
    mutate(Method = as.character(Method)) %>%
    select(Method, gene_id, fdr, rank, tested),
  emp %>% transmute(Method = "GLM_linear_emp", gene_id, fdr = fdr_emp, rank, tested = !is.na(fdr_emp)),
  emp %>% transmute(Method = "GLM_linear_emp_abs", gene_id, fdr = fdr_emp_abs, rank = rank_abs,
                    tested = !is.na(fdr_emp_abs))) %>%
  mutate(tested = dplyr::coalesce(tested, TRUE), is_GS = gene_id %in% gs_ids,
         sig = !is.na(fdr) & fdr < ALPHA_FDR)

# --- 2. Geni significativi e Golden Standard ------------------------------------------
auc_gs <- function(rank, gs) {
  ok <- !is.na(rank); rank <- rank[ok]; gs <- gs[ok]; n1 <- sum(gs); n0 <- sum(!gs)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(rank(-rank)[gs]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
tab_geni <- prova %>% filter(tested) %>% group_by(Method) %>%
  summarise(n_tested = n(), n_sig = sum(sig), pct_sig = round(100 * mean(sig), 1),
            GS_recall = round(100 * sum(sig & is_GS) / max(sum(is_GS), 1), 1),
            GS_fold = round((sum(sig & is_GS) / max(sum(sig), 1)) / (sum(is_GS) / n()), 2),
            AUC_GS = round(auc_gs(rank, is_GS), 3), .groups = "drop")
cat("\n>>> Geni significativi e Golden Standard <<<\n")
print(as.data.frame(tab_geni), row.names = FALSE)

# --- 3. Sovrapposizione con SPATA2 ---------------------------------------------------
tab_sp <- data.frame()
if ("SPATA2" %in% METODI_PROVA) {
sp_ref <- prova %>% filter(Method == "SPATA2", tested)
sig_sp <- sp_ref$gene_id[sp_ref$sig]
tab_sp <- bind_rows(lapply(setdiff(METODI_PROVA, "SPATA2"), function(m) {
  o <- prova %>% filter(Method == m, tested); com <- intersect(o$gene_id, sp_ref$gene_id)
  data.frame(Method = m, n_sig = sum(o$sig),
             pct_SPATA2_anche_qui = round(100 * length(intersect(sig_sp, o$gene_id[o$sig])) / max(length(sig_sp), 1), 1),
             pct_atteso_a_caso = round(100 * mean(o$sig), 1),
             spearman_con_SPATA2 = round(cor(sp_ref$rank[match(com, sp_ref$gene_id)],
                                             o$rank[match(com, o$gene_id)], method = "spearman"), 3))
}))
cat("\n>>> Sovrapposizione con SPATA2 <<<\n")
print(tab_sp, row.names = FALSE)
}
cat("Correlazione di Spearman fra il rank parametrico e quello empirico di GLM_linear:",
    round(cor(emp$rank, order(order(emp$fdr, emp$log_p)), method = "spearman", use = "complete.obs"), 3), "\n")

# --- 4. Pathway (stessa costruzione della sezione 15) ----------------------------------
 
 
PATHWAY_SETS <- list(
  tumor = list(
    hallmark = c("HALLMARK_HYPOXIA", "HALLMARK_GLYCOLYSIS", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                "HALLMARK_ANGIOGENESIS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_E2F_TARGETS",
                "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MTORC1_SIGNALING", "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
                "HALLMARK_P53_PATHWAY", "HALLMARK_APOPTOSIS", "HALLMARK_KRAS_SIGNALING_UP",
                "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INFLAMMATORY_RESPONSE",
                "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_COMPLEMENT"),
    bp_regex = paste0("HYPOXIA|ANGIOGEN|VASCUL|CELL_CYCLE|MITOTIC|DNA_REPLICATION|PROLIFERATION|APOPTO|",
                      "CELL_MIGRATION|EXTRACELLULAR_MATRIX|WOUND_HEALING|IMMUNE|INFLAMMAT|CYTOKINE|",
                      "INTERFERON|CHEMOKINE|LEUKOCYTE|MYELOID")
  ),
  inflammation = list(
    hallmark = c("HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_COMPLEMENT", "HALLMARK_IL2_STAT5_SIGNALING"),
    bp_regex = paste0("INFLAMMAT|CYTOKINE|INTERLEUKIN|CHEMOKINE|NF_KAPPA|TOLL_LIKE|INTERFERON|",
                      "IMMUNE_RESPONSE|LEUKOCYTE|MICROGLIA")
  )
)
pathway_cfg <- PATHWAY_SETS[[ENRICHMENT_CONTEXT]]
 
get_msig <- function(coll, sub = NULL) {
  a <- names(formals(msigdbr::msigdbr))
  if ("collection" %in% a) {
    msigdbr(species = ORGANISM, collection = coll, subcollection = sub)
  } else {
    msigdbr(species = ORGANISM, category = coll, subcategory = sub)
  }
}
msig_h <- get_msig("H"); msig_bp <- get_msig("C5", "GO:BP")
ens_col <- intersect(c("ensembl_gene", "db_ensembl_gene"), names(msig_h))[1]
 
term2gene <- bind_rows(msig_h %>% select(gs_name, gene = all_of(ens_col)),
                       msig_bp %>% select(gs_name, gene = all_of(ens_col))) %>%
  filter(!is.na(gene)) %>% distinct()
 
all_terms <- unique(term2gene$gs_name)
TARGET_PATHWAYS <- unique(c(intersect(pathway_cfg$hallmark, all_terms),
                           grep(pathway_cfg$bp_regex, grep("^GOBP", all_terms, value = TRUE), value = TRUE)))
cat("Pathway di contesto (", ENRICHMENT_CONTEXT, "):", length(TARGET_PATHWAYS), "su", length(all_terms), "testabili\n")
 
clean_label <- function(id) gsub("_", " ", gsub("^(HALLMARK|GOBP)_", "", id))
 
# --- 15.1  Universo e pathway non ridondanti ----------------------------------
# Universo: geni della tabella unificata. Pathway con 10-500 geni nell'universo.
# Ridondanza: un pathway entra solo se il suo Jaccard con ognuno di quelli gia'
# tenuti e' < 0.5 (prima gli Hallmark, poi i piu' grandi). Ogni GSEA (sezione 16)
# lavora poi solo sui geni testati dal metodo.
U_enr <- sort(unique(as.character(all_results$gene_id)))
stopifnot(all(U_enr %in% rownames(logcounts_raw)))
expr_U <- rowMeans(logcounts_raw)[U_enr]
 
t2g_U <- term2gene %>% filter(gene %in% U_enr)
pn    <- unique(t2g_U$gs_name)
M_pw  <- sparseMatrix(i = match(t2g_U$gs_name, pn), j = match(t2g_U$gene, U_enr), x = 1,
                      dims = c(length(pn), length(U_enr)), dimnames = list(pn, U_enr))
M_pw  <- M_pw[rowSums(M_pw) >= 10 & rowSums(M_pw) <= 500, ]
rs <- rowSums(M_pw); Mt <- t(M_pw); keep <- integer(0)
for (i in order(!grepl("^HALLMARK", rownames(M_pw)), -rs)) {
  v <- as.numeric(M_pw %*% Mt[, i])
  if (length(keep) == 0 || all(v[keep] / (rs[keep] + rs[i] - v[keep]) < 0.5)) keep <- c(keep, i)
}
M_pw <- M_pw[keep, ]
rm(Mt); invisible(gc())
 
T2_pw <- t2g_U %>% filter(gs_name %in% rownames(M_pw))
cat("Pathway non ridondanti:", nrow(M_pw), "| di contesto:", sum(rownames(M_pw) %in% TARGET_PATHWAYS),
    "| geni nell'universo:", length(U_enr), "\n")
stopifnot(sum(rownames(M_pw) %in% TARGET_PATHWAYS) >= 5,
          sum(!rownames(M_pw) %in% TARGET_PATHWAYS) >= 5)

# --- 5. GSEA come nella sezione 16 -----------------------------------------------------
pos_score <- function(ids) setNames(rev(seq_along(ids)), ids)
SC <- lapply(setNames(METODI_PROVA, METODI_PROVA), function(m) {
  d <- prova %>% filter(Method == m, tested, !is.na(rank)) %>% arrange(rank)
  pos_score(d$gene_id)
})
SC$Expression_control <- pos_score(U_enr[order(-expr_U)])

# una GSEA per processo figlio, fino a PROVA_THREADS insieme
gsea_prova <- bind_rows(parallel::mclapply(names(SC), mc.cores = min(PROVA_THREADS, length(SC)), function(m) {
  cat("  GSEA:", m, "(", length(SC[[m]]), "geni )\n")
  set.seed(PROVA_SEED)
  g <- as.data.frame(suppressWarnings(GSEA(SC[[m]], TERM2GENE = T2_pw, minGSSize = 10, maxGSSize = 500,
                                           pvalueCutoff = 1, pAdjustMethod = "BH", scoreType = "pos", exponent = 0,
                                           verbose = FALSE, seed = TRUE)))
  data.frame(Method = m, pathway = g$ID, NES = g$NES, p_adjust = g$p.adjust,
             is_context = g$ID %in% TARGET_PATHWAYS, stringsAsFactors = FALSE)
}))
stopifnot(all(names(SC) %in% gsea_prova$Method))     # un figlio fallito non deve passare inosservato
write_csv(gsea_prova, file.path(OUTPUT_DIR, "prova_p_empirico_gsea.csv"))

auc_ctx <- function(nes, ctx) {
  ok <- !is.na(nes); nes <- nes[ok]; ctx <- ctx[ok]; n1 <- sum(ctx); n0 <- sum(!ctx)
  if (n1 < 5 || n0 < 5) return(NA_real_)
  (sum(rank(nes)[ctx]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
nes_expr <- gsea_prova %>% filter(Method == "Expression_control") %>% select(pathway, NES_expr = NES)
tab_gsea <- gsea_prova %>% left_join(nes_expr, by = "pathway") %>% group_by(Method) %>%
  summarise(auc_ctx_nes = round(auc_ctx(NES, is_context), 3),
            nes_ctx_med = round(median(NES[is_context], na.rm = TRUE), 2),
            nes_other_med = round(median(NES[!is_context], na.rm = TRUE), 2),
            cor_con_espressione = round(cor(NES, NES_expr, method = "spearman", use = "pairwise.complete.obs"), 2),
            n_p05 = sum(p_adjust < ALPHA_FDR, na.rm = TRUE),
            n_p05_ctx = sum(p_adjust < ALPHA_FDR & is_context, na.rm = TRUE), .groups = "drop")
cat("\n>>> GSEA (n_p05 solo diagnostica) <<<\n")
print(as.data.frame(tab_gsea), row.names = FALSE)

# --- 6. Tabelle su disco (le legge prova_p_empirico_tutti.R) ---------------------------
write_csv(tab_geni %>% mutate(dataset = SAMPLE_ID), file.path(OUTPUT_DIR, "prova_p_empirico_tab_geni.csv"))
write_csv(tab_gsea %>% mutate(dataset = SAMPLE_ID), file.path(OUTPUT_DIR, "prova_p_empirico_tab_gsea.csv"))
if (nrow(tab_sp) > 0) write_csv(tab_sp %>% mutate(dataset = SAMPLE_ID), file.path(OUTPUT_DIR, "prova_p_empirico_tab_spata2.csv"))
write_csv(data.frame(dataset = SAMPLE_ID, centri_lontani = length(lontani),
                     n_sig_param = sum(emp$fdr < ALPHA_FDR, na.rm = TRUE),
                     n_sig_emp = sum(emp$fdr_emp < ALPHA_FDR, na.rm = TRUE),
                     n_sig_emp_abs = sum(emp$fdr_emp_abs < ALPHA_FDR, na.rm = TRUE),
                     al_minimo_emp = sum(emp$k_sgn == 0, na.rm = TRUE),
                     al_minimo_emp_abs = sum(emp$k_abs == 0, na.rm = TRUE),
                     spearman_zstd_logp = round(cor(abs(emp$z_std), -emp$log_p, method = "spearman", use = "complete.obs"), 3)),
          file.path(OUTPUT_DIR, "prova_p_empirico_tab_riassunto.csv"))

if (length(.prova_args) >= 2 && file.exists(.prova_args[2])) source(.prova_args[2])

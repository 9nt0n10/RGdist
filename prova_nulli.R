# ##############################################################################
# ##   PROVA: DUE NULLI SPAZIALI PER IL GLM LINEARE                           ##
# ##                                                                          ##
# ##   Riparte dal checkpoint della sezione 12b (nessun fit dei GLM rifatto). ##
# ##   Per ogni gene lo Z dello score test lineare (sezione 9c) al centro     ##
# ##   vero viene confrontato con due distribuzioni nulle:                    ##
# ##                                                                          ##
# ##   corr     : Z degli altri spot presi come centro, tenendo solo quelli   ##
# ##              la cui mappa di distanza ha |r| < SOGLIA_COR con quella     ##
# ##              vera (ne' stesso gradiente ne' gradiente opposto).          ##
# ##   spettrale: Z di B_SPET covariate casuali con lo stesso spettro di      ##
# ##              autocorrelazione del raggio (randomizzazione spettrale      ##
# ##              sugli autovettori della matrice dei vicini, kNN).           ##
# ##                                                                          ##
# ##   p = (k + 1) / (n + 1) su |Z|, BH sui geni. Ranking: fdr, p, log p del  ##
# ##   GLM (sezione 7), come la regola della pipeline.                        ##
# ##                                                                          ##
# ##   Controllo di validita' (falsi centri): gli stessi test applicati a     ##
# ##   N_FALSI centri finti, poco correlati con il vero. Un test calibrato    ##
# ##   da' circa il 5% di geni con p < 0.05 e quasi nessun BH-significativo.  ##
# ##                                                                          ##
# ##   Poi Golden Standard e GSEA come la sezione 16 (solo geni testati,      ##
# ##   punteggio N..1, exponent = 0), con SPATA2, SpatialDE e il controllo    ##
# ##   sull'espressione. Tutto in CSV nella cartella del checkpoint.          ##
# ##############################################################################

suppressPackageStartupMessages({
  library(SpatialExperiment); library(SingleCellExperiment); library(Matrix)
  library(readr); library(dplyr); library(tidyr); library(edgeR)
  library(msigdbr); library(clusterProfiler); library(FNN)
})

.nl_args <- commandArgs(trailingOnly = TRUE)
.nl_ck <- if (length(.nl_args) >= 1) .nl_args[1] else if (exists("CHECKPOINT_PATH")) CHECKPOINT_PATH else
  "/home/ammendola/laboratorio/spatial/topo/code/benchmark_v7/heme_1000/checkpoint_pre_grafici.RData"
load(.nl_ck)
CHECKPOINT_PATH <- .nl_ck
cat("\n========== PROVA NULLI:", SAMPLE_ID, "==========\n")
stopifnot(exists("Z1_scan"), exists("glm_lin"), exists("all_results"), exists("spe"))

NL_THREADS <- suppressWarnings(as.integer(Sys.getenv("PROVA_THREADS", "")))
if (is.na(NL_THREADS) || NL_THREADS < 1) NL_THREADS <- if (exists("N_CORES")) N_CORES else 4L
blas <- function(k) if (requireNamespace("RhpcBLASctl", quietly = TRUE)) RhpcBLASctl::blas_set_num_threads(k)
blas(NL_THREADS)
BiocParallel::register(BiocParallel::SerialParam())
cat("Processori:", NL_THREADS, "\n")

# --- parametri -----------------------------------------------------------------------
SOGLIA_COR   <- 0.3      # corr: |r| massimo fra la mappa di distanza del candidato e quella vera
N_MIN_CENTRI <- 50       # corr: sotto questo numero di centri il test non si fa
K_VICINI     <- 6        # spettrale: vicini per la matrice di adiacenza
BANDA        <- 10       # spettrale: mezza larghezza (in autovettori) per smussare lo spettro
B_SPET       <- 1000     # spettrale: campi casuali al centro vero
N_FALSI      <- 10       # falsi centri
B_FALSI      <- 200      # spettrale: campi casuali per ogni falso centro
NL_SEED      <- 20261029
ALPHA        <- if (exists("ALPHA_FDR")) ALPHA_FDR else 0.05

t_tot <- Sys.time()
xy <- cbind(spots$pixel_col, spots$pixel_row); n <- nrow(xy)
D  <- as.matrix(dist(xy))
stopifnot(identical(colnames(Z1_scan), colnames(spe)), ncol(Z1_scan) == n)
geni   <- rownames(Z1_scan)
log_p  <- glm_lin$log_p[match(geni, glm_lin$gene_id)]
z_true <- Z1_scan[, center_idx]

# --- residui e pesi del fit nullo (come la 9c) ---------------------------------------
if (!exists("dge_null")) dge_null <- estimateDisp(dge, matrix(1, n, 1))
phi <- if (!is.null(dge_null$tagwise.dispersion)) dge_null$tagwise.dispersion else
  rep(dge_null$common.dispersion, nrow(dge_null))
mu  <- glmFit(dge_null, matrix(1, n, 1))$fitted.values
Y   <- as.matrix(counts(spe))
den <- 1 + phi * mu
R   <- (Y - mu) / den; W <- mu / den; S0 <- rowSums(W)
rm(Y, mu, den); invisible(gc())
stopifnot(identical(rownames(R), geni))

# p empirico su |Z| per una matrice di covariate X [spot x campi], a blocchi
p_campi <- function(X, z_ref, blocco = 200) {
  k <- numeric(length(z_ref)); nv <- numeric(length(z_ref))
  for (s in seq(1, ncol(X), by = blocco)) {
    Xb <- X[, s:min(s + blocco - 1, ncol(X)), drop = FALSE]
    S1 <- W %*% Xb; V <- W %*% (Xb^2) - S1^2 / S0; V[V <= 0] <- NA
    Zb <- abs((R %*% Xb) / sqrt(V))
    k  <- k + rowSums(Zb >= abs(z_ref), na.rm = TRUE); nv <- nv + rowSums(!is.na(Zb))
  }
  ifelse(is.na(z_ref) | nv == 0, NA_real_, (k + 1) / (nv + 1))
}
# controllo: lo score test ricalcolato al centro vero coincide con Z1_scan
z_chk <- { x <- D[, center_idx] / max(D[, center_idx]); S1 <- W %*% x
           as.numeric((R %*% x) / sqrt(W %*% x^2 - S1^2 / S0)) }
cat("Controllo score test (max |differenza| con Z1_scan):", signif(max(abs(z_chk - z_true), na.rm = TRUE), 3), "\n")

# --- nullo 1: centri poco correlati ---------------------------------------------------
centri_corr <- function(j) { r <- cor(D[, j], D)[1, ]; which(abs(r) < SOGLIA_COR & seq_len(n) != j) }
p_corr <- function(j) {
  cc <- centri_corr(j)
  if (length(cc) < N_MIN_CENTRI) return(list(p = rep(NA_real_, length(geni)), n = length(cc)))
  Zc <- abs(Z1_scan[, cc, drop = FALSE]); zj <- abs(Z1_scan[, j])
  nv <- rowSums(!is.na(Zc))
  list(p = ifelse(is.na(zj) | nv == 0, NA_real_, (rowSums(Zc >= zj, na.rm = TRUE) + 1) / (nv + 1)), n = length(cc))
}

# --- nullo 2: randomizzazione spettrale ---------------------------------------------------
knn <- FNN::get.knn(xy, k = K_VICINI)$nn.index
A   <- matrix(0, n, n); A[cbind(rep(seq_len(n), K_VICINI), as.vector(knn))] <- 1
A   <- pmax(A, t(A))
Ac  <- sweep(A, 1, rowMeans(A)); Ac <- sweep(Ac, 2, colMeans(Ac))
cat("Autovettori della matrice dei vicini (", n, " spot) ...\n", sep = "")
V_ev <- eigen(Ac, symmetric = TRUE)$vectors
rm(Ac); invisible(gc())
moran <- function(x) { xc <- x - mean(x); (n / sum(A)) * sum(xc * (A %*% xc)) / sum(xc^2) }
campi_spettrali <- function(x, B, seed) {
  pw  <- crossprod(V_ev, x - mean(x))[, 1]^2
  pws <- vapply(seq_len(n), function(i) mean(pw[max(1, i - BANDA):min(n, i + BANDA)]), numeric(1))
  set.seed(seed)
  V_ev %*% (matrix(rnorm(n * B), n, B) * sqrt(pws))
}

# --- centro vero ------------------------------------------------------------------------
x_true <- D[, center_idx]
pc <- p_corr(center_idx)
cat("Nullo corr: centri con |r| <", SOGLIA_COR, ":", pc$n, if (pc$n < N_MIN_CENTRI) "(troppo pochi: saltato)" else "", "\n")
X_true <- campi_spettrali(x_true, B_SPET, NL_SEED)
mor_campi <- apply(X_true[, seq_len(min(100, B_SPET)), drop = FALSE], 2, moran)
mor_raggio <- moran(x_true)
cat("Moran's I del raggio:", round(mor_raggio, 3), "| dei campi casuali: media", round(mean(mor_campi), 3),
    "sd", round(sd(mor_campi), 3), "\n")
p_sp <- p_campi(X_true, z_true)
rm(X_true); invisible(gc())

nul <- data.frame(gene_id = geni, z = z_true, p_score = 2 * pnorm(-abs(z_true)),
                  p_corr = pc$p, p_spet = p_sp, log_p = log_p, stringsAsFactors = FALSE) %>%
  left_join(glm_lin %>% select(gene_id, p_QL = p_value, fdr_QL = fdr), by = "gene_id") %>%
  mutate(fdr_corr = p.adjust(p_corr, "BH"), fdr_spet = p.adjust(p_spet, "BH"))
rk <- function(f, p) { r <- order(order(f, p, nul$log_p)); r[is.na(f)] <- NA_integer_; r }
nul$rank_corr <- rk(nul$fdr_corr, nul$p_corr)
nul$rank_spet <- rk(nul$fdr_spet, nul$p_spet)
write_csv(nul, file.path(OUTPUT_DIR, "prova_nulli_pvalori.csv"))
cat("Significativi GLM_linear: QL", sum(nul$fdr_QL < ALPHA, na.rm = TRUE),
    "| corr", sum(nul$fdr_corr < ALPHA, na.rm = TRUE), "| spettrale", sum(nul$fdr_spet < ALPHA, na.rm = TRUE), "\n")
cat("Geni al p minimo: corr", sum(nul$p_corr == min(nul$p_corr, na.rm = TRUE), na.rm = TRUE),
    "| spettrale", sum(nul$p_spet == min(nul$p_spet, na.rm = TRUE), na.rm = TRUE), "\n")

# --- falsi centri -------------------------------------------------------------------------
set.seed(NL_SEED)
cand_f <- centri_corr(center_idx)
if (length(cand_f) < N_FALSI) cand_f <- which(D[, center_idx] > median(D[, center_idx]))
falsi  <- cand_f[sample.int(length(cand_f), min(N_FALSI, length(cand_f)))]
riga_f <- function(j, metodo, p) data.frame(dataset = SAMPLE_ID, falso = colnames(spe)[j], Method = metodo,
                                            frac_p05 = mean(p < 0.05, na.rm = TRUE),
                                            n_bh = sum(p.adjust(p, "BH") < ALPHA, na.rm = TRUE))
tab_falsi <- bind_rows(lapply(seq_along(falsi), function(i) {
  j <- falsi[i]; cat("  falso centro", i, "/", length(falsi), "\n")
  pcj <- p_corr(j)
  psj <- p_campi(campi_spettrali(D[, j], B_FALSI, NL_SEED + i), Z1_scan[, j])
  bind_rows(riga_f(j, "GLM_linear", 2 * pnorm(-abs(Z1_scan[, j]))),
            if (pcj$n >= N_MIN_CENTRI) riga_f(j, "GLM_linear_corr", pcj$p),
            riga_f(j, "GLM_linear_spettrale", psj))
}))
cat("\n>>> Falsi centri: un test calibrato da' frac_p05 vicino a 0.05 e n_bh vicino a 0 <<<\n")
print(as.data.frame(tab_falsi %>% group_by(Method) %>%
                      summarise(frac_p05_media = round(mean(frac_p05), 3), frac_p05_max = round(max(frac_p05), 3),
                                n_bh_media = round(mean(n_bh), 1), n_bh_max = max(n_bh), .groups = "drop")),
      row.names = FALSE)
cat("(GLM_linear qui e' lo score test parametrico della 9c, come nella 10f)\n")
rm(R, W, V_ev, A, D); invisible(gc())
blas(1)

# --- Golden Standard -------------------------------------------------------------------------
METODI <- c("GLM_linear", "GLM_linear_corr", "GLM_linear_spettrale",
            intersect(c("SPATA2", "SpatialDE_orig"), unique(as.character(all_results$Method))))
if (pc$n < N_MIN_CENTRI) METODI <- setdiff(METODI, "GLM_linear_corr")
gs_ids <- unique(all_results$gene_id[all_results$is_GS %in% TRUE])
prova <- bind_rows(
  all_results %>% filter(Method %in% c("GLM_linear", "SPATA2", "SpatialDE_orig")) %>%
    mutate(Method = as.character(Method)) %>% select(Method, gene_id, fdr, rank, tested),
  nul %>% transmute(Method = "GLM_linear_corr", gene_id, fdr = fdr_corr, rank = rank_corr, tested = !is.na(fdr_corr)),
  nul %>% transmute(Method = "GLM_linear_spettrale", gene_id, fdr = fdr_spet, rank = rank_spet, tested = !is.na(fdr_spet)),
  data.frame(Method = "Expression_control", gene_id = geni, fdr = NA_real_,
             rank = order(order(-rowMeans(logcounts_raw)[geni])), tested = TRUE)) %>%
  filter(Method %in% c(METODI, "Expression_control")) %>%
  mutate(tested = dplyr::coalesce(tested, TRUE), is_GS = gene_id %in% gs_ids, sig = !is.na(fdr) & fdr < ALPHA)

auc_gs <- function(rank, gs) {
  ok <- !is.na(rank); rank <- rank[ok]; gs <- gs[ok]; n1 <- sum(gs); n0 <- sum(!gs)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(rank(-rank)[gs]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
tab_geni <- prova %>% filter(tested) %>% group_by(Method) %>%
  summarise(n_tested = n(), n_sig = sum(sig), pct_sig = round(100 * mean(sig), 1),
            GS_recall = round(100 * sum(sig & is_GS) / max(sum(is_GS), 1), 1),
            GS_fold = round((sum(sig & is_GS) / max(sum(sig), 1)) / max(sum(is_GS) / n(), 1e-12), 2),
            AUC_GS = round(auc_gs(rank, is_GS), 3), .groups = "drop")
cat("\n>>> Geni significativi e Golden Standard <<<\n")
print(as.data.frame(tab_geni), row.names = FALSE)

# --- Pathway (stessa costruzione della sezione 15) -------------------------------------------
 
 
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

# --- GSEA come nella sezione 16 ------------------------------------------------------------------
pos_score <- function(ids) setNames(rev(seq_along(ids)), ids)
SC <- lapply(setNames(METODI, METODI), function(m) {
  d <- prova %>% filter(Method == m, tested, !is.na(rank)) %>% arrange(rank); pos_score(d$gene_id)
})
SC$Expression_control <- pos_score(U_enr[order(-expr_U)])
gsea_nul <- bind_rows(parallel::mclapply(names(SC), mc.cores = min(NL_THREADS, length(SC)), function(m) {
  set.seed(NL_SEED)
  g <- as.data.frame(suppressWarnings(GSEA(SC[[m]], TERM2GENE = T2_pw, minGSSize = 10, maxGSSize = 500,
                                           pvalueCutoff = 1, pAdjustMethod = "BH", scoreType = "pos", exponent = 0,
                                           verbose = FALSE, seed = TRUE)))
  data.frame(Method = m, pathway = g$ID, NES = g$NES, p_adjust = g$p.adjust,
             is_context = g$ID %in% TARGET_PATHWAYS, stringsAsFactors = FALSE)
}))
stopifnot(all(names(SC) %in% gsea_nul$Method))

auc_ctx <- function(nes, ctx) {
  ok <- !is.na(nes); nes <- nes[ok]; ctx <- ctx[ok]; n1 <- sum(ctx); n0 <- sum(!ctx)
  if (n1 < 5 || n0 < 5) return(NA_real_)
  (sum(rank(nes)[ctx]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
nes_expr <- gsea_nul %>% filter(Method == "Expression_control") %>% select(pathway, NES_expr = NES)
tab_gsea <- gsea_nul %>% left_join(nes_expr, by = "pathway") %>% group_by(Method) %>%
  summarise(auc_ctx_nes = round(auc_ctx(NES, is_context), 3),
            nes_ctx_med = round(median(NES[is_context], na.rm = TRUE), 2),
            nes_other_med = round(median(NES[!is_context], na.rm = TRUE), 2),
            cor_con_espressione = round(cor(NES, NES_expr, method = "spearman", use = "pairwise.complete.obs"), 2),
            n_p05 = sum(p_adjust < ALPHA, na.rm = TRUE),
            n_p05_ctx = sum(p_adjust < ALPHA & is_context, na.rm = TRUE), .groups = "drop")
cat("\n>>> GSEA (n_p05 solo diagnostica) <<<\n")
print(as.data.frame(tab_gsea), row.names = FALSE)

# --- su disco ---------------------------------------------------------------------------------------
write_csv(prova %>% select(Method, gene_id, fdr, rank, tested, is_GS, sig) %>% mutate(dataset = SAMPLE_ID),
          file.path(OUTPUT_DIR, "prova_nulli_ranks.csv"))
write_csv(gsea_nul %>% mutate(dataset = SAMPLE_ID), file.path(OUTPUT_DIR, "prova_nulli_gsea.csv"))
write_csv(tab_geni %>% mutate(dataset = SAMPLE_ID), file.path(OUTPUT_DIR, "prova_nulli_tab_geni.csv"))
write_csv(tab_gsea %>% mutate(dataset = SAMPLE_ID), file.path(OUTPUT_DIR, "prova_nulli_tab_gsea.csv"))
write_csv(tab_falsi, file.path(OUTPUT_DIR, "prova_nulli_falsi_centri.csv"))
write_csv(nul %>% select(gene_id, p_QL, p_score, p_corr, p_spet) %>% mutate(dataset = SAMPLE_ID),
          file.path(OUTPUT_DIR, "prova_nulli_pvalori_long.csv"))
write_csv(data.frame(dataset = SAMPLE_ID, n_spot = n, centri_corr = pc$n,
                     moran_raggio = round(mor_raggio, 3), moran_campi = round(mean(mor_campi), 3),
                     n_sig_QL = sum(nul$fdr_QL < ALPHA, na.rm = TRUE),
                     n_sig_corr = sum(nul$fdr_corr < ALPHA, na.rm = TRUE),
                     n_sig_spet = sum(nul$fdr_spet < ALPHA, na.rm = TRUE)),
          file.path(OUTPUT_DIR, "prova_nulli_riassunto.csv"))
cat("\nFatto in", round(as.numeric(difftime(Sys.time(), t_tot, units = "mins")), 1), "min\n")

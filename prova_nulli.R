# ##############################################################################
# ##   PROVA DI DUE NULLI SPAZIALI PER IL GLM LINEARE - FILE UNICO            ##
# ##                                                                          ##
# ##   Si lancia da solo (source o Rscript): carica i checkpoint della        ##
# ##   sezione 12b dei dataset, li analizza in parallelo e scrive tabelle e   ##
# ##   grafici in una sola cartella (OUT_DIR). Nessun altro script, niente    ##
# ##   scritto nelle cartelle dei dataset.                                    ##
# ##                                                                          ##
# ##   Per ogni gene lo Z dello score test lineare (sezione 9c) al centro     ##
# ##   vero si confronta con due distribuzioni nulle:                         ##
# ##     corr     : Z degli altri spot presi come centro, solo quelli la cui  ##
# ##                mappa di distanza ha |r| < SOGLIA_COR con quella vera.    ##
# ##     spettrale: Z di B_SPET covariate casuali con lo stesso spettro di    ##
# ##                autocorrelazione del raggio (autovettori della matrice    ##
# ##                dei vicini).                                              ##
# ##   p = (k + 1) / (n + 1) su |Z|, BH sui geni. Ranking: fdr, p, log p del  ##
# ##   GLM (sezione 7). GLM_linear resta il p parametrico di edgeR.           ##
# ##     spettrale_sigma: dagli stessi campi si stima, gene per gene, la sd   ##
# ##                di Z sotto il nullo (sigma) e p = 2 * pnorm(-|Z| / sigma).##
# ##                Continuo: nessun limite di risoluzione.                   ##
# ##     spettrale_10k: come spettrale ma con B_VERIFICA campi, solo per i    ##
# ##                dataset in DATASET_VERIFICA (verifica della risoluzione). ##
# ##                                                                          ##
# ##   Controllo (falsi centri): gli stessi test su N_FALSI centri finti. Un  ##
# ##   test affidabile da' circa il 5% di geni con p < 0.05 e quasi nessun    ##
# ##   gene BH-significativo. Non corregge niente: verifica i test.          ##
# ##                                                                          ##
# ##   Poi Golden Standard (boxplot, ROC, AUC con IC) e GSEA come la sezione  ##
# ##   16, con SPATA2, SpatialDE e il controllo sull'espressione.            ##
# ##############################################################################

suppressPackageStartupMessages({
  library(SpatialExperiment); library(SingleCellExperiment); library(Matrix)
  library(readr); library(dplyr); library(tidyr); library(ggplot2); library(edgeR)
  library(msigdbr); library(clusterProfiler); library(FNN)
})

# --- da impostare ----------------------------------------------------------------------
OUT_DIR  <- "/home/ammendola/laboratorio/spatial/test_arricchimenti/prova_nulli"
N_CPU    <- 40
MAX_JOBS <- 6        # dataset in parallelo: abbassalo se la RAM non basta (1 = uno alla volta)
CHECKPOINTS <- c(
  "/home/ammendola/laboratorio/spatial/topo/code/benchmark_v7/heme_1000/checkpoint_pre_grafici.RData",
  "/home/ammendola/laboratorio/spatial/brain_injury/code/stab_inj1/benchmark_v7/stab_inj1/checkpoint_pre_grafici.RData",
  "/home/ammendola/laboratorio/spatial/brain_injury/code/stab_inj2/benchmark_v7/stab_inj2/checkpoint_pre_grafici.RData",
  "/home/ammendola/laboratorio/spatial/glioblastoma_mice/code/benchmark_v7/gbm/checkpoint_pre_grafici.RData",
  "/home/ammendola/laboratorio/spatial/periodontitis/code/benchmark_v7/perio/checkpoint_pre_grafici.RData",
  "/home/ammendola/laboratorio/spatial/tumor/full_tissue/code/benchmark_v7/head_and_neck/checkpoint_pre_grafici.RData"
)
SOGLIA_COR   <- 0.3      # corr: |r| massimo fra la mappa di distanza del candidato e quella vera
N_MIN_CENTRI <- 50       # corr: sotto questo numero di centri il test non si fa
K_VICINI     <- 6        # spettrale: vicini per la matrice di adiacenza
BANDA        <- 10       # spettrale: mezza larghezza (in autovettori) per smussare lo spettro
B_SPET       <- 1000     # spettrale: campi casuali al centro vero
N_FALSI      <- 10       # falsi centri
B_FALSI      <- 200      # spettrale: campi casuali per ogni falso centro
B_VERIFICA   <- 10000    # spettrale_10k: campi casuali al centro vero
DATASET_VERIFICA <- c("stab_inj1", "stab_inj2")   # dove fare spettrale_10k (piu' lento)
NL_SEED      <- 20261029
# ------------------------------------------------------------------------------------------

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
# BLAS a 1 thread: i dataset girano in processi figli (fork) e il BLAS multithread
# dentro un fork puo' bloccarsi. Il parallelismo e' fra dataset e fra GSEA.
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) RhpcBLASctl::blas_set_num_threads(1)

# ==========================================================================================
# Analisi di un dataset. Gira dentro l'ambiente del checkpoint (vede i suoi oggetti).
# Restituisce le tabelle; l'output a schermo va nel log del dataset in OUT_DIR.
# ==========================================================================================
analizza <- function(threads) {
  BiocParallel::register(BiocParallel::SerialParam())
  stopifnot(exists("Z1_scan"), exists("glm_lin"), exists("all_results"), exists("spe"))
  ALPHA <- if (exists("ALPHA_FDR")) ALPHA_FDR else 0.05
  t_tot <- Sys.time()
  cat("========== PROVA NULLI:", SAMPLE_ID, "==========\n")

  xy <- cbind(spots$pixel_col, spots$pixel_row); n <- nrow(xy)
  D  <- as.matrix(dist(xy))
  stopifnot(identical(colnames(Z1_scan), colnames(spe)), ncol(Z1_scan) == n)
  geni   <- rownames(Z1_scan)
  log_p  <- glm_lin$log_p[match(geni, glm_lin$gene_id)]
  z_true <- Z1_scan[, center_idx]

  # residui e pesi del fit nullo (come la 9c)
  dn  <- if (exists("dge_null")) dge_null else estimateDisp(dge, matrix(1, n, 1))
  phi <- if (!is.null(dn$tagwise.dispersion)) dn$tagwise.dispersion else rep(dn$common.dispersion, nrow(dn))
  mu  <- glmFit(dn, matrix(1, n, 1))$fitted.values
  Y   <- as.matrix(counts(spe))
  den <- 1 + phi * mu
  R   <- (Y - mu) / den; W <- mu / den; S0 <- rowSums(W)
  rm(Y, mu, den, dn); invisible(gc())
  stopifnot(identical(rownames(R), geni))

  # p empirico su |Z| e sigma (sd di Z sotto il nullo, attorno a 0) per ogni gene.
  # X = covariate casuali [spot x campi]; si lavora a blocchi per la memoria.
  p_campi <- function(X, z_ref, blocco = 200) {
    k <- numeric(length(z_ref)); nv <- numeric(length(z_ref)); q <- numeric(length(z_ref))
    for (s in seq(1, ncol(X), by = blocco)) {
      Xb <- X[, s:min(s + blocco - 1, ncol(X)), drop = FALSE]
      S1 <- W %*% Xb; V <- W %*% (Xb^2) - S1^2 / S0; V[V <= 0] <- NA
      Zb <- (R %*% Xb) / sqrt(V)
      k  <- k + rowSums(abs(Zb) >= abs(z_ref), na.rm = TRUE); nv <- nv + rowSums(!is.na(Zb))
      q  <- q + rowSums(Zb^2, na.rm = TRUE)
    }
    sigma <- sqrt(q / pmax(nv, 1))
    list(p = ifelse(is.na(z_ref) | nv == 0, NA_real_, (k + 1) / (nv + 1)),
         p_sigma = ifelse(is.na(z_ref) | nv == 0 | sigma <= 0, NA_real_, 2 * pnorm(-abs(z_ref) / sigma)),
         sigma = sigma)
  }
  z_chk <- { x <- D[, center_idx] / max(D[, center_idx]); S1 <- W %*% x
             as.numeric((R %*% x) / sqrt(W %*% x^2 - S1^2 / S0)) }
  cat("Controllo score test (max |differenza| con Z1_scan):", signif(max(abs(z_chk - z_true), na.rm = TRUE), 3), "\n")

  # nullo 1: centri poco correlati
  centri_corr <- function(j) { r <- cor(D[, j], D)[1, ]; which(abs(r) < SOGLIA_COR & seq_len(n) != j) }
  p_corr <- function(j) {
    cc <- centri_corr(j)
    if (length(cc) < N_MIN_CENTRI) return(list(p = rep(NA_real_, length(geni)), n = length(cc)))
    Zc <- abs(Z1_scan[, cc, drop = FALSE]); zj <- abs(Z1_scan[, j]); nv <- rowSums(!is.na(Zc))
    list(p = ifelse(is.na(zj) | nv == 0, NA_real_, (rowSums(Zc >= zj, na.rm = TRUE) + 1) / (nv + 1)), n = length(cc))
  }

  # nullo 2: randomizzazione spettrale
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

  # centro vero
  x_true <- D[, center_idx]
  pc <- p_corr(center_idx)
  cat("Nullo corr: centri con |r| <", SOGLIA_COR, ":", pc$n, if (pc$n < N_MIN_CENTRI) "(troppo pochi: saltato)" else "", "\n")
  X_true <- campi_spettrali(x_true, B_SPET, NL_SEED)
  mor_raggio <- moran(x_true)
  mor_campi  <- apply(X_true[, seq_len(min(100, B_SPET)), drop = FALSE], 2, moran)
  cat("Moran's I del raggio:", round(mor_raggio, 3), "| dei campi casuali: media", round(mean(mor_campi), 3),
      "sd", round(sd(mor_campi), 3), "\n")
  p_sp <- p_campi(X_true, z_true)
  rm(X_true); invisible(gc())
  cat("sigma di Z sotto il nullo spettrale: mediana", round(median(p_sp$sigma, na.rm = TRUE), 2),
      "| 5%-95%:", paste(round(quantile(p_sp$sigma, c(0.05, 0.95), na.rm = TRUE), 2), collapse = " - "),
      "(1 = nessuna inflazione)\n")
  p_10k <- rep(NA_real_, length(geni))
  if (SAMPLE_ID %in% DATASET_VERIFICA) {
    cat("Verifica risoluzione:", B_VERIFICA, "campi ...\n")
    p_10k <- p_campi(campi_spettrali(x_true, B_VERIFICA, NL_SEED + 999), z_true)$p
    invisible(gc())
  }

  nul <- data.frame(gene_id = geni, z = z_true, p_score = 2 * pnorm(-abs(z_true)),
                    p_corr = pc$p, p_spet = p_sp$p, p_sigma = p_sp$p_sigma, sigma = p_sp$sigma, p_10k = p_10k,
                    log_p = log_p, stringsAsFactors = FALSE) %>%
    left_join(glm_lin %>% select(gene_id, p_QL = p_value, fdr_QL = fdr), by = "gene_id") %>%
    mutate(fdr_corr = p.adjust(p_corr, "BH"), fdr_spet = p.adjust(p_spet, "BH"),
           fdr_sigma = p.adjust(p_sigma, "BH"), fdr_10k = p.adjust(p_10k, "BH"))
  rk <- function(f, p) { r <- order(order(f, p, nul$log_p)); r[is.na(f)] <- NA_integer_; r }
  nul$rank_corr <- rk(nul$fdr_corr, nul$p_corr)
  nul$rank_spet <- rk(nul$fdr_spet, nul$p_spet)
  nul$rank_sigma <- rk(nul$fdr_sigma, nul$p_sigma)
  nul$rank_10k   <- rk(nul$fdr_10k, nul$p_10k)
  cat("Significativi GLM_linear: QL", sum(nul$fdr_QL < ALPHA, na.rm = TRUE),
      "| corr", sum(nul$fdr_corr < ALPHA, na.rm = TRUE), "| spettrale", sum(nul$fdr_spet < ALPHA, na.rm = TRUE),
      "| spettrale_sigma", sum(nul$fdr_sigma < ALPHA, na.rm = TRUE),
      if (any(!is.na(nul$p_10k))) paste("| spettrale_10k", sum(nul$fdr_10k < ALPHA, na.rm = TRUE)) else "", "\n")

  # falsi centri
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
              riga_f(j, "GLM_linear_spettrale", psj$p),
              riga_f(j, "GLM_linear_spettrale_sigma", psj$p_sigma))
  }))
  rm(R, W, V_ev, A, D); invisible(gc())

  # Golden Standard
  METODI <- c("GLM_linear", "GLM_linear_corr", "GLM_linear_spettrale", "GLM_linear_spettrale_sigma",
              if (any(!is.na(nul$p_10k))) "GLM_linear_spettrale_10k",
              intersect(c("SPATA2", "SpatialDE_orig"), unique(as.character(all_results$Method))))
  if (pc$n < N_MIN_CENTRI) METODI <- setdiff(METODI, "GLM_linear_corr")
  gs_ids <- unique(all_results$gene_id[all_results$is_GS %in% TRUE])
  prova <- bind_rows(
    all_results %>% filter(Method %in% c("GLM_linear", "SPATA2", "SpatialDE_orig")) %>%
      mutate(Method = as.character(Method)) %>% select(Method, gene_id, fdr, rank, tested),
    nul %>% transmute(Method = "GLM_linear_corr", gene_id, fdr = fdr_corr, rank = rank_corr, tested = !is.na(fdr_corr)),
    nul %>% transmute(Method = "GLM_linear_spettrale", gene_id, fdr = fdr_spet, rank = rank_spet, tested = !is.na(fdr_spet)),
    nul %>% transmute(Method = "GLM_linear_spettrale_sigma", gene_id, fdr = fdr_sigma, rank = rank_sigma, tested = !is.na(fdr_sigma)),
    nul %>% transmute(Method = "GLM_linear_spettrale_10k", gene_id, fdr = fdr_10k, rank = rank_10k, tested = !is.na(fdr_10k)),
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

  # pathway: stessa costruzione della sezione 15
   
   
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

  # GSEA come la sezione 16
  pos_score <- function(ids) setNames(rev(seq_along(ids)), ids)
  SC <- lapply(setNames(METODI, METODI), function(m) {
    d <- prova %>% filter(Method == m, tested, !is.na(rank)) %>% arrange(rank); pos_score(d$gene_id)
  })
  SC$Expression_control <- pos_score(U_enr[order(-expr_U)])
  gsea_nul <- bind_rows(parallel::mclapply(names(SC), mc.cores = max(1, min(threads, length(SC))), function(m) {
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
  cat("Fatto in", round(as.numeric(difftime(Sys.time(), t_tot, units = "mins")), 1), "min\n")

  ds <- SAMPLE_ID
  list(ranks = prova %>% select(Method, gene_id, fdr, rank, tested, is_GS, sig) %>% mutate(dataset = ds),
       gsea  = gsea_nul %>% mutate(dataset = ds),
       geni  = tab_geni %>% mutate(dataset = ds),
       tgsea = tab_gsea %>% mutate(dataset = ds),
       falsi = tab_falsi,
       pv    = nul %>% select(gene_id, p_QL, p_score, p_corr, p_spet, p_sigma, p_10k, sigma) %>% mutate(dataset = ds),
       riass = data.frame(dataset = ds, n_spot = n, centri_corr = pc$n,
                          moran_raggio = round(mor_raggio, 3), moran_campi = round(mean(mor_campi), 3),
                          n_sig_QL = sum(nul$fdr_QL < ALPHA, na.rm = TRUE),
                          n_sig_corr = sum(nul$fdr_corr < ALPHA, na.rm = TRUE),
                          n_sig_spet = sum(nul$fdr_spet < ALPHA, na.rm = TRUE),
                          n_sig_sigma = sum(nul$fdr_sigma < ALPHA, na.rm = TRUE),
                          n_sig_10k = if (any(!is.na(nul$p_10k))) sum(nul$fdr_10k < ALPHA, na.rm = TRUE) else NA_integer_,
                          sigma_mediana = round(median(nul$sigma, na.rm = TRUE), 2)))
}

# ==========================================================================================
# Un processo per dataset: carica il checkpoint in un ambiente proprio e lancia analizza()
# ==========================================================================================
presenti <- CHECKPOINTS[file.exists(CHECKPOINTS)]
for (ck in setdiff(CHECKPOINTS, presenti)) cat("MANCA:", ck, "\n")
stopifnot(length(presenti) > 0)
n_jobs  <- max(1L, min(MAX_JOBS, length(presenti)))
threads <- max(1L, floor(N_CPU / n_jobs))
cat("Dataset:", length(presenti), "| in parallelo:", n_jobs, "| processori per dataset:", threads, "\n")

uno <- function(ck) {
  nome <- basename(dirname(ck))
  log_f <- file.path(OUT_DIR, paste0("log_", nome, ".txt"))
  zz <- file(log_f, open = "wt"); sink(zz); sink(zz, type = "message")
  on.exit({ sink(type = "message"); sink(); close(zz) }, add = TRUE)
  t0 <- Sys.time()
  res <- tryCatch({
    e <- new.env(parent = globalenv())
    load(ck, envir = e)
    f <- analizza; environment(f) <- e
    f(threads)
  }, error = function(err) { cat("ERRORE:", conditionMessage(err), "\n"); NULL })
  list(nome = nome, res = res, min = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1))
}
t_tot <- Sys.time()
esiti <- if (n_jobs > 1) parallel::mclapply(presenti, uno, mc.cores = n_jobs, mc.preschedule = FALSE) else lapply(presenti, uno)
for (e in esiti) {
  if (inherits(e, "try-error")) { cat("\nProcesso fallito:", as.character(e), "\n"); next }
  cat("\n#################", e$nome, "(", e$min, "min ) #################\n")
  cat(readLines(file.path(OUT_DIR, paste0("log_", e$nome, ".txt"))), sep = "\n")
}
cat("\nTempo totale:", round(as.numeric(difftime(Sys.time(), t_tot, units = "mins")), 1), "min\n")

ok <- Filter(function(e) !inherits(e, "try-error") && !is.null(e$res), esiti)
if (length(ok) == 0) stop("Nessun dataset ha prodotto risultati: vedi i log in ", OUT_DIR)
unisci  <- function(k) bind_rows(lapply(ok, function(e) e$res[[k]]))
R_ranks <- unisci("ranks"); R_gsea <- unisci("gsea")
geni <- unisci("geni"); gsea <- unisci("tgsea"); falsi <- unisci("falsi"); pv <- unisci("pv"); riass <- unisci("riass")

# ==========================================================================================
# Tabelle
# ==========================================================================================
larga  <- function(d, col) as.data.frame(d %>% select(Method, dataset, all_of(col)) %>%
                                           pivot_wider(names_from = dataset, values_from = all_of(col)))
mostra <- function(titolo, d) { cat("\n>>>", titolo, "<<<\n"); print(d, row.names = FALSE) }
mostra("Riassunto per dataset", as.data.frame(riass))
fs <- falsi %>% group_by(dataset, Method) %>%
  summarise(frac_p05 = round(mean(frac_p05), 3), n_bh = round(mean(n_bh), 1), .groups = "drop")
mostra("Falsi centri: frazione di geni con p < 0.05 (affidabile = 0.05)", larga(fs, "frac_p05"))
mostra("Falsi centri: geni BH-significativi, media (affidabile = circa 0)", larga(fs, "n_bh"))
mostra("% geni significativi", larga(geni, "pct_sig"))
mostra("Golden Standard recuperato (%)", larga(geni, "GS_recall"))
mostra("Golden Standard: fold sul caso", larga(geni, "GS_fold"))
mostra("AUC sul Golden Standard", larga(geni, "AUC_GS"))
mostra("GSEA: AUC contesto vs altri sul NES", larga(gsea, "auc_ctx_nes"))
mostra("GSEA: correlazione del NES con il controllo sull'espressione", larga(gsea, "cor_con_espressione"))
mostra("GSEA: pathway con FDR < 0.05 (solo diagnostica)", larga(gsea, "n_p05"))
mostra("GSEA: pathway di contesto con FDR < 0.05 (solo diagnostica)", larga(gsea, "n_p05_ctx"))
for (x in list(list(riass, "riassunto"), list(geni, "tab_geni"), list(gsea, "tab_gsea"),
               list(falsi, "falsi_centri"), list(pv, "pvalori"), list(R_gsea, "gsea_pathway")))
  write_csv(x[[1]], file.path(OUT_DIR, paste0("nulli_", x[[2]], ".csv")))

# --- 3. Grafici della calibrazione ------------------------------------------------------------
COL_T <- c(GLM_linear = "#7570b3", GLM_linear_corr = "#1b9e77", GLM_linear_spettrale = "#e6ab02",
           GLM_linear_spettrale_sigma = "#a6761d", GLM_linear_spettrale_10k = "#66a61e")
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
  print(pf); ggsave(file.path(OUT_DIR, "nulli_00a_falsi_centri_frac.png"), pf, width = 11, height = 5.5, dpi = 150)
  pb <- falsi %>% mutate(Method = factor(Method, levels = rev(names(COL_T)))) %>%
    ggplot(aes(x = n_bh + 1, y = Method, colour = Method)) +
    geom_point(size = 1.8, alpha = 0.8, position = position_jitter(height = 0.15, width = 0, seed = 1)) + scale_x_log10() +
    facet_wrap(~ dataset, ncol = 3) + scale_colour_manual(values = COL_T, guide = "none") + tema0 +
    labs(title = "Falsi centri: geni BH-significativi", subtitle = "scala log (n + 1) | 1 = nessun gene, come atteso con un test calibrato",
         x = "Geni con FDR < 0.05 (+1)", y = NULL)
  print(pb); ggsave(file.path(OUT_DIR, "nulli_00b_falsi_centri_bh.png"), pb, width = 11, height = 5.5, dpi = 150)
}
if (nrow(pv) > 0) {
  ph <- pv %>% pivot_longer(c(p_QL, p_corr, p_spet, p_sigma), names_to = "Test", values_to = "p") %>% filter(!is.na(p)) %>%
    mutate(Test = factor(recode(Test, p_QL = "GLM parametrico (QL)", p_corr = "nullo corr", p_spet = "spettrale",
                                p_sigma = "spettrale sigma"),
                         levels = c("GLM parametrico (QL)", "nullo corr", "spettrale", "spettrale sigma"))) %>%
    ggplot(aes(x = p)) + geom_histogram(breaks = seq(0, 1, 0.05), fill = "grey55", colour = "white") +
    facet_grid(dataset ~ Test, scales = "free_y") + tema0 + theme(strip.text.y = element_text(angle = 0)) +
    labs(title = "P-value al centro vero", subtitle = "barra a sinistra = geni con segnale | piatto altrove = test calibrato sui nulli",
         x = "p", y = "geni")
  print(ph); ggsave(file.path(OUT_DIR, "nulli_00c_istogrammi_p.png"), ph, width = 10, height = 1.6 * n_distinct(pv$dataset) + 1, dpi = 150)
  psg <- ggplot(pv %>% filter(is.finite(sigma)), aes(x = sigma)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_histogram(bins = 60, fill = "#a6761d", colour = "white") + scale_x_log10() +
    facet_wrap(~ dataset, ncol = 3, scales = "free_y") + tema0 +
    labs(title = "Inflazione di Z per gene sotto il nullo spettrale (sigma)",
         subtitle = "1 = come prevede il test parametrico | 3 = Z tre volte piu' grandi del previsto anche senza un gradiente vero",
         x = "sigma (scala log)", y = "geni")
  print(psg); ggsave(file.path(OUT_DIR, "nulli_00d_sigma.png"), psg, width = 11, height = 5.5, dpi = 150)
}

# --- 4. Golden Standard e GSEA (come la pipeline) -----------------------------------------------
ALPHA_G <- 0.05; B_GS <- 1000; B_CTX <- 1000; SEED_G <- 20261029
LIV <- c("GLM_linear", "GLM_linear_corr", "GLM_linear_spettrale", "GLM_linear_spettrale_sigma",
         "GLM_linear_spettrale_10k", "SPATA2", "SpatialDE_orig", "Expression_control")
COL <- c(COL_T, SPATA2 = "#d95f02", SpatialDE_orig = "#e7298a", Expression_control = "grey50")
rk_all <- R_ranks
gs_all <- R_gsea
if (nrow(rk_all) == 0) stop("Nessun dataset ha prodotto risultati: vedi i log")
DS  <- unique(rk_all$dataset)
LIV <- intersect(LIV, unique(c(rk_all$Method, gs_all$Method)))
rk_all <- rk_all %>% filter(tested, !is.na(rank)) %>%
  mutate(Method = factor(Method, levels = LIV), dataset = factor(dataset, levels = DS))
gs_all <- gs_all %>% mutate(Method = factor(Method, levels = LIV), dataset = factor(dataset, levels = DS))

salva_g <- function(p, nome, w, h) {
  print(p); ggsave(file.path(OUT_DIR, paste0(nome, ".png")), p, width = w, height = h, dpi = 150)
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
  write_csv(auc_ci, file.path(OUT_DIR, "nulli_auc_gs_bootstrap.csv"))
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
  write_csv(ctx, file.path(OUT_DIR, "nulli_auc_contesto_bootstrap.csv"))
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
cat("Grafici di tutti i dataset in", OUT_DIR, "\n")

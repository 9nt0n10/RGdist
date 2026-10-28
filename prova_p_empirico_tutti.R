# ##############################################################################
# ##   PROVA P EMPIRICO SU TUTTI I DATASET                                    ##
# ##                                                                          ##
# ##   Per ogni checkpoint lancia prova_p_empirico.R (e i suoi grafici) in    ##
# ##   un processo R separato: ogni dataset parte pulito e la memoria del     ##
# ##   checkpoint precedente viene liberata. Poi legge le tabelle scritte da  ##
# ##   ogni run (nella cartella del checkpoint) e le mette a confronto.       ##
# ##############################################################################

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })

SCRIPT_DIR  <- "/home/ammendola/laboratorio/spatial/RGdist"     # cartella con gli script prova_p_empirico*.R
CONFRONTO_DIR <- "/home/ammendola/laboratorio/spatial/prova_p_empirico_confronto"
RIFAI       <- TRUE      # FALSE = non rilancia, legge solo le tabelle gia' scritte
N_CPU       <- 40        # processori in tutto
MAX_JOBS    <- 6         # dataset in parallelo: abbassalo se la RAM non basta (ogni run carica un checkpoint)

CHECKPOINTS <- c(
  "/home/ammendola/laboratorio/spatial/topo/code/benchmark_v7/heme_1000/checkpoint_pre_grafici.RData",
  "/home/ammendola/laboratorio/spatial/brain_injury/code/stab_inj1/benchmark_v7/stab_inj1/checkpoint_pre_grafici.RData",
  "/home/ammendola/laboratorio/spatial/brain_injury/code/stab_inj2/benchmark_v7/stab_inj2/checkpoint_pre_grafici.RData",
  "/home/ammendola/laboratorio/spatial/glioblastoma_mice/code/benchmark_v7/gbm/checkpoint_pre_grafici.RData",
  "/home/ammendola/laboratorio/spatial/periodontitis/code/benchmark_v7/perio/checkpoint_pre_grafici.RData",
  "/home/ammendola/laboratorio/spatial/tumor/full_tissue/code/benchmark_v7/head_and_neck/checkpoint_pre_grafici.RData"
)

core    <- file.path(SCRIPT_DIR, "prova_p_empirico.R")
grafici <- file.path(SCRIPT_DIR, "prova_p_empirico_grafici.R")
stopifnot(file.exists(core), file.exists(grafici))
dir.create(CONFRONTO_DIR, recursive = TRUE, showWarnings = FALSE)
rscript <- file.path(R.home("bin"), "Rscript")       # stesso R di questa sessione

# --- 1. Un processo per dataset, fino a MAX_JOBS insieme -------------------------------
presenti <- CHECKPOINTS[file.exists(CHECKPOINTS)]
for (ck in setdiff(CHECKPOINTS, presenti)) cat("MANCA:", ck, "\n")
n_jobs  <- max(1L, min(MAX_JOBS, length(presenti)))
threads <- max(1L, floor(N_CPU / n_jobs))
cat("Dataset:", length(presenti), "| in parallelo:", n_jobs, "| processori per dataset:", threads, "\n")

lancia <- function(ck) {
  t0  <- Sys.time()
  env <- c(paste0("PROVA_THREADS=", threads), "OMP_NUM_THREADS=1", "OPENBLAS_NUM_THREADS=1", "MKL_NUM_THREADS=1")
  out <- suppressWarnings(system2(rscript, c(shQuote(core), shQuote(ck), shQuote(grafici)),
                                  stdout = TRUE, stderr = TRUE, env = env))
  writeLines(out, file.path(dirname(ck), "prova_p_empirico_log.txt"))
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

# --- 2. Tabelle di tutti i dataset ----------------------------------------------------
leggi <- function(f) bind_rows(lapply(CHECKPOINTS, function(ck) {
  p <- file.path(dirname(ck), f)
  if (file.exists(p)) suppressMessages(read_csv(p, show_col_types = FALSE)) else NULL
}))
riass <- leggi("prova_p_empirico_tab_riassunto.csv")
geni  <- leggi("prova_p_empirico_tab_geni.csv")
gsea  <- leggi("prova_p_empirico_tab_gsea.csv")
spt   <- leggi("prova_p_empirico_tab_spata2.csv")
if (nrow(geni) == 0) stop("Nessuna tabella trovata: controlla i log prova_p_empirico_log.txt")

larga <- function(d, col) as.data.frame(d %>% select(Method, dataset, all_of(col)) %>%
                                          pivot_wider(names_from = dataset, values_from = all_of(col)))
mostra <- function(titolo, d) { cat("\n>>>", titolo, "<<<\n"); print(d, row.names = FALSE) }

mostra("Riassunto per dataset", as.data.frame(riass))
mostra("% geni significativi", larga(geni, "pct_sig"))
mostra("Golden Standard recuperato (%)", larga(geni, "GS_recall"))
mostra("Golden Standard: fold sul caso", larga(geni, "GS_fold"))
mostra("AUC sul Golden Standard", larga(geni, "AUC_GS"))
mostra("GSEA: AUC contesto vs altri sul NES", larga(gsea, "auc_ctx_nes"))
mostra("GSEA: correlazione del NES con il controllo sull'espressione", larga(gsea, "cor_con_espressione"))
mostra("GSEA: pathway con FDR < 0.05 (solo diagnostica)", larga(gsea, "n_p05"))
mostra("GSEA: pathway di contesto con FDR < 0.05 (solo diagnostica)", larga(gsea, "n_p05_ctx"))
if (nrow(spt) > 0) {
  mostra("SPATA2: % dei suoi significativi ritrovati", larga(spt, "pct_SPATA2_anche_qui"))
  mostra("SPATA2: % attesa a caso", larga(spt, "pct_atteso_a_caso"))
  mostra("SPATA2: Spearman dei rank", larga(spt, "spearman_con_SPATA2"))
}
write_csv(geni, file.path(CONFRONTO_DIR, "tab_geni_tutti.csv"))
write_csv(gsea, file.path(CONFRONTO_DIR, "tab_gsea_tutti.csv"))
write_csv(riass, file.path(CONFRONTO_DIR, "tab_riassunto_tutti.csv"))
if (nrow(spt) > 0) write_csv(spt, file.path(CONFRONTO_DIR, "tab_spata2_tutti.csv"))

# --- 3. Grafico riassuntivo -------------------------------------------------------------
LIV <- c("GLM_linear", "GLM_linear_emp", "GLM_linear_emp_abs", "SPATA2", "SpatialDE_orig", "Expression_control")
lunga <- bind_rows(
  geni %>% select(Method, dataset, `% geni significativi` = pct_sig, `AUC Golden Standard` = AUC_GS,
                  `Golden Standard: fold` = GS_fold) %>%
    pivot_longer(-c(Method, dataset), names_to = "Misura", values_to = "Valore"),
  gsea %>% select(Method, dataset, `GSEA: AUC contesto` = auc_ctx_nes,
                  `GSEA: correlazione con espressione` = cor_con_espressione) %>%
    pivot_longer(-c(Method, dataset), names_to = "Misura", values_to = "Valore")) %>%
  filter(is.finite(Valore)) %>%
  mutate(Method = factor(Method, levels = rev(intersect(LIV, unique(Method)))))
p <- ggplot(lunga, aes(x = Valore, y = Method, colour = dataset)) +
  geom_point(size = 2, alpha = 0.85, position = position_jitter(height = 0.15, width = 0, seed = 1)) +
  facet_wrap(~ Misura, scales = "free_x", ncol = 3) +
  theme_minimal(base_size = 9) + theme(legend.position = "bottom") +
  labs(title = "Prova p empirico: tutti i dataset", x = NULL, y = NULL, colour = NULL)
print(p)
ggsave(file.path(CONFRONTO_DIR, "confronto_dataset.png"), p, width = 12, height = 7, dpi = 150)
cat("\nTabelle e grafico in", CONFRONTO_DIR, "\n")

# --- 4. Grafici come nella pipeline, tutti i dataset insieme --------------------------------
f_graf <- file.path(SCRIPT_DIR, "prova_p_empirico_grafici_tutti.R")
if (file.exists(f_graf)) source(f_graf)

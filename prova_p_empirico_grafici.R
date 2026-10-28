# ##############################################################################
# ##   GRAFICI DELLA PROVA SUL P EMPIRICO                                     ##
# ##                                                                          ##
# ##   Da lanciare DOPO prova_p_empirico.R, nella stessa sessione: usa gli    ##
# ##   oggetti rimasti in memoria (prova, gsea_prova, PLOT_DIR, ALPHA_FDR).   ##
# ##   Nessun calcolo pesante, nessuna GSEA rifatta.                         ##
# ##     1. boxplot dei rank dei geni del Golden Standard                     ##
# ##     2. curve ROC sul Golden Standard                                    ##
# ##     3. AUC contesto/altri sul NES, con IC bootstrap                     ##
# ##     4. NES dei pathway di contesto e degli altri (boxplot)              ##
# ##     5. NES dei pathway di contesto piu' alti, metodo per metodo         ##
# ##     6. pathway con FDR < ALPHA (solo diagnostica)                       ##
# ##############################################################################

if (!exists("prova") || !exists("gsea_prova"))
  stop("Lancia prima prova_p_empirico.R nella stessa sessione.")
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })

LIV  <- c("GLM_linear", "GLM_linear_emp", "GLM_linear_emp_abs", "SPATA2", "SpatialDE_orig", "Expression_control")
COL  <- c(GLM_linear = "#7570b3", GLM_linear_emp = "#1b9e77", GLM_linear_emp_abs = "#66a61e", SPATA2 = "#d95f02",
          SpatialDE_orig = "#e7298a", Expression_control = "grey50")
B    <- 1000
SEED <- 20261028
salva <- function(p, nome, w = 8, h = 5) {
  print(p); ggsave(file.path(PLOT_DIR, paste0("prova_emp_", nome, ".png")), p, width = w, height = h, dpi = 150)
}
gr <- prova %>% filter(tested, !is.na(rank)) %>% mutate(Method = factor(Method, levels = LIV))

# --- 1. Boxplot dei rank del Golden Standard -------------------------------------------
gs <- gr %>% filter(is_GS)
med <- gs %>% group_by(Method) %>% summarise(m = median(rank), .groups = "drop")
p1 <- ggplot(gs, aes(x = rank, y = Method, fill = Method)) +
  geom_boxplot(outlier.size = 0.6, width = 0.6) +
  geom_text(data = med, aes(x = m, y = Method, label = round(m)), vjust = -1.4, size = 3, inherit.aes = FALSE) +
  scale_fill_manual(values = COL, guide = "none") + scale_y_discrete(limits = rev) +
  theme_minimal(base_size = 10) +
  labs(title = "Rank dei geni del Golden Standard", subtitle = paste0(n_distinct(gs$gene_id), " geni | piu' a sinistra = meglio"),
       x = "Rank", y = NULL)
salva(p1, "01_boxplot_gs")

# --- 2. Curve ROC ------------------------------------------------------------------------
roc <- gr %>% group_by(Method) %>% arrange(rank, .by_group = TRUE) %>%
  mutate(TPR = cumsum(is_GS) / sum(is_GS), FPR = cumsum(!is_GS) / sum(!is_GS)) %>%
  select(Method, FPR, TPR) %>% ungroup()
roc <- bind_rows(roc, data.frame(Method = factor(unique(roc$Method), levels = LIV), FPR = 0, TPR = 0)) %>%
  arrange(Method, FPR, TPR)
auc_tab <- roc %>% group_by(Method) %>%
  summarise(AUC = sum(diff(FPR) * (head(TPR, -1) + tail(TPR, -1)) / 2), .groups = "drop")
lab <- setNames(paste0(auc_tab$Method, " (AUC ", sprintf("%.3f", auc_tab$AUC), ")"), auc_tab$Method)
p2 <- ggplot(roc, aes(x = FPR, y = TPR, colour = Method)) +
  geom_abline(linetype = "dashed", colour = "grey60") + geom_line(linewidth = 0.7) +
  scale_colour_manual(values = COL, labels = lab) + coord_equal() +
  theme_minimal(base_size = 10) +
  labs(title = "ROC sul Golden Standard", x = "False positive rate", y = "True positive rate", colour = NULL)
salva(p2, "02_roc", w = 7, h = 6)

# --- 3. AUC contesto/altri sul NES con IC bootstrap -------------------------------------
auc_r <- function(r, i, n1, n0) (sum(r[i]) - n1 * (n1 + 1) / 2) / (n1 * n0)
ctx <- bind_rows(lapply(split(gsea_prova, gsea_prova$Method), function(d) {
  d <- d[!is.na(d$NES), ]; y <- d$is_context; n1 <- sum(y); n0 <- sum(!y)
  a <- auc_r(rank(d$NES), which(y), n1, n0)
  set.seed(SEED)
  bt <- replicate(B, auc_r(rank(d$NES[c(sample(which(y), replace = TRUE), sample(which(!y), replace = TRUE))]),
                           seq_len(n1), n1, n0))
  data.frame(Method = d$Method[1], auc = a, lo = quantile(bt, 0.025), hi = quantile(bt, 0.975))
})) %>% mutate(Method = factor(Method, levels = LIV))
p3 <- ggplot(ctx, aes(x = auc, y = Method, colour = Method)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.25, orientation = "y") + geom_point(size = 2.5) +
  scale_colour_manual(values = COL, guide = "none") + scale_y_discrete(limits = rev) +
  theme_minimal(base_size = 10) +
  labs(title = "Pathway di contesto contro gli altri, sul NES",
       subtitle = paste0(sum(gsea_prova$is_context[gsea_prova$Method == LIV[1]]), " pathway di contesto | 0.5 = nessuna preferenza | IC 95% bootstrap"),
       x = "AUC (contesto vs altri)", y = NULL)
salva(p3, "03_auc_contesto")

# --- 4. NES: contesto contro altri --------------------------------------------------------
p4 <- gsea_prova %>% mutate(Method = factor(Method, levels = LIV),
                            Pathway = ifelse(is_context, "Contesto", "Altri")) %>%
  ggplot(aes(x = NES, y = Method, fill = Pathway)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_boxplot(outlier.size = 0.4, linewidth = 0.3) +
  scale_fill_manual(values = c(Contesto = "#d95f02", Altri = "grey75")) + scale_y_discrete(limits = rev) +
  theme_minimal(base_size = 10) +
  labs(title = "NES di tutti i pathway", subtitle = "nessun filtro di significativita'", x = "NES", y = NULL, fill = NULL)
salva(p4, "04_nes_boxplot")

# --- 5. NES dei pathway di contesto piu' alti ----------------------------------------------
top_pw <- gsea_prova %>% filter(is_context, Method != "Expression_control") %>%
  group_by(pathway) %>% summarise(m = mean(NES, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(m)) %>% slice_head(n = 30) %>% pull(pathway)
p5 <- gsea_prova %>% filter(pathway %in% top_pw) %>%
  mutate(Method = factor(Method, levels = LIV),
         Label = factor(gsub("_", " ", gsub("^(HALLMARK|GOBP)_", "", pathway)),
                        levels = rev(gsub("_", " ", gsub("^(HALLMARK|GOBP)_", "", top_pw)))),
         sig = p_adjust < ALPHA_FDR) %>%
  ggplot(aes(x = Method, y = Label, colour = NES, size = -log10(p_adjust))) +
  geom_point() + geom_point(data = function(d) d[d$sig, ], shape = 1, colour = "black", stroke = 0.6) +
  scale_colour_viridis_c() + theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(title = "NES dei 30 pathway di contesto piu' alti (media sui metodi)",
       subtitle = paste0("cerchio nero = FDR < ", ALPHA_FDR), x = NULL, y = NULL, size = "-log10(FDR)")
salva(p5, "05_nes_contesto_top", w = 9, h = 8)

# --- 6. Pathway con FDR < ALPHA (solo diagnostica) -----------------------------------------
p6 <- gsea_prova %>% group_by(Method) %>%
  summarise(Tutti = sum(p_adjust < ALPHA_FDR, na.rm = TRUE),
            Contesto = sum(p_adjust < ALPHA_FDR & is_context, na.rm = TRUE), .groups = "drop") %>%
  pivot_longer(c(Tutti, Contesto), names_to = "Pathway", values_to = "n") %>%
  mutate(Method = factor(Method, levels = LIV), Pathway = factor(Pathway, levels = c("Tutti", "Contesto"))) %>%
  ggplot(aes(x = n, y = Method, fill = Method)) + geom_col() +
  geom_text(aes(label = n), hjust = -0.2, size = 3) +
  facet_wrap(~ Pathway, scales = "free_x") +
  scale_fill_manual(values = COL, guide = "none") + scale_y_discrete(limits = rev) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_minimal(base_size = 10) +
  labs(title = paste0("Pathway con FDR < ", ALPHA_FDR, " (solo diagnostica)"), x = "Numero di pathway", y = NULL)
salva(p6, "06_pathway_fdr", w = 9, h = 4.5)

cat("Grafici salvati in", PLOT_DIR, "con prefisso prova_emp_\n")

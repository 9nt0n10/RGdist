# Questo codice calcola la correlazione tra il fold change dei conteggi genici e la distanza dai punti di partenza, espressa in anelli concentrici. 
# La mediana delle correlazioni viene poi visualizzata in un grafico per ogni punto di partenza considerato uno alla volta.
# In dettaglio, il codice esegue i seguenti passaggi:
# 1. Per ogni punto di partenza, calcola la distanza euclidea di tutti gli altri punti di un dataset.
# 2. Classifica queste distanze in intervalli (anelli concentrici) usando i quantili.
# 3. Calcola la correlazione tra il fold change dei conteggi genici e l'assegnazione agli anelli per ogni gene.
# 4. Per ogni punto di partenza, calcola la mediana delle correlazioni assolute dei geni significativi (p-value < 0.05).
# 5. Visualizza la mediana delle correlazioni in un grafico, dove i colori rappresentano i valori di correlazione mediana in base alla posizione degli spot.
# Calcola e visualizza la mediana delle correlazioni come descritto sopra

library(SpatialExperiment)
library(rtracklayer)
library(lobstr)
library(readr)
library(scater)
library(sva)
library(plyr)
library(DropletUtils)
library(ggplot2)
library(foreach)
library(doParallel)

# Importa i dati come oggetto SpatialExperiment
spe <- SpatialExperiment::read10xVisium(
  samples = "path/to/samples/heme_1000/outs/",
  sample_id = "heme_1000",
  type = c("HDF5", "sparse"),
  data = "filtered",
  images = "lowres", 
  load = TRUE
)

# Variabili continue per numero di UMI e geni per ogni campione
spe$sum_umi <- colSums(counts(spe))
spe$sum_gene <- colSums(counts(spe) > 0)

# Allineamento al genoma di riferimento
gtf <- rtracklayer::import("path/to/gencode.vM23.annotation.gtf")
gtf <- gtf[gtf$type == "gene"]
gtf$gene_id <- gsub("\\..*", "", gtf$gene_id)
names(gtf) <- gtf$gene_id

match_genes <- match(rownames(spe), gtf$gene_id)
spe <- spe[!is.na(match_genes), ] 
match_genes <- match_genes[!is.na(match_genes)]
mcols(gtf) <- mcols(gtf)[, c("source", "type", "gene_id", "gene_name", "gene_type")]
rowRanges(spe) <- gtf[match_genes]

rowData(spe)$gene_search <- paste0(
  rowData(spe)$gene_name, "; ", rowData(spe)$gene_id
)

# Informazioni sull'anatomia e coordinate del sito di iniezione
spots_heme_1000 <- read_csv("path/to/tissue_positions_list.csv", col_names = FALSE)
names(spots_heme_1000) <- c("Barcode", "in_tissue", "array_row", "array_col", "pixel_row", "pixel_col")
spots_heme_1000 <- subset(spots_heme_1000, in_tissue == 1)

anatomy <- merge(read_csv("path/to/heme_1000_anatomy.csv"),
                 read_csv("path/to/heme_1000_injection_site.csv"),
                 by = "Barcode")

spots_heme_1000 <- merge(spots_heme_1000, anatomy, by = "Barcode")

# Aggiungi informazioni all'oggetto spe
spe$anatomy <- spots_heme_1000$anatomy
spe$inj_site <- spots_heme_1000$injection_site

# heme_1000lo qualità e filtraggio dei dati
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
spe <- addPerCellQC(spe, subsets = list(mito = is_mito))

hist(colData(spe)$detected, breaks = 20)

spe <- spe[!is_mito, ]

# Rimuove gli spot non assegnati a nessun cluster anatomico
NA_spot <- which(is.na(spots_heme_1000$anatomy))
spe <- spe[, -NA_spot, drop = FALSE]
spots_heme_1000 <- spots_heme_1000[-NA_spot, ]

# Rimuove i geni non espressi
no_expr <- which(rowSums(counts(spe)) == 0)
spe <- spe[-no_expr, , drop = FALSE]

summary(spe$sum)

if (any(spe$sum == 0)) {
  spots_no_counts <- which(spe$sum == 0)
  print(length(spots_no_counts))
  print(length(spots_no_counts) / ncol(spe) * 100)
  spe <- spe[, -spots_no_counts, drop = FALSE]
}

# Rimuove i geni espressi in meno di 3 spot
h1000_counts <- as.matrix(assays(spe)$counts)
no_rel <- which(rowSums(h1000_counts != 0) <= 2)
spe <- spe[-no_rel, , drop = FALSE]

# Normalizzazione dei conteggi
spe <- scuttle::logNormCounts(spe)

h1000_logcounts <- as.matrix(assays(spe)$logcounts)

gene_counts_median <- apply(h1000_logcounts, 1, function(x) median(x[x != 0]))
gene_counts_median <- matrix(gene_counts_median, nrow = length(gene_counts_median), ncol = 1, byrow = TRUE,
                             dimnames = list(rownames(h1000_logcounts)))

h1000_ratio <- matrix(NA, nrow = nrow(h1000_logcounts), ncol = ncol(h1000_logcounts),
                      dimnames = list(rownames(h1000_logcounts), colnames(h1000_logcounts)))

h1000_logcounts[h1000_logcounts == 0] <- NA

for (i in 1:nrow(h1000_logcounts)) {
  no_na <- !is.na(h1000_logcounts[i, ])
  h1000_ratio[i, no_na] <- h1000_logcounts[i, no_na] / gene_counts_median[i, 1]
}

logratio <- h1000_ratio
logratio[!is.na(logratio)] <- log2(logratio[!is.na(logratio)])

# Step 7: Filter rows with less than 70% NA values
percentuale_na <- apply(logratio, 1, function(x) sum(is.na(x)) / length(x) * 100)
logratio_filtered <- logratio[percentuale_na < 70, ]

# Imposta il numero di core da utilizzare
num_cores <- 100  # Ad esempio, utilizza 4 core

# Inizializza il backend parallelo
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Iterazione su tutti i punti
cor_results <- foreach(j = 1:nrow(spots_heme_1000), .combine = rbind, .packages = c('foreach', 'doParallel')) %dopar% {
  distances <- sqrt((spots_heme_1000[, 5] - spots_heme_1000[j, 5])^2 + (spots_heme_1000[, 6] - spots_heme_1000[j, 6])^2)
  
  cor_before <- numeric(nrow(logratio_filtered))
  pval_before <- numeric(nrow(logratio_filtered))
  
  # Calcola gli intervalli di distanza usando quantili
  quantiles <- quantile(distances, probs = seq(0, 1, length.out = 10))  # 10 intervalli, puoi cambiare il numero
  interval_indices <- .bincode(distances, breaks = quantiles, include.lowest = TRUE)
  
  # Assegna gli intervalli di distanza ai punti di spots
  spots_heme_1000$ring <- interval_indices
  
  # Calcola la correlazione per ogni gene
  for (i in 1:nrow(logratio_filtered)) {
    zero <- which(is.na(logratio_filtered[i, ]))    
    if (length(zero) > 0) {
      cor_before[i] <- cor(spots_heme_1000$ring[-zero], logratio_filtered[i, -zero], method = "spearman")
      test <- cor.test(spots_heme_1000$ring[-zero], logratio_filtered[i, -zero], method = "spearman")
      pval_before[i] <- test$p.value
    } else {
      cor_before[i] <- cor(spots_heme_1000$ring, logratio_filtered[i,], method = "spearman")
      test <- cor.test(spots_heme_1000$ring, logratio_filtered[i,], method = "spearman")
      pval_before[i] <- test$p.value
    }
  }
  
  # Combina i risultati
  data.frame(gene_id = rownames(logratio_filtered), cor_before = cor_before, pval_before = pval_before)
}

stopCluster(cl)

# Filtra il dataframe

cor_nodup <- aggregate(. ~ gene_id, data = cor_results, FUN = mean)
cor_filtered <- subset(cor_nodup, cor_nodup[, 3] <= 0.05)
# Passo 2: Ordina il dataframe in base al valore assoluto della colonna 2 in ordine decrescente
cor_nodup_main <- cor_nodup[order(-abs(cor_nodup[, 2])), ]

# Tieni solo i primi top geni
#cor_nodup_main <- head(cor_nodup_main, 200)
#cor_nodup_main <- head(cor_nodup_main, 100)
#cor_nodup_main <- head(cor_nodup_main, 50)
#cor_nodup_main <- head(cor_nodup_main, 10)
# Passo 1: Imposta la prima colonna come rownames
rownames(cor_nodup_main) <- cor_nodup_main[, 1]

# Rimuovi la prima colonna dal dataframe
cor_nodup_heme_1000 <- cor_nodup_main[, -1]
logratio_filtered_heme_1000 <- logratio_filtered[rownames(cor_nodup_heme_1000), , drop = FALSE]

num_cores <- 100  # Ad esempio, utilizza 4 core

# Inizializza il backend parallelo
cl <- makeCluster(num_cores)
registerDoParallel(cl)

ring_median_heme_1000 <- foreach(j = 1:nrow(spots_heme_1000), .combine = rbind, .packages = c('foreach', 'doParallel')) %dopar% {
  distances <- sqrt((spots_heme_1000[, 5] - spots_heme_1000[j, 5])^2 + (spots_heme_1000[, 6] - spots_heme_1000[j, 6])^2)
  
  # Inizializza i vettori per la correlazione e il p-value
  cor_before <- numeric(nrow(logratio_filtered_heme_1000))
  pval_before <- numeric(nrow(logratio_filtered_heme_1000))
  
  # Calcola gli intervalli di distanza usando quantili
  quantiles <- quantile(distances, probs = seq(0, 1, length.out = 10))  # 10 intervalli, puoi cambiare il numero
  interval_indices <- .bincode(distances, breaks = quantiles, include.lowest = TRUE)
  
  # Assegna gli intervalli di distanza ai punti di spots
  spots_heme_1000$ring <- interval_indices
  
  # Calcola la correlazione per ogni gene
  for (i in 1:nrow(logratio_filtered_heme_1000)) {
    zero <- which(is.na(logratio_filtered_heme_1000[i, ]))
    
    if (length(zero) > 0) {
      cor_before[i] <- cor(spots_heme_1000$ring[-zero], logratio_filtered_heme_1000[i, -zero], method = "spearman")
      test <- cor.test(spots_heme_1000$ring[-zero], logratio_filtered_heme_1000[i, -zero], method = "spearman")
      pval_before[i] <- test$p.value
    } else {
      cor_before[i] <- cor(spots_heme_1000$ring, logratio_filtered_heme_1000[i,], method = "spearman")
      test <- cor.test(spots_heme_1000$ring, logratio_filtered_heme_1000[i,], method = "spearman")
      pval_before[i] <- test$p.value
    }
  }
  
  # Filtra i geni in base al p-value
  significant_genes <- which(pval_before < 0.05)  # Filtra per p-value < 0.05
  
  # Calcola la correlazione assoluta solo per i geni significativi
  median_cor <- median(abs(cor_before[significant_genes]))
  
  # Restituisci il risultato per il punto j
  data.frame(arrayrow = spots_heme_1000[j, 5], arraycolumn = spots_heme_1000[j, 6], median_cor = median_cor)
}

# Arresta il cluster parallelo
stopCluster(cl)

# Crea un dataframe per i dati della mappa del tessuto e della correlazione media
map_data <- data.frame(
  x = spots_heme_1000$pixel_col,
  y = spots_heme_1000$pixel_row,
  median_correlation = ring_median_heme_1000$median_cor
)

# Crea il grafico
p <- ggplot(map_data, aes(x = x, y = y, fill = median_correlation)) +
  geom_point(shape = 21, size = 3) +
  scale_fill_gradientn(
    colors = c("blue", "cyan", "lightgreen", "gold", "darkorange", "red"), 
    limits = c(0.1, 0.5),
    oob = scales::squish
  ) +
  labs(x = "Colonna", y = "Riga", fill = "Mediana della Correlazione") +
  theme_minimal() +
  ylim(max(map_data$y), min(map_data$y))

# Salva il grafico con le dimensioni specificate
ggsave("main 10 ring cor heme_1000.png", plot = p, width = 600 / 72, height = 531 / 72, units = "in")

p_no_legend <- p + theme(legend.position = "none")

# Salva il grafico senza la legenda
ggsave("main 10 ring cor heme_1000_senza_legenda.png", plot = p_no_legend, width = 600 / 72, height = 531 / 72, units = "in", dpi = 300)


#PROVA RINGS
for (j in 1697) {
  distances <- sqrt((spots_heme_1000[, 5] - spots_heme_1000[j, 5])^2 + (spots_heme_1000[, 6] - spots_heme_1000[j, 6])^2)
}

quantiles <- quantile(distances, probs = seq(0, 1, length.out = 10))  # 10 intervalli, puoi cambiare il numero
interval_indices <- .bincode(distances, breaks = quantiles, include.lowest = TRUE)

# Assegna gli intervalli di distanza ai punti di spots
spots_heme_1000$ring <- interval_indices

map_data <- data.frame(
  x = spots_heme_1000$pixel_col,
  y = spots_heme_1000$pixel_row,
  ring = spots_heme_1000$ring
)

ggplot(map_data, aes(x = x, y = - y, color = factor(ring))) +
  geom_point(size = 2) +
  scale_color_manual(values = c(
    "1" = "#d73027",  # rosso
    "2" = "#f46d43",  # arancione-rosso
    "3" = "#fdae61",  # arancione
    "4" = "#fee08b",  # giallo
    "5" = "#d9ef8b",  # giallo-verde
    "6" = "#a6d96a",  # verde chiaro
    "7" = "#66bd63",  # verde
    "8" = "#3288bd",  # blu
    "9" = "#4575b4"   # blu scuro
  )) +
  theme_minimal() +
  labs(title = "Distances in rings",
       x = "Coordinata X",
       y = "Coordinata Y",
       color = "Ring")

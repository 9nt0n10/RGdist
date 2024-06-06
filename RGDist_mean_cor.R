---
  title: "Analysis - Heme 1000"
output: html_notebook
---
  
  # Data Preprocessing
  
  ```{r}

library(SpatialExperiment)
library(rtracklayer)
library(lobstr)
library(readr)
library(scater)
library(ggspavis)
library(sva)
library(plyr)

```

## Import data as SpatialExperiment object

```{r}

spe <- SpatialExperiment::read10xVisium(
  samples = "C:/Users/Antonio/Laboratorio/R_wd/SpatialTranscriptomics_Unisannio/SpatialTranscriptomics/samples/heme_1000/outs/",
  sample_id = "heme_1000",
  type = c("HDF5", "sparse"),
  data = "filtered",
  images = "lowres", 
  load = TRUE
)

```


Preso da SpatialLIBD (vedere bene!)

```{r}
# continuous variable containing total number of counts for each sample prior to filtering any genes
spe$sum_umi <- colSums(counts(spe))

# continuous variable containing the number of genes that have at least 1 count
spe$sum_gene <- colSums(counts(spe) > 0)
```


Reference genome alignment 

```{r}
# gene annotation - reference genome GRCm38.p6 http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.annotation.gtf.gz

gtf <-                                           
  rtracklayer::import(
    "C:/Users/Antonio/Laboratorio/R_wd/SpatialTranscriptomics_Unisannio/SpatialTranscriptomics/gencode.vM23.annotation/gencode.vM23.annotation.gtf"
  )

gtf <- gtf[gtf$type == "gene"]                   
gtf$gene_id <- gsub("\\..*", "", gtf$gene_id)    
names(gtf) <- gtf$gene_id  

match_genes <- match(rownames(spe), gtf$gene_id) 
table(is.na(match_genes))
spe <- spe[!is.na(match_genes), ]                
match_genes <- match_genes[!is.na(match_genes)]
mcols(gtf) <- mcols(gtf)[, c("source", "type", "gene_id", "gene_name", "gene_type")]  
rowRanges(spe) <- gtf[match_genes]               

rowData(spe)$gene_search <- paste0(
  rowData(spe)$gene_name, "; ", rowData(spe)$gene_id
)

```


## Anatomy info and injection site coords 

```{r}

spots <- read_csv("C:/Users/Antonio/Laboratorio/R_wd/SpatialTranscriptomics_Unisannio/SpatialTranscriptomics/samples/heme_1000/outs/spatial/tissue_positions_list.csv", col_names = FALSE)
names(spots) <- c("Barcode", "in_tissue", "arrary_row", "array_col", "pixel_row", "pixel_col")
spots <- subset(spots, in_tissue == 1) # here we remove spots not covered by tissue


anatomy <- merge(read_csv("C:/Users/Antonio/Laboratorio/R_wd/SpatialTranscriptomics_Unisannio/SpatialTranscriptomics/samples/heme_1000/heme_1000_anatomy.csv"),
                 read_csv("C:/Users/Antonio/Laboratorio/R_wd/SpatialTranscriptomics_Unisannio/SpatialTranscriptomics/samples/heme_1000/heme_1000_injection_site.csv"),
                 by = "Barcode")

spots <- merge(spots, anatomy, by = "Barcode")


# Add info to spe object
spe$anatomy <- spots$anatomy
spe$inj_site <- spots$injection_site

```


## Quality Control and data filtering

Mitochondrial genes

```{r}
# identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
table(is_mito)

# calculate per-spot QC metrics and store in colData
spe <- addPerCellQC(spe, subsets = list(mito = is_mito))

# histogram of numbers of expressed genes
hist(colData(spe)$detected, breaks = 20)

# Remove mitochondrial genes
spe <- spe[!is_mito, ]
```


NAs spots: spots that are not assigned to any anatomical cluster. This clustering was manually made by the dataset authors, some spots were not assigned to any cluster and since this information is crucial for our further analysis, we decided to remove them. 

```{r}
# we have some NAs in our anatomy info so we need to remove them
NA_spot <- c(which(is.na(spots[,"anatomy"])))  
spe <- spe[, -NA_spot, drop = FALSE]
spots <- spots[-NA_spot, ]

```


We also need to remove genes that are not expressed at all across the whole tissue. 

```{r}
# remove not expressed genes
no_expr <- which(rowSums(counts(spe)) == 0)           
length(no_expr) / nrow(spe) * 100  # percentage of not expressed genes                     
spe <- spe[-no_expr, , drop = FALSE]
```


```{r}
summary(spe$sum)

if (any(spe$sum == 0)) {
  spots_no_counts <- which(spe$sum == 0)
  ## Number of spots with no counts
  print(length(spots_no_counts))
  ## Percent of spots with no counts
  print(length(spots_no_counts) / ncol(spe) * 100)
  spe <- spe[, -spots_no_counts, drop = FALSE]
}
```

In order to compute correlations we need to have at least three counts for each gene. We consider as no relevant all those genes that have less than 3 counts.

```{r}
# Check the number of spots in which a gene is expressed and removal of genes expressed in less than 3 spots

h1000_counts <- as.matrix(assays(spe)$counts)

no_rel <- c()
n <- 0

for(i in 1:nrow(h1000_counts)){
  for(j in 1:ncol(h1000_counts)){
    if(h1000_counts[i,j] != 0){
      n = n + 1
    }
  }
  if(n <= 2){
    no_rel = c(no_rel, i)
  }
  n = 0
}

remove(h1000_counts)
spe <- spe[-no_rel, , drop = FALSE]

```


## Normalization of counts

In our object spe the only available assay is "counts" (assayNames(spe)). Normalization is an essential step in an RNA-Seq analysis, in which the read count matrix is transformed to allow for meaningful comparisons of counts across samples.

```{r}
spe <- scuttle::logNormCounts(spe)
```

h1000_counts <- as.matrix(assays(spe)$counts)


## Dimensionality Reduction: UMAP 

```{r}
# No corrected data UMAP
spe <- runUMAP(spe, 
               exprs_values = "logcounts", 
               name = "UMAP_noBC")

plotReducedDim(object = spe, 
               dimred = "UMAP_noBC", 
               colour_by = "anatomy")


## Distances from injection site

```{r}
gene_table <- data.frame(gtf$gene_id, gtf$gene_name)

rownames(spots) <- NULL

h1000_logcounts <- as.matrix(assays(spe)$logcounts)

library(foreach)
library(doParallel)

# Imposta il numero di core da utilizzare (puoi impostarlo a qualsiasi valore in base alla tua CPU)
num_cores <- 50  # Ad esempio, utilizza 4 core

# Inizializza il backend parallelo
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Inizializza i vettori risultato
# spot_mean <- data.frame(arrayrow = numeric(), arraycolumn = numeric(), mean_cor = numeric())

# Iterazione su tutti i punti
spot_mean <- foreach(j = 1:nrow(spots), .combine = rbind, .packages = c('foreach', 'doParallel')) %dopar% {
  distances <- sqrt((spots[, 5] - spots[j, 5])^2 + (spots[, 6] - spots[j, 6])^2)
  
  # Inizializza i vettori per la correlazione e il p-value
  cor_before <- numeric(nrow(h1000_logcounts))
  pval_before <- numeric(nrow(h1000_logcounts))
  
  # Iterazione su tutti i geni
  for (i in 1:nrow(h1000_logcounts)) {
    zero <- which(h1000_logcounts[i,] == 0)
    if (length(zero) > 0) {
      cor_before[i] <- cor(distances[-zero], h1000_logcounts[i, -zero], method = "spearman")
      test <- cor.test(distances[-zero], h1000_logcounts[i, -zero], method = "spearman")
      pval_before[i] <- test$p.value
    } else {
      cor_before[i] <- cor(distances, h1000_logcounts[i,], method = "spearman")
      test <- cor.test(distances, h1000_logcounts[i,], method = "spearman")
      pval_before[i] <- test$p.value
    }
  }
  
  # Filtra i geni in base al p-value
  significant_genes <- which(pval_before < 0.05)  # Filtra per p-value < 0.05
  
  # Calcola la correlazione assoluta solo per i geni significativi
  mean_cor <- mean(abs(cor_before[significant_genes]))
  
  # Restituisci il risultato
  data.frame(arrayrow = spots[j, 5], arraycolumn = spots[j, 6], mean_cor = mean_cor)
}

# Arresta il cluster parallelo
stopCluster(cl)

# Stampa i risultati
print(spot_mean)


# Crea un dataframe per i dati della mappa del tessuto e della correlazione media
map_data <- data.frame(
  x = spots$pixel_col,
  y = spots$pixel_row,
  mean_correlation = spots_mean$mean_cor
)

# Crea il grafico
ggplot(map_data, aes(x = x, y = y, fill = mean_correlation)) +
  geom_point(shape = 21, size = 3) +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(x = "Colonna", y = "Riga", fill = "Correlazione Media") +
  theme_minimal() +
  ylim(max(map_data$y), min(map_data$y))

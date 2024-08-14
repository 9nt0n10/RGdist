library(SpatialExperiment)
library(rtracklayer)
library(lobstr)
library(readr)
library(scater)
library(ggspavis)
library(sva)
library(plyr)



## Import data as SpatialExperiment object



spe <- SpatialExperiment::read10xVisium(
  samples = "C:/Users/Antonio/Laboratorio/R_wd/SpatialTranscriptomics_Unisannio/SpatialTranscriptomics/samples/heme_1000/outs/",
  sample_id = "heme_1000",
  type = c("HDF5", "sparse"),
  data = "filtered",
  images = "lowres", 
  load = TRUE
)





# continuous variable containing total number of counts for each sample prior to filtering any genes
spe$sum_umi <- colSums(counts(spe))

# continuous variable containing the number of genes that have at least 1 count
spe$sum_gene <- colSums(counts(spe) > 0)



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




## Anatomy info and injection site coords 



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




## Quality Control and data filtering


# identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
table(is_mito)

# calculate per-spot QC metrics and store in colData
spe <- addPerCellQC(spe, subsets = list(mito = is_mito))

# histogram of numbers of expressed genes
hist(colData(spe)$detected, breaks = 20)

# Remove mitochondrial genes
spe <- spe[!is_mito, ]




# we have some NAs in our anatomy info so we need to remove them
NA_spot <- c(which(is.na(spots[,"anatomy"])))  
spe <- spe[, -NA_spot, drop = FALSE]
spots <- spots[-NA_spot, ]



# remove not expressed genes
no_expr <- which(rowSums(counts(spe)) == 0)           
length(no_expr) / nrow(spe) * 100  # percentage of not expressed genes                     
spe <- spe[-no_expr, , drop = FALSE]




summary(spe$sum)

if (any(spe$sum == 0)) {
  spots_no_counts <- which(spe$sum == 0)
  ## Number of spots with no counts
  print(length(spots_no_counts))
  ## Percent of spots with no counts
  print(length(spots_no_counts) / ncol(spe) * 100)
  spe <- spe[, -spots_no_counts, drop = FALSE]
}



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




## Normalization of counts


spe <- scuttle::logNormCounts(spe)


h1000_counts <- as.matrix(assays(spe)$counts)


## Dimensionality Reduction: UMAP 


# No corrected data UMAP
spe <- runUMAP(spe, 
               exprs_values = "logcounts", 
               name = "UMAP_noBC")

plotReducedDim(object = spe, 
               dimred = "UMAP_noBC", 
               colour_by = "anatomy")


## Distances from injection site


gene_table <- data.frame(gtf$gene_id, gtf$gene_name)

rownames(spots) <- NULL

h1000_logcounts <- as.matrix(assays(spe)$logcounts)

# Inizializzazione delle variabili
cor_before <- numeric(nrow(spots))
pval_before <- numeric(nrow(spots))
spot_median_gs <- data.frame(arrayrow = numeric(), arraycolumn = numeric(), median_rank = numeric())

# Iterazione su tutti i punti
for (j in 1697) {
  distances <- sqrt((spots[, 5] - spots[j, 5])^2 + (spots[, 6] - spots[j, 6])^2)
  
  # Calcola gli intervalli di distanza usando quantili
  quantiles <- quantile(distances, probs = seq(0, 1, length.out = 10))  # 10 intervalli, puoi cambiare il numero
  interval_indices <- .bincode(distances, breaks = quantiles, include.lowest = TRUE)
  
  # Assegna gli intervalli di distanza ai punti di spots
  spots$ring <- interval_indices
  
  # Iterazione su tutti i geni
  for (i in 1:nrow(h1000_logcounts)) {
    zero <- which(h1000_logcounts[i,] == 0)
    if (length(zero) > 0) {
      cor_before[i] <- cor(spots$ring[-zero], h1000_logcounts[i, -zero], method = "spearman")
      test <- cor.test(spots$ring[-zero], h1000_logcounts[i, -zero], method = "spearman")
      pval_before[i] <- test$p.value
    } else {
      cor_before[i] <- cor(spots$ring, h1000_logcounts[i,], method = "spearman")
      test <- cor.test(spots$ring, h1000_logcounts[i,], method = "spearman")
      pval_before[i] <- test$p.value
    }
  }
}

names(cor_before) <- rownames(rowData(spe))
names(pval_before) <- rownames(rowData(spe))
padjust_bf <- p.adjust(pval_before, method="fdr")

rowData(spe)$cor_before <- cor_before
rowData(spe)$pvalue_before <- pval_before
rowData(spe)$pvalue_fdr_before <- padjust_bf

rowData(spe)$rank_before <- rank(padjust_bf, ties.method = "min")

cor_genes_bf <- which(cor_before >= 0.5 & padjust_bf < 0.01)
anticor_genes_bf <- which(cor_before <= -0.5 & padjust_bf < 0.01)

col <- rep("gray", nrow(rowData(spe)))
col[cor_genes_bf] <- "pink"
col[anticor_genes_bf] <- "green"

gene_name <- rowData(spe)$gene_name

plot(x = cor_before,
     y = -log(padjust_bf),
     xlab = "Correlation",
     ylab = "-log p.values",
     type = "p",
     col = col,
     pch = 20,
     main = "Before batch correction"
)
abline(v=-0.5, col="black", lty=2)
abline(v=0.5, col="black", lty=2)
abline(h=-log(0.01), col="black", lty=2)
x = c(cor_before[anticor_genes_bf], cor_before[cor_genes_bf])
y = c(padjust_bf[anticor_genes_bf], padjust_bf[cor_genes_bf])
text(x, -log(y), 
     labels = c(gene_name[anticor_genes_bf], gene_name[cor_genes_bf]),
     cex = 0.4)

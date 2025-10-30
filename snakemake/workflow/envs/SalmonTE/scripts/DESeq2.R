library(DESeq2)
library(readr)
library(dplyr)
library(ggplot2)
library(tibble)
library(BiocManager)
library(EnhancedVolcano)
library(tidyverse)
library(ggpubr) 
library(ggrepel)
library(ashr)
library(MASS) 
library(RColorBrewer)
library(pheatmap)
library(PoiClaClu) 
library(openxlsx)

args <- commandArgs(trailingOnly = TRUE)
meta <- args[1]
contrast <- args[2]
count <- args[3]
outdir <- args[4]
max_fdr <- as.numeric(args[5])
min_lfc <- as.numeric(args[6])
indfilter <- as.logical(args[7])
cookscutoff <- as.logical(args[8])
blackSamples <- args[9]
anno <- args[10]
nclass <- 4

# Prepare output dir
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

meta.df <- read_delim(meta, delim =",", comment = '#') 
meta.df <- meta.df[!(meta.df[, 1] %in% blackSamples), ]
contrast.df <- read_delim(contrast, delim =",", comment = '#') 
df <- read_delim(count, delim =",", comment = '#')   
names(df)[1] <- "GENE_ID"
df <- df %>% mutate(GENE_ID = gsub("\\.[^.]*$", "", GENE_ID))
if (sum(duplicated(df$GENE_ID)) > 0) {
  warning("count table gene_id not unique, has made it unique to run, but duplicated gene_id not valid for further analysis, please manually check duplicated gene_id")
  df$GENE_ID <- make.unique(df$GENE_ID)
  

}
df <- column_to_rownames(df, var =  colnames(df)[1]) 
df[, 1:ncol(df)] <- sapply(df[, 1:ncol(df)], as.integer)
cts <- as.matrix(df[, 1:ncol(df)])

# Skip the filtering of cts step
group_count <- plyr::count(meta.df[[2]])
meta.df <- meta.df[match(colnames(cts), meta.df[[1]]), ]
dim(cts)

# Design matrix
meta_idx <- match(colnames(cts), meta.df[[1]])
meta_idx <- meta_idx[!is.na(meta_idx)]
meta.df <- meta.df[meta_idx, ]
cts <- cts[, colnames(cts) %in% meta.df[[1]]]

sample <- factor(meta.df[[1]])
batch <- factor(meta.df[[3]])
group <- factor(meta.df[[2]])

coldata <- data.frame(row.names=colnames(cts), 
                      sample,
                      group,
                      batch)
coldata

# Module fitting
if (length(levels(batch)) > 1){
  dds <- DESeqDataSetFromMatrix(countData = cts, 
                                colData = coldata, 
                                design = ~  0 + group + batch)
  } else {
    dds <- DESeqDataSetFromMatrix(countData = cts, 
                                  colData = coldata, 
                                  design = ~  0 + group)
    }
dds
dds_res <- try(DESeq(dds), silent = TRUE)
if (inherits(dds_res, "try-error")) {
  cat("using manual dispersion estimation\n")
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersionsGeneEst(dds)
  dispersions(dds) <- mcols(dds)$dispGeneEst
  dds <- nbinomWaldTest(dds)
} else {
  dds <- dds_res
}
resultsNames(dds)

# Save DESeq2 normalized counts
normalized_counts <- counts(dds, normalized = TRUE)
normalized_counts <- data.frame(normalized_counts)
write.csv(normalized_counts, file.path(outdir, 'DESeq2NormalizedCounts.csv'))
print("saved in DESeq2NormalizedCounts.csv")

# Save dispersion plot
pdf(file.path(outdir, "dispersion_plot.pdf"))
plotDispEsts(dds, main = "Dispersion plot")
dev.off()

# Save PCA/heatmap plot:
plotQC_PCA <- function(dds, outdir = NULL, fname = NULL) {
  vsd <- varianceStabilizingTransformation(dds) # fixed num < num(rowsum>5)
  pcaData <- plotPCA(vsd, intgroup = 'group', returnData=TRUE) # labeling fixed
  percentVar <- round(100 * attr(pcaData, 'percentVar'), 1)
  if (length(levels(batch)) > 1){
    ggplot(pcaData, aes(PC1, PC2, color = group, shape = batch)) +
      geom_point(size = 3) +
      xlab(paste0("PC1: ", percentVar[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar[2], "% variance")) +
      geom_label_repel(aes(label = sample),
                       box.padding = 0.35,
                       point.padding = 1,
                       segment.color = 'grey50',
                       segment.alpha = 0.5,
                       show.legend = FALSE) + # if TRUE, legend display might not be correct
      theme_classic()
    ggsave(file.path(outdir, fname), width = 10, height = 10)
  }else{
    ggplot(pcaData, aes(PC1, PC2, color = group)) +
      geom_point(size = 3) +
      xlab(paste0("PC1: ", percentVar[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar[2], "% variance")) +
      geom_label_repel(aes(label = sample),
                       box.padding = 0.35,
                       point.padding = 1,
                       segment.color = 'grey50',
                       segment.alpha = 0.5,
                       max.overlaps = 20,
                       show.legend = FALSE) + # if TRUE, legend display might not be correct
      theme_classic()
    ggsave(file.path(outdir, fname), width = 10, height = 10)
  }
}
if (max(unlist(lapply(meta.df[[1]], nchar))) <= 36){
  plotQC_PCA(dds, outdir = outdir, "sample_PCA.labeled.pdf")
}

PoissonDistanceHeatmap <-  function(){
  # dist matrix
  poisd <- PoissonDistance(t(counts(dds)))
  samplePoisDistMatrix <- as.matrix( poisd$dd ) 
  rownames(samplePoisDistMatrix) <- coldata$sample
  colnames(samplePoisDistMatrix) <- NULL
  
  # sample_col
  my_sample_col <- data.frame(group = coldata[,'group'])
  row.names(my_sample_col) <- coldata[,'sample']
  
  pheatmap(samplePoisDistMatrix,
           annotation_row = my_sample_col,
           clustering_distance_rows=poisd$dd,
           clustering_distance_cols=poisd$dd,
           clustering_method='complete', 
           legend = T, 
           col=colorRampPalette( rev(brewer.pal(9, "Blues")) )(255), 
           filename = file.path(outdir, "sample_poisson_distance.pdf"))
}
PoissonDistanceHeatmap()

# DE analysis
zscore <- function(matrix){
  return( t(scale(t(matrix))))
}
rename_num_vector_by_order <- function(l){
  # l have to be a vector of numbers
  # output vector of roman numbers ordered by appearance in the input vector
  # e.g. c(2,3,3,2,1) -> c(I, II, II, I, III)
  # test rename_num_vector_by_order(c(2,3,3,2,1))
  u <- unique(l)
  n=0
  for (i in u){
    n = n+1; 
    l <- replace(l, l==i, as.character(as.roman(n)))
  }
  return(l)
}
Heatmap <- function(df, nclass=2, outdir = NULL, fname="heatmap", main="title"){
  num_sig <- dim(df)[1]
  
  if (num_sig > 10000){
    warning("sorry, too many significant genes for heatmap, heatmap skipped")
    warning(num_sig)
    return (0)
  }
  
  if (num_sig < 20){CELL_HEIGHT = 20}else{CELL_HEIGHT = NA} # for plotting
  
  # prep group label
  my_sample_col <- data.frame(group = coldata[,'group'])
  row.names(my_sample_col) <- coldata[,'sample']
  
  if (num_sig < 1){
    print(paste("No sig DEG to create Heatmap for:", main) )
    return (0)
  }
  
  if (num_sig <= nclass){
    print(paste("Num-DEG less than num-class, there is no cut-tree for:", main))
    plot1 <- pheatmap(df,
                      annotation_col = my_sample_col,
                      main = main,
                      cluster_cols = F,
                      cluster_rows = F,
                      border_color = NA,
                      show_rownames = F,
                      breaks = seq(-3, 3, length.out = 100),
                      cellheight = CELL_HEIGHT,
                      filename = file.path(outdir, paste(fname, "pdf", sep = ".")))
    plot2 <- pheatmap(df,
                      annotation_col = my_sample_col,
                      main = main,
                      cluster_cols = T,
                      cluster_rows = F,
                      border_color = NA,
                      show_rownames = F,
                      breaks = seq(-3, 3, length.out = 100),
                      cellheight = CELL_HEIGHT,
                      filename = file.path(outdir, paste(fname, "v2.pdf", sep = ".")))
  }
  
  if (num_sig > nclass){
    p <- pheatmap(df, cluster_cols = F, cutree_rows = nclass) # Pre-plot to get hclust for cut-tree
    # cutree manually, print out gene classification (https://www.biostars.org/p/287512/)
    gene_classes <- sort(cutree(p$tree_row, k=nclass))
    gene_classes <- data.frame(gene_classes) 
    classification <- merge(gene_classes, df, by = 0)
    row.names(classification) <- classification$Row.names
    # Re-order original data (genes) to match ordering in heatmap (top-to-bottom)
    idx <- rownames(df[p$tree_row[["order"]],])
    classification <- classification[idx,] 
    # rename gene classes
    classification$gene_classes <- rename_num_vector_by_order(classification$gene_classes)
    write.xlsx(classification, file.path(outdir, paste(fname,"gene_class.xlsx", sep = ".")), 
               keepNA = T, na.string = 'NA')
    # get class label
    annotation_row = data.frame(class = classification$gene_classes)
    row.names(annotation_row) <- row.names(classification)
    annotation_row$class <- as.character(annotation_row$class)
    # output final plot
    plot1 <- pheatmap(df,
                      annotation_col = my_sample_col,
                      main = main,
                      border_color = NA,
                      cluster_cols = F,
                      cutree_rows = nclass,
                      breaks = seq(-3, 3, length.out = 100),
                      show_rownames = F,
                      annotation_row = annotation_row,
                      cellheight = CELL_HEIGHT,
                      filename = file.path(outdir, paste(fname, "pdf", sep = ".")))
    
    plot2 <- pheatmap(df,  # sample clustering
                      annotation_col = my_sample_col,
                      main = main,
                      border_color = NA,
                      cluster_cols = T,
                      cutree_rows = nclass,
                      breaks = seq(-3, 3, length.out = 100),
                      show_rownames = F,
                      annotation_row = annotation_row,
                      cellheight = CELL_HEIGHT,
                      filename = file.path(outdir, paste(fname, "v2.pdf", sep = ".")))
  }
}
parse_name <- function(name){
  name <- gsub(" ", "", name)
  name <- gsub(";$", "", name)
  names <- strsplit(name, ";") [[1]]
  name <- gsub(";", ".", name)
  return (list(name=name, names=names))
}
process_deseq_res <- function(res = "lfcshrink.res", 
                              res2 = "results.res", 
                              outdir = NULL,
                              name = 'name', 
                              anno='anno.df', 
                              norm_exp = "tpm.df", 
                              nclass = 4) {
  sig_idx <- res$padj<max_fdr & abs(res$log2FoldChange) > min_lfc
  sig_idx[is.na(sig_idx)] <- FALSE
  res_sig <- res[sig_idx,]
  print(table(sig_idx))
  
  up_idx <- res$padj<max_fdr & res$log2FoldChange > min_lfc
  up_idx[is.na(up_idx)] <- FALSE
  res_sig <- res[up_idx,]
  print(table(up_idx))
  
  down_idx <- res$padj<max_fdr & res$log2FoldChange < -min_lfc
  down_idx[is.na(down_idx)] <- FALSE
  res_sig <- res[down_idx,]
  print(table(down_idx))
  
  res.df <- as.data.frame(res)
  names(res.df)[2] <- "log2FoldChange_shrunken"
  names(res.df)[3] <- "lfcSE_shrunken"
  
  res2.df <- as.data.frame(res2)
  names(res2.df)[2] <- "log2FoldChange_raw"
  names(res2.df)[3] <- "lfcSE_raw"
  res2.df <- res2.df[, c(2, 3)]
  
  resdata <- merge(res.df, res2.df, by=0, sort=F, all.x=T)
  resdata <- merge(resdata, norm_exp, by.x=1, by.y=0, all.x=T, sort=F)
  head(resdata)
  sig_idx <- resdata$padj<max_fdr & abs(resdata$log2FoldChange_shrunken) > min_lfc # important to put this line right before output sig.xlsx
  sig_idx[is.na(sig_idx)] <- FALSE
  resdata.sig <- resdata[sig_idx,]
  head(resdata.sig)
  
  write.csv(resdata, file.path(outdir, paste(name, 'deseq2.csv', sep = '.')))
  write.csv(resdata.sig, 
             file.path(paste(name, 'deseq2.sig.FDR', max_fdr, 'LFC', min_lfc, 'csv', sep = '.')))
  
  ##  Plots
  ggplot(data.frame(res), aes(x=pvalue))+
    geom_histogram(color="darkblue", fill="lightblue")
  ggsave(file.path(paste0(name, '.pvalue.pdf')))
  
  ggplot(data.frame(res), aes(x=padj))+
    geom_histogram(color="darkblue", fill="lightblue")
  ggsave(file.path(paste0(name, '.fdr.pdf')))
  
  # Skip the following plots as there is no matching annotation file for TE
  # maplot(res, anno, paste0(name, '.maplot.shrunken_lfc.pdf'))
  # maplot(res2, anno, paste0(name, '.maplot.raw_lfc.pdf'))
  # 
  # volcanoplot(res, anno, lfcthresh=min_lfc, sigthresh=max_fdr,
  #             textcx=.8,  name= paste(name, "LFC_shrunken", sep="."))
  # volcanoplot(res2, anno, lfcthresh=min_lfc, sigthresh=max_fdr,
  #             textcx=.8,name= paste(name, "LFC_raw", sep="."))
  
  n1 <- dim(resdata.sig)[2]
  n2 <- dim(norm_exp)[2]
  zscore.df <- zscore(resdata.sig[, (n1-n2+1):n1])
  rownames(zscore.df) <- resdata.sig[,1]
  colnames(zscore.df) <- gsub (":TPM", "", colnames(zscore.df))
  colnames(zscore.df) <- gsub (":DESeq2NormalizedCount", "", colnames(zscore.df))
  Heatmap(zscore.df, nclass = nclass,
          outdir = outdir,
          fname = paste(name, "heatmap", sep="."),
          main = paste(name, "LFC >", min_lfc, "FDR <", max_fdr ))
}

for (i in 1:dim(contrast.df)[2]){
  name1 <- parse_name(contrast.df[1,i])
  name2 <- parse_name(contrast.df[2, i])
  name <- paste(name1$name, name2$name, sep = "_vs_")
  if (nchar(name) > 100) {name = paste0('contrast', i)}
  print(paste(">>>", i, name))
  
  # check group_count for all groups in each contrast
  groups <- unique(c(name1$names, name2$names))
  with_sample <- groups %in% group_count[[1]]  # count >= 1
  if (!all(with_sample)){
    print(paste('some groups in', name, 'has no filtered samples, skipped !!!'))
    print(group_count)
    next
  }
  group_with_reps <-  group_count[[1]][group_count$freq >= 2]
  with_reps1 <- any(name1$names %in% group_with_reps)
  with_reps2 <- any(name2$names %in% group_with_reps)
  if (!any(with_reps1, with_reps2)) {
    print(paste('no groups in', name, 'has reps, skipped !!!'))
    print(group_count)
    next
  }
  gsub("group", "", resultsNames(dds))
  poss <- match(name1$names, gsub("group", "", resultsNames(dds)))
  negs <- match(name2$names, gsub("group", "", resultsNames(dds)))
  contrast <- rep(0, length(resultsNames(dds)))
  contrast[poss] <- 1/length(poss)
  contrast[negs] <- -1/length(negs)
  print(data.frame(resNames=gsub("group", "", resultsNames(dds)), 
                   contrast=contrast))
  
  res <- lfcShrink(dds, contrast = contrast, type = 'ashr')
  res2 <- results(dds, contrast = contrast, 
                  independentFilter = indfilter, cooksCutoff = cookscutoff)
  
  process_deseq_res(res = res, res2 = res2, outdir = outdir, name = name, anno = anno, norm_exp = normalized_counts, nclass = nclass)
}
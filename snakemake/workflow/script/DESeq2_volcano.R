# Standalone volcano plot from a pre-computed DESeq2 results table (no DESeq2 re-run)
# Search for "NEED EDIT" and update as needed
# Adapted from DESeq2.Rmd and ScienceMachine

setwd("/Users/kaihu/Projects/AvijitMallick/RNAseq_06042026")
library(EnhancedVolcano)
library(tidyverse)
library(ggplot2)
library(openxlsx)
library(readxl)

gene_symbols <- sapply(df[[3]], function(x) {
  parts <- strsplit(x, ";")[[1]]
  if (length(parts) >= 2) {
    return(trimws(parts[2]))
  } else {
    return(NA_character_)
  }
})

readTable <- function(fname){
  if (str_ends(fname, '.txt')){
    df <- read_delim(fname, delim = "\t", comment = '#')
  } else if (str_ends(fname, '.csv')){
    df <- read_delim(fname, delim = ",", comment = '#')
  } else if (str_ends(fname, '.xlsx')){
    df <- read.xlsx(fname, na.strings = 'NA', sheet = 1)  # na term important
  } else {
    stop("file format not supported")
  }
  return(df)
}

volcanoplot <- function(res, name = 'name',
                        lfc_col = "log2FoldChange", padj_col = "padj", name_col = NULL,
                        lfcthresh = 0.585, sigthresh = 0.05,
                        xlim = 100, ylim = 1000) {
  res <- data.frame(res, check.names = FALSE)
  
  # gene-label column: use name_col, else first non-numeric column, else row names
  if (is.null(name_col)) {
    char_cols <- names(res)[sapply(res, function(x) is.character(x) || is.factor(x))]
    name_col <- if (length(char_cols) > 0) char_cols[1] else NA
  }
  if (is.na(name_col) || !(name_col %in% names(res))) {
    res$.label <- rownames(res)
    name_col <- ".label"
    warning("No gene-name column found, using row names as labels")
  }
  
  # remove NA padj
  res[[padj_col]][is.na(res[[padj_col]])] <- 1
  # cap x, y to set plot limits
  res[[padj_col]][res[[padj_col]] < 10^(-ylim)] <- 10^(-ylim)
  res[[lfc_col]][res[[lfc_col]] >  xlim] <-  xlim
  res[[lfc_col]][res[[lfc_col]] < -xlim] <- -xlim
  
  # count up / down significant genes
  pos.n <- sum(res[[padj_col]] < sigthresh & res[[lfc_col]] >  lfcthresh)
  neg.n <- sum(res[[padj_col]] < sigthresh & res[[lfc_col]] < -lfcthresh)
  
  EnhancedVolcano(res,
                  lab = res[[name_col]],
                  x = lfc_col,
                  y = padj_col,
                  title = "",
                  subtitle = paste0("Up:", pos.n, ", Down:", neg.n),
                  xlab = bquote(~Log[2]~ "Fold Change"),
                  ylab = bquote(~-Log[10]~italic(FDR)),
                  pCutoff = max_fdr,
                  FCcutoff = min_lfc,
                  cutoffLineType = 'twodash',
                  cutoffLineWidth = 0.8,
                  legendLabels = c('NS', expression(Log[2]~FC),
                                   "FDR", expression(FDR~and~Log[2]~FC)),
                  caption = paste0('Total = ', nrow(res), ' genes'),
                  legendPosition = 'right',
                  legendLabSize = 10,
                  axisLabSize = 10,
                  legendIconSize = 3.0)
  ggsave(paste(name, "pdf", sep = "."), width = 8, height = 6)
}


# Parameters
max_fdr <- 0.05
min_lfc <- 0.585
xlim    <- 100   # cap |log2FoldChange| at this value
ylim    <- 1000  # cap -log10(padj) at this value
lfc_col  <- "log2FoldChange_shrunken"
padj_col <- "padj"
name_col <- "Name"

# Pre-computed DE results table (OneStopRNAseq *.deseq2.xlsx / .csv / .txt)
# deFile <- "results/OneStopRNAseq/analysis_1/DESeq2/WT_oligomycin_vs_WT_DMSO.deseq2.xlsx" #NEED EDIT
# res <- readTable(deFile)
# 
# geneSet <- "data/from_aviji/Endosomal transport.xlsx" #NEED EDIT
# plotName <- "results/volcano/Endosomal_transport_WT_oligomycin_vs_WT_DMSO.LFC_shrunken" #NEED EDIT
# geneSet <- "data/from_aviji/extracellular structure organization.xlsx" #NEED EDIT
# plotName <- "results/volcano/Extracellular_structure_organization_WT_oligomycin_vs_WT_DMSO.LFC_shrunken" #NEED EDIT
# geneSet <- "data/from_aviji/Retrograde transport.xlsx" #NEED EDIT
# plotName <- "results/volcano/Retrograde_transport_WT_oligomycin_vs_WT_DMSO.LFC_shrunken" #NEED EDIT

deFile <- "results/OneStopRNAseq/analysis_1/DESeq2/ATF4_DMSO_vs_WT_DMSO.deseq2.xlsx" #NEED EDIT
res <- readTable(deFile)
geneSet <- "data/from_aviji/Endosomal transport.xlsx" #NEED EDIT
plotName <- "results/volcano/Endosomal_transport_ATF4_DMSO_vs_WT_DMSO.LFC_shrunken" #NEED EDIT
geneSet <- "data/from_aviji/extracellular structure organization.xlsx" #NEED EDIT
plotName <- "results/volcano/Extracellular_structure_organization_ATF4_DMSO_vs_WT_DMSO.LFC_shrunken" #NEED EDIT
geneSet <- "data/from_aviji/Retrograde transport.xlsx" #NEED EDIT
plotName <- "results/volcano/Retrograde_transport_ATF4_DMSO_vs_WT_DMSO.LFC_shrunken" #NEED EDIT

# deFile <- "results/OneStopRNAseq/analysis_1/DESeq2/ATF5_DMSO_vs_WT_DMSO.deseq2.xlsx" #NEED EDIT
# res <- readTable(deFile)
# geneSet <- "data/from_aviji/Endosomal transport.xlsx" #NEED EDIT
# plotName <- "results/volcano/Endosomal_transport_ATF5_DMSO_vs_WT_DMSO.LFC_shrunken" #NEED EDIT
# geneSet <- "data/from_aviji/extracellular structure organization.xlsx" #NEED EDIT
# plotName <- "results/volcano/Extracellular_structure_organization_ATF5_DMSO_vs_WT_DMSO.LFC_shrunken" #NEED EDIT
# geneSet <- "data/from_aviji/Retrograde transport.xlsx" #NEED EDIT
# plotName <- "results/volcano/Retrograde_transport_ATF5_DMSO_vs_WT_DMSO.LFC_shrunken" #NEED EDIT

df_geneset <- read_excel(geneSet, col_names = FALSE, sheet = 1)
gene_symbol <- sapply(df_geneset[[3]], function(x) {
  parts <- strsplit(x, ";")[[1]]
  if (length(parts) >= 2) {
    return(trimws(parts[2]))
  } else {
    return(NA_character_)
  }
})

res_subset <- res |> dplyr::filter(Name %in% gene_symbol)

volcanoplot(res_subset, name = plotName,
            lfc_col = lfc_col, padj_col = padj_col, name_col = name_col,
            lfcthresh = min_lfc, sigthresh = max_fdr,
            xlim = xlim, ylim = ylim)


#' Visualize PAM50 Prediction
#'
#' This is to generate a heatmap to check the result of PAM50 prediction.
#'
#' @param es dataframe, rows are genes columns are samples. Row names has to be EntrezID.
#' @param PAM50.pred array of each sample's PAM50 label
#'   (you can get the label using PAM50 function)
#' @param add_labels Additional argument for adding extra labels (dataframe) to show in
#'   the heatmap. Labels have to be the same sample order as PAM50.pred
#' @param cluster_rows bool (defaul = False), indicate whether to reorder the genes.
#'   Default uses Parker et al. paper's gene order.
#' @param cluster_cols bool (defaul = False), indicate whether to reorder the samples.
#' @param scale string (default = "row"). Options are "none", "row", or "column".
#'   This is to rescale data (calculate zscore) based on row, column or no rescaling.
#' @param ... other parameters passed onto pheatmap, for instance, main = "Treatment Arms".
#'
#' @import pheatmap
#' @import RColorBrewer
#'
#' @export

plotPAM50 <- function(es, PAM50.pred, add_labels = NULL, cluster_rows = F, cluster_cols = F, scale = "row", ...){

  parker_gene_lst = c("UBE2C", "PTTG1", "MYBL2", "BIRC5", "CCNB1", "TYMS", "MELK", "CEP55", "NDC80", "UBE2T",
                      "RRM2", "CDC6", "ANLN", "ORC6", "KIF2C", "EXO1", "NUF2", "CENPF", "CCNE1", "MKI67",
                      "CDC20", "MMP11", "GRB7", "ERBB2", "TMEM45B", "BAG1", "PGR", "MAPT", "NAT1", "GPR160",
                      "FOXA1", "BLVRA", "CXXC5", "ESR1", "SLC39A6", "KRT17", "KRT5", "SFRP1", "BCL2", "KRT14",
                      "MLPH", "MDM2", "FGFR4", "MYC", "MIA", "FOXC1", "ACTR3B", "PHGDH", "CDH3", "EGFR")

  df.p50 <- extract_PAM50(es)

  sidx <- sort.int(PAM50.pred, index.return=T)
  df_in_sort <- df.p50[parker_gene_lst, colnames(es)[sidx$ix]]

  annot.c <- data.frame(PAM50 = sidx$x)
  ## if additional labels are supplied
  if(!is.null(add_labels)){
    annot.c <- cbind(annot.c, add_labels[sidx$ix, ])
  }

  rownames(annot.c) <- colnames(df_in_sort)

  breaksList <- seq(-3, 3, by = 0.5)
  pheatmap(df_in_sort,
           scale = scale,
           color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(length(breaksList)),
           breaks = breaksList,
           annotation_col = annot.c,
           clustering_distance_rows = "correlation",
           clustering_distance_cols = "correlation",
           clustering_method = "ward.D2",  # complete, ward.D2, single, average
           show_colnames = F,
           cluster_rows = cluster_rows,
           cluster_cols = cluster_cols,
           fontsize_col = 8,
           ...)
}


#' Plot Validation Genes
#'
#' This plots distribution of validation genes (4 by default) for validating
#' the prediction result.
#'
#' Default 4 genes:
#'
#' "ERBB2": Her2 gene; Her2 group should have the highest ERBB2 value.
#'
#' "ESR1": ER gene; LumA and LumB groups should have the highest ESR1 value.
#'
#' "MKI67": tumor growth factor; LumA should have the lowest MKI67 value because
#' LumA patients usually have the best prognosis.
#'
#' "PGR": PR gene; LumA should have the highest PGR value.
#'
#'
#' @param es dataframe, rows are genes columns are samples. Row names has to be EntrezID.
#' @param PAM50.pred array of each sample's PAM50 label
#' @param genes (optional) vector of validation gene names
#'
#' @import tidyverse
#'
#' @export

## Plot Validation Genes for Checking the Prediction Result
plotValidationGenes <- function(es, PAM50.pred,
                                genes = c("ESR1","PGR","MKI67","ERBB2")) {

  df.p50 <- extract_PAM50(es)
  idx <- genes %in% rownames(df.p50)

  if (sum(idx) == length(genes)) {
    df.p50 <- df.p50[genes,] ## subset only those 4 genes
  } else {
    cat("Missing:", genes[!idx])
    stop("Please check the gene name.")
  }

  es.t <- data.frame(t(df.p50))
  colnames(es.t) <- rownames(df.p50)
  es.t$SampleID <- rownames(es.t)
  es.t$Subtype <- PAM50.pred

  ## create labels with count
  lvls <- levels(as.factor(PAM50.pred))
  lb.count <- table(PAM50.pred)
  labels <- paste0(lvls, " (N=", lb.count, ")")

  ## boxplot
  es.t %>%
    gather(Symbol, value, -SampleID, -Subtype) %>%
    ggplot(aes(Subtype, value, fill = Subtype)) +
    geom_boxplot() +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    facet_wrap(~Symbol, scales = "free_y") +
    scale_fill_discrete(name="Subtype", labels=labels)
}

#' Convert expression values to Rank
#'
#' @description
#' `getPAM50rank(df)` first extracts PAM50 genes based on the EntrezID. Therefore,
#' the input dataframe `df` (usually the assayData from expressionSet;
#' rows are genes, and columns are samples) has to have rownames assigned to
#' EntrezID (Gene ID). After extracting 50 genes, it will then rank those 50 genes
#' based on the values per sample, and then finally normalized the rank values
#' by the number of genes in the dataframe (usually 50 genes) .
#'
#' @param df dataframe, rows are genes columns are samples. Row names has to be EntrezID.
#' @author Andrew Chang
#' @export
#'
getPAM50rank <- function(df){

  ## rank and normalize
  df.rank <- data.frame(apply(df, 2,
                             function(x) rank(x, ties.method = "first")))/dim(df)[1]
  rownames(df.rank) <- annot$EntrezID

  return(df.rank)
}

#' PAM50 prediction via PAMR
#'
#' @description
#' After getting the rank values, data will be passed into `getPAM50` which uses
#' a trained PAMR model for prediction. This will give a single PAM50 subtype prediction
#' for each sample, it also provides the prediction probability for all 5 subtypes
#' per patient.
#'
#' @param df dataframe, rows are genes columns are samples. Row names has to be EntrezID.
#' @importFrom pamr pamr.predict
#' @export

getPAM50 <- function(df){

  threshold <- 0.3
  pamr.pred <- pamr.predict(fit, as.matrix(df),
                            threshold,
                            type = "posterior")
  pamr.pred <- data.frame(pamr.pred)
  pamr.pred$PAM50 <- as.character(pamr.predict(fit, as.matrix(df), threshold))

  return(pamr.pred)
}


#' PAM50 main function
#'
#' @description
#' This is what you need for getting PAM50 prediction for each sample. It takes
#' RNA-seq values (raw or CPM) and returns PAM50 subtype for each sample.
#'
#' @param df dataframe, rows are genes columns are samples. Row names has to be EntrezID.
#' @param cutoff numeric (0-1, default = 0), threshold value for calling the prediction.
#'     0 means no cutoff, predicted subtype is whichever subtype gives the highest probability.
#'     0.8 means the predicted subtype has to have probability > 0.8, otherwise 'NA' is the label.
#'
#' @import tidyverse
#' @export
#'
#' @examples
#' library(Biobase)
#' library(tidyverse)
#' library(pamr)
#'
#' df <- data.frame(assayData(es)$exprs)
#' df <- df[,1:10] ## select 10 samples
#' annotation <- fData(es)
#'
#' rownames(df) <- annotation$EntrezID
#'
#' df_pred <- PAM50(df)
#' knitr::kable(df_pred)

PAM50 <- function(df, cutoff = 0){

  df <- prepare_df(df)

  df.rank <- getPAM50rank(df)
  df.pred <- getPAM50(df.rank)
  rownames(df.pred) <- colnames(df)

  ## prediction probability <= cutoff will be "NA"
  df.pred <- df.pred %>%
    mutate(SampleID = rownames(.)) %>%
    gather(Subtype, Prob, -c("PAM50", "SampleID")) %>%
    group_by(SampleID) %>%
    mutate(PAM50 = ifelse(max(Prob)<=cutoff, "NA", PAM50)) %>%
    ungroup() %>%
    spread(Subtype, Prob) %>%
    arrange(match(SampleID, colnames(df))) %>% ## make sure the row order is the same as input expression sample order
    select(-SampleID)

  return(df.pred)
}

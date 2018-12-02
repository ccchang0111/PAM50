#' Extract PAM50 genes
#'
#' @param df expression dataframe, rownames has to be EntrezID (Gene ID)
#' @param an annotation dataframe. Optional, if passed in, those 50 gene info will be extracted.
#'
#' @return es.50 dataframe with extracted PAM50 gene symbols as rownames
#' @return an.50 (optional) if an is passed in, it will return an.50 dataframe that contains extracted PAM50 gene info
#'
#' @export
#'
extract_PAM50 <- function(df, an = NULL){

  #df <- prepare_df(df)
  idx_sym <- rownames(df) %in% annot$Symbol
  idx_gid <- rownames(df) %in% annot$EntrezID

  if (sum(idx_sym) == 50) {
    es.50 <- df[idx_sym,]
    idx <- idx_sym

  } else if (sum(idx_gid) == 50) {
    es.50 <- df[idx_gid,]
    idx <- idx_gid

    ## re-order the PAM50 data based on the input data's Gene ID
    annot_tmp <- data.frame(annot)
    rownames(annot_tmp) <- annot_tmp$EntrezID
    annot_tmp <- annot_tmp[rownames(es.50),]

    ## rename the rownames to Gene symbol (PAM50 genes are unique)
    rownames(es.50) <- annot_tmp$Symbol

  } else {
    cat("only", sum(idx_sym) + sum(idx_gid), "genes are found. Please check your data.")
  }

  #es.50 <- df[idx,]

  if (!is.null(an)){
    an.50 <- an[idx,]

    return(list("es50" = es.50,
                "an50" = an.50))
  } else {
    return(es.50)
  }
}


##' Check input data frame for main gPAM50 function
##'
##' @param df input data frame that is used for prediction
##' @author Paul Paczuski
prepare_df <- function(df) {

  testthat::expect_true(is.data.frame(df))

  ## the following are checks for getPAM50rank
  sym_order <- annot$Symbol
  gid_order <- annot$EntrezID

  ## check if rownames are symbol or GeneID
  ## if it is symbol, convert to GeneID
  ## if GeneID, then procede
  ## else, throw a warning

  rnames <- rownames(df)
  idx_sym <- sym_order %in% rnames
  idx_gid <- gid_order %in% rnames

  if (sum(idx_sym) == 50) {
    #message("Found all 50 genes.")

    ## convert rownames to geneid
    df <- df[sym_order,]
    rownames(df) <- gid_order

  } else if (sum(idx_gid) == 50) {
    #message("Found all 50 genes.")
    df <- df[gid_order,]

  } else {
    warning("Missing Gene IDs:", gid_order[!idx_gid])
    warning("Missing Gene Symbols:", annot$Symbol[!idx_gid])
    warning("This may cause problem in prediction.")
    stop()
  }

  return(df)
}


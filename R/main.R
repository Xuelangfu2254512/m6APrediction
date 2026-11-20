#' Encode DNA sequence into features
#'
#' This function encodes DNA sequences into a data frame of features.
#' @param dna_strings A character vector of DNA sequences.
#' @return A data frame containing the encoded features.
#' @export
#' @import randomForest
dna_encoding <- function(dna_strings) {
  nn <- nchar(dna_strings[1])
  seq_list <- strsplit(dna_strings, "")
  seq_df <- as.data.frame(do.call(rbind, seq_list))
  colnames(seq_df) <- paste0("V", 1:nn)
  return(seq_df)
}

#' Multiple Samples Prediction for m6A Sites
#'
#' This function takes a data frame of features and predicts m6A probability.
#' @param ml_fit A trained machine learning model.
#' @param feature_df A data frame containing features for multiple samples.
#' @param positive_threshold A numeric value for classification threshold.
#' @return A data frame with predicted probabilities and status.
#' @examples
#' \dontrun{
#' model <- readRDS(system.file("extdata", "rf_fit.rds", package="m6APrediction"))
#' df <- read.csv(system.file("extdata", "m6A_input_example.csv", package="m6APrediction"))
#' prediction_multiple(model, df, positive_threshold = 0.6)
#' }
#' @export
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5) {
  pred_prob <- predict(ml_fit, feature_df, type = "prob")[, "Positive"]
  pred_status <- ifelse(pred_prob >= positive_threshold, "Positive", "Negative")

  feature_df$predicted_m6A_prob <- pred_prob
  feature_df$predicted_m6A_status <- pred_status
  return(feature_df)
}

#' Single Sample Prediction for m6A Sites
#'
#' This function takes individual feature values and uses a pre-trained model.
#' @param ml_fit A trained machine learning model.
#' @param gc_content A numeric value for the GC content.
#' @param RNA_type A character string for RNA type.
#' @param RNA_region A character string for RNA region.
#' @param exon_length A numeric value for exon length.
#' @param distance_to_junction A numeric value for distance to junction.
#' @param evolutionary_conservation A numeric value for conservation.
#' @param DNA_5mer A character string representing the 5-mer.
#' @param positive_threshold A numeric value between 0 and 1.
#' @return A named vector with predicted probability and status.
#' @examples
#' \dontrun{
#' model <- readRDS(system.file("extdata", "rf_fit.rds", package="m6APrediction"))
#' prediction_single(model, gc_content = 0.6, RNA_type = "mRNA", RNA_region = "CDS",
#' exon_length = 12, distance_to_junction = 5, evolutionary_conservation = 0.8,
#' DNA_5mer = "GGACA", positive_threshold = 0.5)
#' }
#' @export
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length,
                              distance_to_junction, evolutionary_conservation, DNA_5mer,
                              positive_threshold = 0.5) {
  single_df <- data.frame(
    gc_content = gc_content,
    RNA_type = RNA_type,
    RNA_region = RNA_region,
    exon_length = exon_length,
    distance_to_junction = distance_to_junction,
    evolutionary_conservation = evolutionary_conservation,
    DNA_5mer = DNA_5mer,
    stringsAsFactors = FALSE
  )

  res_df <- prediction_multiple(ml_fit, single_df, positive_threshold)

  out_vec <- c(
    predicted_m6A_prob = res_df$predicted_m6A_prob[1],
    predicted_m6A_status = res_df$predicted_m6A_status[1]
  )
  return(out_vec)
}

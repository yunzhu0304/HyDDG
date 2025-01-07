#' Calculate Fold Changes for Experimental Samples
#'
#' @description Calculating fold changes for experimental samples by dividing
#'     each experimental sample's values by the mean of the control group values.
#'
#' @param data A normalized data frame where rows represent probes/genes and columns represent samples.
#' Each element of the data frame is expected to be on a log2 scale.
#' @param group.list A character factor vector indicating group information for each column of the `data`
#' data frame. The factor should contain exactly two levels: the first level corresponds to the control
#' group, and the second level corresponds to the experimental group. The order must match the column
#' names of `data`.
#'
#' @return A data frame of the same number of rows as `data` and columns corresponding to the experimental
#' group samples. The values represent fold changes (experimental group values divided by the control
#' group mean) for each sample. Row names correspond to the row names of `data`, and column names
#' correspond to the experimental group sample names.
#'
#' @details The function performs the following steps:
#' 1. Converts `data` from log2 scale to the original scale using a base-2 exponential transformation.
#' 2. Computes the mean of the control group samples for each row (probe/gene).
#' 3. Calculates fold changes by dividing each experimental sample's values by the control group mean.
#' 4. Outputs a data frame with the calculated fold changes for the experimental group.
#'
#' @examples
#' # Install and load the CLL package
#' # BiocManager::install("CLL")
#' library(CLL)
#' data("CLLbatch")
#'
#' # Normalize data using RMA
#' CLLrma <- rma(CLLbatch)
#'
#' # Convert normalized data to a data frame
#' ourData <- as.data.frame(exprs(CLLrma))
#'
#' # Define group information: first 12 columns as control, last 12 as treatment
#' group <- factor(c(rep("Control", 12), rep("Treat", 12)),levels = c("Control", "Treat"))
#'
#' # Apply the hyfit function
#' fit <- hyfit(data = ourData, group.list = group)
#'
#' @export
hyfit <- function(data,group.list){

  data <- 2^data

  control_mean <- rowMeans(data[, group.list == levels(group.list)[1], drop = FALSE])


  treatment_data <- data[, group.list == levels(group.list)[2], drop = FALSE]
  result <- sweep(treatment_data, 1, control_mean, FUN = "/")


  rownames(result) <- rownames(data)
  colnames(result) <- colnames(data)[group.list == levels(group.list)[2]]

  return(result)
}

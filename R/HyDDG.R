#' Identify Differentially Expressed Genes Using Hypergeometric Distribution
#'
#' Identifing differentially expressed genes (DEGs) by calculating
#' up- and down-regulated genes using a Hypergeometric Distribution-based method.
#' It computes the significance of each gene based on fold change and outputs a
#' data frame with statistical results.
#'
#' @param data A normalized data frame where rows represent probes/genes and columns represent samples.
#' The data should be in log2 scale.
#' @param group.list A character factor vector indicating group information for each column of the `data`
#' data frame. The factor should contain exactly two levels: the first level corresponds to the control
#' group, and the second level corresponds to the experimental group. The order must match the column
#' names of `data`.
#' @param adj.p.method A character string specifying the method for p-value adjustment.
#' Default is "BH" (Benjamini-Hochberg). Other options include: "holm", "hochberg",
#' "hommel", "bonferroni", "BH", "BY", "fdr", "none". See \code{\link[stats]{p.adjust}} for details.
#' @param BV A numeric value representing the boundary value for fold change determination.
#' Default is 1. Genes with fold change greater or smaller than this value will
#' be considered for up- or down-regulation, respectively.
#' @param ... Additional arguments passed to other internal functions if needed.
#'
#' @return A data frame named `HyDDGResult` with the following columns:
#' \itemize{
#'   \item \code{ID}: Gene or probe IDs corresponding to row names of the input `data`.
#'   \item \code{ave.expr}: Average expression value of the gene across all samples.
#'   \item \code{lg.p}: Log-transformed p-values, with a positive sign for up-regulation
#'   and a negative sign for down-regulation.
#'   \item \code{p.value}: The minimum p-value from hypergeometric and binomial tests.
#'   \item \code{p.adj}: Adjusted p-values based on the specified method in \code{adj.p.method}.
#'   \item \code{FC}: Average fold change for each gene, calculated as the mean of fold changes
#'   in experimental samples relative to controls.
#'   \item \code{logFC}: Log2-transformed fold change values.
#' }
#'
#' @details The function performs the following:
#' 1. Calls \code{hyfit} to calculate fold changes for experimental samples.
#' 2. For each gene (row), computes the hypergeometric p-values and binomial p-values
#'    for up- and down-regulated genes based on the boundary value (\code{BV}).
#' 3. Identifies the direction of regulation (up or down) and calculates the corresponding
#'    log-transformed p-value (\code{lg.p}).
#' 4. Adjusts p-values using the method specified in \code{adj.p.method}.
#' 5. Recommends genes with \code{abs(lg.p) ≥ 3} and \code{FDR < 0.05} as significant.
#'
#' @examples
#' # Load example data
#' # BiocManager::install("CLL")
#' library(CLL)
#' data("CLLbatch")
#'
#' # Normalize data using RMA
#' CLLrma <- rma(CLLbatch)
#'
#' # Convert normalized data to a data frame
#' ourData <- as.data.frame(exprs(CLLrma))[c(1:100),]
#'
#' # Define group information: first 12 columns as control, last 12 as treatment
#' group <- factor(c(rep("Control", 12), rep("Treat", 12)),
#'                 levels = c("Control", "Treat"))
#'
#' # Run the HyDDG function
#' HyDDGResult <- HyDDG(data = ourData, group.list = group)
#'
#' # Access the result
#' head(HyDDGResult)
#'
#' @export
HyDDG <- function(data, group.list, adj.p.method = "BH", BV = 1, ...){
  hyfitR <- hyfit(data,group.list)

  pb <- txtProgressBar(min = 0, max = nrow(hyfitR), style = 3)

  m_i <- ncol(hyfitR)
  N <- nrow(hyfitR) * ncol(hyfitR)


  HyDDGResult <- data.frame()

  for (j in 1:nrow(hyfitR)) {

    current_row <- hyfitR[j, ]

    up_values <- current_row[current_row > BV]
    down_values <- current_row[current_row < BV]

    k_up <- length(up_values)
    if (k_up > 0) {
      sorted_up <- sort(up_values, decreasing = TRUE)
      P_up_HGvalues <- sapply(1:k_up, function(i) {
        M_l <- sum(hyfitR >= sorted_up[i])
        sum(dhyper(i:m_i, M_l, N - M_l, m_i))
      })
      P_up_HG <- min(P_up_HGvalues)
      P_up_binom <- pbeta(0.5, k_up, m_i - k_up + 1)
    } else {
      P_up_HG <- Inf
      P_up_binom <- Inf
    }

    k_down <- length(down_values)
    if (k_down > 0) {
      sorted_down <- sort(down_values, decreasing = FALSE)
      P_down_HGvalues <- sapply(1:k_down, function(i) {
        M_l <- sum(hyfitR <= sorted_down[i])
        sum(dhyper(i:m_i, M_l, N - M_l, m_i))
      })
      P_down_HG <- min(P_down_HGvalues)
      P_down_binom <- pbeta(0.5, k_down, m_i - k_down + 1)
    } else {
      P_down_HG <- Inf
      P_down_binom <- Inf
    }

        P_min <- min(P_up_binom, P_down_binom, P_up_HG, P_down_HG)

    P_values <- c(P_up_binom, P_down_binom, P_up_HG, P_down_HG)

    if (which.min(P_values) %in% c(1, 3)) {
      lg_p <- -log10(P_min)
    } else {
      lg_p <- -log10(P_min) * -1
    }

    HyDDGResult <- rbind(HyDDGResult, data.frame(
      ID = rownames(hyfitR)[j],
      ave.expr = rowMeans(data)[j],
      lg.p = lg_p,
      p.value = P_min
    ))

     setTxtProgressBar(pb, j)
  }
  HyDDGResult$p.adj <- p.adjust(HyDDGResult$p.value, method = adj.p.method)
  HyDDGResult$FC <- rowMeans(hyfitR)
  HyDDGResult$logFC <- log2(HyDDGResult$FC)
  return(HyDDGResult)
  close(pb)

}

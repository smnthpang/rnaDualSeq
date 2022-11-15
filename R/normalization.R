# normalization.R

# normalization

#'Normalize data by using trimmed mean of M-values to convert size of data to an effective library size
#'
#' @param data dRNA-seq tab-delimited input csv file with read counts
#' @param pheno csv file containing sample names and time/condition of each sample
#'
#' @return data.norm DGEList
#'
norm_TMM <- function(data, pheno) {
  time <- pheno[colnames(data), "temporal"]
  dge <-
    edgeR::DGEList(counts = data,
                   group = time,
                   genes = rownames(data))
  data.norm <-
    edgeR::cpm(edgeR::calcNormFactors(dge, method = "TMM"), log = TRUE)
  return (data.norm)
}

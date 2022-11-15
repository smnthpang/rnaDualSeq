# normalization.R

# normalization

#'Normalize data by using TMM log to convert size of data to an effective library size
#'
#' @param data
#' @param pheno
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

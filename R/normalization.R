# normalization.R

# normalization

#'Normalize data by using trimmed mean of M-values to convert size of data to an effective library size
#'
#' @param data dRNA-seq tab-delimited input data frame with read counts
#' @param pheno data frame containing sample names and time/condition of each sample
#'
#' @return data.norm DGEList
#'
#' @examples
#' \dontrun{
#'  norm_TMM(host_example, phenotype_example)
#' }
#' @references
#'
#' Robinson MD, McCarthy DJ and Smyth GK (2010). edgeR: a Bioconductor package
#' for differential expression analysis of digital gene expression data.
#' Bioinformatics 26, 139-140
#'
#' @import edgeR
#' @export
#'
norm_TMM <- function(data, pheno) {

  time <- pheno[colnames(data), "groups"]
  dge <-
    edgeR::DGEList(counts = data,
                   group = time,
                   genes = rownames(data))
  data.norm <-
    edgeR::cpm(edgeR::calcNormFactors(dge, method = "TMM"), log = TRUE)
  return (data.norm)
}




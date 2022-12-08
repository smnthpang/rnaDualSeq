# differentialExpression.R


#' Separate data by categorizing them by hour groups (ie. WT.02_h, WT.04_h, etc.) to analyze
#' differential expression at different time periods
#'
#' The limma package is used to find differentially expressed (DE) genes, in which expression data
#' is described by linear modelling
#'
#' @param data dRNA-seq tab-delimited input csv file with read counts
#' @param pheno csv file containing sample names and time/condition of each sample
#'
#' @return list of differentially expressed genes compared with the baseline (0 hours)
#'
#' @examples
#' \dontrun{
#' norm <- norm_TMM(host_example, phenotype_example)
#' DE <- identifyDE(norm, phenotype_example)
#' }
#'
#' @references
#'
#' Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK (2015).
#' “limma powers differential expression analyses for RNA-sequencing and
#' microarray studies.” Nucleic Acids Research, 43(7), e47. doi: 10.1093/nar/gkv007
#'
#' Dinarvand, M., Kock, F., Al Mouiee, D., Vuong, K., Vijayan, A., Tanzim, A.
#' F., Azad, A. K. M., Penesyan, A., Castaño-Rodríguez, N., & Vafaee, F. (2022).
#' DSEQSB: A Systems Biology Approach to decipher dynamics of host-pathogen
#' interactions using temporal dual RNA-Seq Data. https://doi.org/10.1101/2022.02.28.482417
#'
#' Macho Rendón, J., Lang, B., Ramos Llorens, M., Tartaglia, G.G., and
#' Torrent Burgas, M. (2021). DualSeqDB: a database to assess the relevance of
#' bacterial genes during host infection. Nucleic Acids Res. 49, D687–D693.
#'
#' @export
#' @import limma
#' @import stats
#'

identifyDE <- function(data, pheno){

  time <- as.factor(pheno[colnames(data), "groups"])
  design <- model.matrix(~ 0 + time)

  colnames(design) <- levels(time)
  rownames(design) <- colnames(data)

  fit <- limma::lmFit(data, design = design)

  # Construct the contrast matrix using time slots as set of parameters from time to baseline
  cont.matrix <- limma::makeContrasts(WT.02_h - WT.00_h,
                                      WT.04_h - WT.00_h,
                                      WT.08_h - WT.00_h,
                                      WT.16_h - WT.00_h,
                                      WT.24_h - WT.00_h,
                                      levels = design)

  fit.cont <- limma::eBayes(limma::contrasts.fit(fit, cont.matrix))
  summa.fit <- limma::decideTests(fit.cont)

  # statistically significant DE genes
  # false discovery rate (FDR) cut-off ≤ 0.05
  # log fold change (LFC) of at least 2-fold

  significant <- function(fit.cont, coefficient = 1, lfc = 1, adjP = 0.05){
    topgenes <- limma::topTable(fit.cont, n = Inf, coef = coefficient)
    upderegulation = with(topgenes, logFC > lfc & adj.P.Val < adjP)
    downderegulation = with(topgenes, logFC < lfc & adj.P.Val < adjP)
    return (topgenes)
  }

  differentialExpression <- list()

  # for each different time period
  for (n in 1:5){
    differentialExpression[[n]] <- significant(fit.cont, coefficient = n)
  }

  return (differentialExpression)
}

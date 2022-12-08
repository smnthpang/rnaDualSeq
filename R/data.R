#' dRNA-seq read count matrix from the host with 6 time points
#' (0 h, 2 h, 4 h, 8 h, 16 h, and 24 h)
#'
#' @source NCBI’s Gene Expression Omnibus (accession number GSE60144)
#'
#' @format a data frame with 40375 rows and 18 variables
#' \describe{
#' \item{HeLa.S3.WT.00_h.R1}{round 1 hour 0 read counts}
#' \item{HeLa.S3.WT.00_h.R2}{round 2 hour 0 read counts}
#' \item{HeLa.S3.WT.00_h.R3}{round 3 hour 0 read counts}
#' \item{HeLa.S3.WT.02_h.R1}{round 1 hour 2 read counts}
#' \item{HeLa.S3.WT.02_h.R2}{round 2 hour 2 read counts}
#' \item{HeLa.S3.WT.02_h.R3}{round 3 hour 2 read counts}
#' \item{HeLa.S3.WT.04_h.R1}{round 1 hour 4 read counts}
#' \item{HeLa.S3.WT.04_h.R2}{round 2 hour 4 read counts}
#' \item{HeLa.S3.WT.04_h.R3}{round 3 hour 4 read counts}
#' \item{HeLa.S3.WT.08_h.R1}{round 1 hour 8 read counts}
#' \item{HeLa.S3.WT.08_h.R2}{round 2 hour 8 read counts}
#' \item{HeLa.S3.WT.08_h.R3}{round 3 hour 8 read counts}
#' \item{HeLa.S3.WT.16_h.R1}{round 1 hour 16 read counts}
#' \item{HeLa.S3.WT.16_h.R2}{round 2 hour 16 read counts}
#' \item{HeLa.S3.WT.16_h.R3}{round 3 hour 16 read counts}
#' \item{HeLa.S3.WT.24_h.R1}{round 1 hour 24 read counts}
#' \item{HeLa.S3.WT.24_h.R2}{round 2 hour 24 read counts}
#' \item{HeLa.S3.WT.24_h.R3}{round 3 hour 24 read counts}
#' }
#'
#' @examples
#' \dontrun{
#'
#' }
#'
"host_example"

#' dRNA-seq read count matrix from the bacterium with 6 time points
#' (0 h, 2 h, 4 h, 8 h, 16 h, and 24 h)
#'
#' @source NCBI’s Gene Expression Omnibus (accession number GSE60144)
#'
#' @format a data frame with 1153 rows and 18 variables
#' \describe{
#' \item{HeLa.S3.WT.00_h.R1}{round 1 hour 0 read counts}
#' \item{HeLa.S3.WT.00_h.R2}{round 2 hour 0 read counts}
#' \item{HeLa.S3.WT.00_h.R3}{round 3 hour 0 read counts}
#' \item{HeLa.S3.WT.02_h.R1}{round 1 hour 2 read counts}
#' \item{HeLa.S3.WT.02_h.R2}{round 2 hour 2 read counts}
#' \item{HeLa.S3.WT.02_h.R3}{round 3 hour 2 read counts}
#' \item{HeLa.S3.WT.04_h.R1}{round 1 hour 4 read counts}
#' \item{HeLa.S3.WT.04_h.R2}{round 2 hour 4 read counts}
#' \item{HeLa.S3.WT.04_h.R3}{round 3 hour 4 read counts}
#' \item{HeLa.S3.WT.08_h.R1}{round 1 hour 8 read counts}
#' \item{HeLa.S3.WT.08_h.R2}{round 2 hour 8 read counts}
#' \item{HeLa.S3.WT.08_h.R3}{round 3 hour 8 read counts}
#' \item{HeLa.S3.WT.16_h.R1}{round 1 hour 16 read counts}
#' \item{HeLa.S3.WT.16_h.R2}{round 2 hour 16 read counts}
#' \item{HeLa.S3.WT.16_h.R3}{round 3 hour 16 read counts}
#' \item{HeLa.S3.WT.24_h.R1}{round 1 hour 24 read counts}
#' \item{HeLa.S3.WT.24_h.R2}{round 2 hour 24 read counts}
#' \item{HeLa.S3.WT.24_h.R3}{round 3 hour 24 read counts}
#' }
#'
#' @examples
#' \dontrun{
#'
#' }
#'
"pathogen_example"

#' phenotype data frame that contains names of the columns headers corresponding
#' expression (read count) data frame, and group name which is time/condition

#'
#' @source NCBI’s Gene Expression Omnibus (accession number GSE60144)
#'
#' @format a data frame with 18 rows and 1 variables
#' \describe{
#' \item{groups}{time or condition}
#' }
#'
#' @examples
#' \dontrun{
#'
#' }
#'
"phenotype_example"

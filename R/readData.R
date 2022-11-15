# readData.R

# read_csv.R

#' Read csv files to access data
#'
#' @param path file path
#'
#' @return data frame (data.frame) containing a representation of the data in the file
#'
#' @import utils
#' @export
#'

read_csv <- function(path) {
  return(read.csv(path, row.names = 1))
}

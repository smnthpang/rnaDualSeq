#' Launch Shiny App for rnaDualSeq
#'
#' A function that launches the Shiny app for rnaDualSeq
#'
#' @return No return value but open up a Shiny page.
#'
#' @examples
#' \dontrun{
#'
#' rnaDualSeq::runRnaDualSeq()
#' }
#'
#' @references
#' Grolemund, G. (2015). Learn Shiny - Video Tutorials.
#'     \href{https://shiny.rstudio.com/tutorial/}{Link}
#'
#' @export
#' @importFrom shiny runApp

runRnaDualSeq <- function() {
  appDir <- system.file("shiny-scripts",
                        package = "rnaDualSeq")
  shiny::runApp(appDir, display.mode = "normal")
  return()
}
# [END]

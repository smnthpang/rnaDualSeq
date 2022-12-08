# volcanoPlot.R

#' Creates 5 volcano plots denoting each time period in data set visualing
#' significant upregulation and downregulation genes
#'
#' @param name what gene is being graphed
#' @param de_norm A nested list containing differentially expressed results
#' @param timeperiods a vector containing labels for each time period
#'
#' @return graphs of volcano plots per time period and saves pngs images
#' into a results folder
#'
#' @examples
#' \dontrun{
#' norm <- norm_TMM(host_example, phenotype_example)
#' DE <- identifyDE(norm, phenotype_example)
#' volcanoPlot("Host Example", DE)
#' }
#'
#' @references
#'
#' Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis.
#' Springer-Verlag New York. ISBN 978-3-319-24277-4, https://ggplot2.tidyverse.org.
#'
#'
#'
#' @import ggplot2
#' @export
#'
volcanoPlot <- function(name, de_norm, timeperiods = c("2h", "4h", "8h", "16h", "24h")) {
  myplots <- list()
  for (n in 1:length(timeperiods)){
    de_norm[[n]]$diffexpressed <- "NO"
    de_norm[[n]]$diffexpressed[de_norm[[n]]$logFC > .58 & de_norm[[n]]$adj.P.Val < 0.01] <- "UP"
    de_norm[[n]]$diffexpressed[de_norm[[n]]$logFC < .58 & de_norm[[n]]$adj.P.Val < 0.01] <- "DOWN"

    de_norm[[n]]$de_normlabel <- NA

    graph <- ggplot2::ggplot(data=de_norm[[n]],ggplot2::aes(x=logFC,y=-log10(P.Value),col=diffexpressed, label=de_normlabel))+
      ggplot2::ggtitle(paste(name,"_",timeperiods[n],sep="")) +
      ggplot2::geom_point()+
      ggplot2::theme_minimal()+
      #geom_text_repel()+
      ggplot2::scale_color_manual(values=c('blue', 'black', 'red'))+
      ggplot2::theme(text=ggplot2::element_text(size=20))

    # save plot to list of plots
    myplots[[n]] <- graph

    check_if_directory_exists <- function(dir_path){
      if(!dir.exists(dir_path)){
        dir.create(dir_path, recursive = TRUE)
      }
    }
    check_if_directory_exists("Results/VolcanoPlots")
    ggplot2::ggsave(paste(name,"_",timeperiods[n],".png",sep=""), graph, path = "Results/VolcanoPlots")

  }
  return(myplots)
}

# volcanoPlot.R

#' Creates 5 volcano plots denoting each time period in data set visualing
#' significant upregulation and downregulation genes
#'
#' @param name what gene is being graphed
#' @param de_norm A nested list containing differentially expressed results
#' @param timeperiods a vector containing labels for each time period
#'
#' @return images of volcano plots per time period
#'
#' @import ggplot2
#' @export
#'
volcanoPlot <- function(name, de_norm, timeperiods = c("2h", "4h", "8h", "16h", "24h")) {
  for (n in 1:5){
    de_norm[[n]]$diffexpressed <- "NO"
    de_norm[[n]]$diffexpressed[de_norm[[n]]$logFC > .58 & de_norm[[n]]$adj.P.Val < 0.01] <- "UP"
    de_norm[[n]]$diffexpressed[de_norm[[n]]$logFC < .58 & de_norm[[n]]$adj.P.Val < 0.01] <- "DOWN"

    de_norm[[n]]$de_normlabel <- NA

    graph <- ggplot(data=de_norm[[n]],aes(x=logFC,y=-log10(P.Value),col=diffexpressed, label=de_normlabel))+
      geom_point()+
      theme_minimal()+
      #geom_text_repel()+
      scale_color_manual(values=c('blue', 'black', 'red'))+
      theme(text=element_text(size=20))

    check_if_directory_exists <- function(dir_path){
      if(!dir.exists(dir_path)){
        dir.create(dir_path, recursive = TRUE)
      }
    }
    check_if_directory_exists("Results/VolcanoPlots")
    ggsave(paste(name,"_",timeperiods[n],".pdf",sep=""), graph, path = "Results/VolcanoPlots")

  }
}

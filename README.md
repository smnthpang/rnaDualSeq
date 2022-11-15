rnaDualSeq
================

## Description

rnaDualSeq is an R package to house functions that help facilitate and
visualize differentially expressed genes within host-pathogen infection
studies.

## Installation

To run latest version of package:

``` bash
require("devtools")
devtools::install_github("smnthpang/rnaDualSeq", build_vignettes = TRUE)
library("rnaDualSeq")
```

To run shinyApp: Under construction

## Overview

``` python
ls("package:rnaDualSeq")
data(package = "rnaDualSeq")
```

rnaDualSeq has functions that

-   read in an input data file containing the read counts of a genes
    where samples are matched to different time periods in sections of
    hours
-   read in a phenotype file that describes the sample to groups
-   normalises the raw counts within the data file by using the trimmed
    mean of M-values (TMM) to correct for RNA composition differences.
    After that, the normalised data is log2 transformed. (see function
    norm_TMM)
-   identifies differentially expressed genes using pathogen and host
    genome as inputs and sorts them by time period (see function
    identifyDE)
-   visualizes the differentially expressed genes by plotting each time
    period on a volcano plot where upregulated/downregulated genes are
    labeled (see function volcanoPlot)

See vignettes for tutorial on package:

``` bash
browseVignettes("<rnaDualSeq>")
```

An overview of the package is illustrated below:

## Contributing

The author of this package is Samantha Pang. The

The obtainGTF(), obtainDNA() and obtainCDNA() functions makes use of
rvest package to web scrape from online database. The txdbObj function
makes use of makeTxDbFromGFF() function from GenomicFeatures package to
create a txdb object from GTF file. The generateMatrix() function uses
pheatmap function from pheatmap R package to plot the heatmap. The stats
and ggplot2 R packages are used for principle component analysis in
plotPCA() function. The rstudioapi CRAN package is used to send commands
to terminal for downloading Salmon software in installSalmon() function,
as well as running Salmon commands for indexing and quantification in
quantification() function.

## References

Dinarvand, M., Kock, F., Al Mouiee, D., Vuong, K., Vijayan, A., Tanzim,
A. F., Azad, A. K. M., Penesyan, A., Castaño-Rodríguez, N., & Vafaee, F.
(2022). DSEQSB: A Systems Biology Approach to decipher dynamics of
host-pathogen interactions using temporal dual RNA-Seq Data.
<https://doi.org/10.1101/2022.02.28.482417>

Doyle, M. (2022, October 18). Galaxy training: Visualization of RNA-seq
results with volcano plot in R. Galaxy Training Network. Retrieved
November 15, 2022, from
<https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/rna-seq-viz-with-volcanoplot-r/tutorial.html>

Love, M. I., Anders, S., Kim, V., & Huber, W. (n.d.). RNA-seq workflow:
gene-level exploratory analysis and differential expression. RNA-seq
workflow: Gene-level exploratory analysis and differential expression.
Retrieved November 15, 2022, from
<https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html>

Walch, P., Selkrig, J., Knodler, L. A., Rettel, M., Stein, F.,
Fernandez, K., Viéitez, C., Potel, C. M., Scholzen, K., Geyer, M.,
Rottner, K., Steele-Mortimer, O., Savitski, M. M., Holden, D. W., &
Typas, A. (2021). Global mapping of salmonella enterica-host
protein-protein interactions during infection. Cell Host & Microbe,
29(8). <https://doi.org/10.1016/j.chom.2021.06.004>

Westermann, A. J., Förstner, K. U., Amman, F., Barquist, L., Chao, Y.,
Schulte, L. N., Müller, L., Reinhardt, R., Stadler, P. F., & Vogel, J.
(2016). Dual RNA-seq unveils noncoding RNA functions in host–pathogen
interactions. Nature, 529(7587), 496–501.
<https://doi.org/10.1038/nature16547>

Macho Rendón, J., Lang, B., Ramos Llorens, M., Tartaglia, G.G., and
Torrent Burgas, M. (2021). DualSeqDB: a database to assess the relevance
of bacterial genes during host infection. Nucleic Acids Res. 49,
D687–D693.

Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK (2015).
“limma powers differential expression analyses for RNA-sequencing and
microarray studies.” Nucleic Acids Research, 43(7), e47. doi:
10.1093/nar/gkv007. Wickham H (2016). ggplot2: Elegant Graphics for Data
Analysis. Springer-Verlag New York. ISBN 978-3-319-24277-4,
<https://ggplot2.tidyverse.org>.

Robinson MD, McCarthy DJ and Smyth GK (2010). edgeR: a Bioconductor
package for differential expression analysis of digital gene expression
data. Bioinformatics 26, 139-140

Robinson MD, McCarthy DJ and Smyth GK (2010). edgeR: a Bioconductor
package for differential expression analysis of digital gene expression
data. Bioinformatics 26, 139-140

## Acknowledgements

This package was developed as part of an assessment for 2022 BCB410H:
Applied Bioinformatics course at the University of Toronto, Toronto,
CANADA. <PackageName>welcomes issues, enhancement requests, and other
contributions. To submit an issue, use the GitHub issues

## Including Code

You can include R code in the document as follows:

``` r
summary(cars)
```

    ##      speed           dist       
    ##  Min.   : 4.0   Min.   :  2.00  
    ##  1st Qu.:12.0   1st Qu.: 26.00  
    ##  Median :15.0   Median : 36.00  
    ##  Mean   :15.4   Mean   : 42.98  
    ##  3rd Qu.:19.0   3rd Qu.: 56.00  
    ##  Max.   :25.0   Max.   :120.00

## Including Plots

You can also embed plots, for example:

![](README_files/figure-gfm/pressure-1.png)<!-- -->

Note that the `echo = FALSE` parameter was added to the code chunk to
prevent printing of the R code that generated the plot.

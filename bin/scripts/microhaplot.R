

packages <- c("dplyr", "openssl", "httr", "devtools", "feather", "shiny", "shinyBS", "shinyWidgets", "httr", "plotly", "DT", "ggplot2")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())),repos='http://cran.rstudio.com/')  
}


devtools::install_github("ngthomas/microhaplot", build_vignettes = TRUE)

microhaplot::mvHaplotype("~/results/Final_results/Microhaplot/Shiny")

library(microhaplot)

# to access package sample case study dataset of rockfish
run.label <- "snakemake_haps_microhaplot"
sam.path <- "~/results/Final_results/Microhaplot/"
label.path <- "~/results/Final_results/Microhaplot/label.txt"
vcf.path <- "~/results/Final_results/Microhaplot/snps.vcf"
app.path <- "~/results/Final_results/Microhaplot/Shiny/microhaplot"

# for your dataset: customize the following paths
# sam.path <- "~/microhaplot/extdata/"
# label.path <- "~/microhaplot/extdata/label.txt"
# vcf.path <- "~/microhaplot/extdata/sebastes.vcf"
# app.path <- "~/Shiny/microhaplot"

haplo.read.tbl <- runHaplot(run.label = run.label,
           sam.path=sam.path,
           label.path=label.path,
           vcf.path=vcf.path,
           app.path=app.path)


           

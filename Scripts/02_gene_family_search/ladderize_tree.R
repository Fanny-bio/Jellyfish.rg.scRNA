##Rscript .R -i input.noquote.newick.tre -o output_prefix 

suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidytree))
suppressPackageStartupMessages(library(phylotools))
suppressPackageStartupMessages(library(phangorn))
suppressPackageStartupMessages(library(treeio))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(flextable))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggplotify))
suppressPackageStartupMessages(library(eoffice))
suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("-i","--input"), type="character", default=NULL,
              help="inputTree file' [default %default]",
              dest="inputTree"),
  make_option(c("-o","--output"), type="character", default="output_ggtree",
              help="output file prefix' [default %default]",
              dest="outputPrefix"),
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print out all parameter settings [default]")
)  

options(error=traceback)

parser <- OptionParser(usage = "%prog -i input.noquote.tre -o out.prefix [options]",option_list=option_list)
opt = parse_args(parser)


mytree <- read.tree(opt$inputTree)


  out_tree <- ladderize(mytree)
  write.tree(out_tree,file=paste0(opt$outputPrefix,".ladderized.tree.nwk"))


  
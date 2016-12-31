#!/usr/bin/env Rscript

## Copyright (c) 2016 by Joseph Tran. This work is made available under the terms of the Creative Commons Attribution-ShareAlike 3.0 license, http://creativecommons.org/licenses/by-sa/3.0/.

################################################
##
## @script: plot_heatmap_rnaseq.r
##
## @author: Joseph Tran (Joseph.Tran@ips2.universite-paris-saclay.fr)
##
## @version: 0.0.1.0
##
## @date: 2016-09-05
##
## @description: This R script helps plotting a heatmap on differentially expressed genes list (SARTools output). 
##
###############################################
rm(list=ls())                                        # remove all the objects from the R session

version <- "0.0.1.0"

copyright <- "Copyright (c) 2016 by Joseph Tran. This work is made available under the terms of the Creative Commons Attribution-ShareAlike 3.0 license, http://creativecommons.org/licenses/by-sa/3.0/."

## execution time: start
T1<-Sys.time()

###
### FUNCTIONS ###
###

###
### check for installed package
###
is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1])

###
### Capturing warnings/errors
###
tryCatch.W.E <- function(expr)
{
  W <- NULL
  w.handler <- function(w)
  { # warning handler
    W <<- w
    invokeRestart("muffleWarning")
  }
  list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                   warning = w.handler),
       warning = W)
}

##
## Options and usage
##
MIRROR <- "http://mirrors.ebi.ac.uk/CRAN/"

### optparse
if (is.element('optparse', installed.packages()[,1]))
{
  suppressPackageStartupMessages(require('optparse'));
} else
{
  install.packages('optparse', repos=MIRROR);
  suppressPackageStartupMessages(library('optparse'));
}

### pheatmap
if (is.element('pheatmap', installed.packages()[,1]))
{
  suppressPackageStartupMessages(require('pheatmap'));
} else
{
  install.packages('pheatmap', repos=MIRROR);
  suppressPackageStartupMessages(library('pheatmap'));
}

### log4r
if (is.element('log4r', installed.packages()[,1]))
{
  suppressPackageStartupMessages(require('log4r'));
} else
{
  install.packages('log4r', repos=MIRROR);
  suppressPackageStartupMessages(library('log4r'));
}

### tools
suppressPackageStartupMessages(library(tools))

### RColorBrewer
if (is.element('RColorBrewer', installed.packages()[,1]))
{
  suppressPackageStartupMessages(require('RColorBrewer'));
} else
{
  install.packages('RColorBrewer', repos=MIRROR);
  suppressPackageStartupMessages(library('RColorBrewer'));
}

## OPTIONS
option_list <- list(
  # general
  make_option(c('-v', '--verbosity-level'), type = "integer", dest = "verbosity", default = 1, help = "Verbosity threshold (5=DEBUG, 4=INFO 3=WARN, 2=ERROR, 1=FATAL)"),
  make_option(c("-i", "--input"), type="character", dest="input", help="The DE genes input file path (cf. SARTools output) [mandatory]", metavar="deg_infile"),
  make_option(c("-d", "--working-directory"), type="character", dest="workDir", default=getwd(), help="The working directory [default %default]", metavar="work_dir"),
  make_option(c("-M", "--main"), type="character", dest="main", default="group", help="The first factor (main effect) name [default %default]", metavar="main_effect"),
  make_option(c("-m","--main_levels"), type="character", dest="mainLevels", help="The factor1 levels separated by comma", metavar="main_levels"),
  make_option(c("-B","--batch"), type="character", dest="batch", default="patient", help="The second factor (batch effect) name [default %default]", metavar="batch_effect"),
  make_option(c("-b","--batch_levels"), type="character", dest="batchLevels", help="The factor2 levels separated by comma", metavar="batch_levels"),
  make_option(c("-c", "--genes_count"), type="double", dest="genesCount", default=-1, help="The number of top genes, with the highest variance/fold change, to plot [default %default meaning all genes]", metavar="genes_to_plot"),
  make_option(c("-s", "--sort-DEG"), type="character", dest="sortDEG", default="log2FC_padj", help="The DEG sort mode. 2 possible modes: i. log2FC_padj, sort DEG by log2foldchange (increasing if down-regulated, or decreasing if up-regulated), then adjusted pvalue (increasing); ii. padj, only by adjusted pvalue (increasing) [log2FC_padj, padj]", metavar="sort_deg_mode")
) 

parser <- OptionParser(usage = "usage:  %prog mandatory [options]", option_list=option_list)
opt <- parse_args(parser, positional_arguments = FALSE, print_help_and_exit=FALSE)

if (opt$help) {
  cat(paste("Plot a heatmap on differentially expressed genes list", sep=""))
  cat(paste("\nauthor:  Joseph Tran <Joseph.Tran@ips2.universite-paris-saclay.fr>\nversion: ", version, "\n"))
  print_help(parser)
  cat(paste("notes:"))
  cat(paste("\n", copyright, "\n", sep=""))
  quit()
}

##
## Load and install libraries
##

##
## Logger
##
logger <- create.logger(logfile = "", level = verbosity(opt$verbosity))

null_values <- c("NULL", "null")

###
### Check user input parameters
###

## General

if (file.access(opt$input, mode=4) == -1) {
  error(logger, paste("DE genes input file ", opt$input, " does not exist or is not readable.", sep=""))
  stop(sprintf("DE genes input file ( %s ) does not exist or is not readable.", opt$input))
} else
{
  debug(logger, paste("OK DE genes input file ", opt$input, " does exist and is readable", sep=""))
  degfn <- basename(opt$input)
  deg <- file_path_as_absolute(opt$input)
}

if (file.access(opt$workDir, mode=2) == -1) {
  error(logger, paste("Working directory ", opt$workDir, " does not exist or is not writeable.", sep=""))
  stop(sprintf("Given working directory ( %s ) does not exist or is not writeable.", opt$workDir))      
} else
{
  workDir <- opt$workDir
  setwd(workDir)
  debug(logger, paste("OK Working directory ", workDir, " exists.", sep=""))
}

if (is.null(opt$main)) {
  error(logger, paste("Main effect string cannot be null. It is mandatory to set a non null main effect string.", sep=""))
  stop(sprintf("Main effect string ( %s ) is null. It is mandatory to set a non null main effect string.", opt$main))
} else 
{
  mainEff <- opt$main
  debug(logger, paste("Main effect string is now set to ", mainEff, sep=""))
}

if (is.null(opt$batch)) {
  error(logger, paste("Batch effect string cannot be null. It is mandatory to set a non null batch effect string.", sep=""))
  stop(sprintf("Batch effect string ( %s ) is null. It is mandatory to set a non null batch effect string.", opt$batch))
} else 
{
  batchEff <- opt$batch
  debug(logger, paste("Batch effect string is now set to ", batchEff, sep=""))
}

if (is.null(opt$mainLevels)) {
  error(logger, paste("Main effect levels string cannot be null. It is mandatory to set a non null main effect levels string. Each level is separated by a comma.", sep=""))
  stop(sprintf("Main effect levels string ( %s ) is null. It is mandatory to set a non null main effect levels string. Each level is separated by a comma.", opt$mainLevels))
} else 
{
  mainEffLevels <- trimws(unlist(strsplit(opt$mainLevels,split = ",")))
  debug(logger, paste("Main effect levels list is now set to (", paste(mainEffLevels, collapse = ","), ")", sep=""))
}

if (is.null(opt$batchLevels)) {
  error(logger, paste("Batch effect levels string cannot be null. It is mandatory to set a non null batch effect levels string. Each level is separated by a comma.", sep=""))
  stop(sprintf("Batch effect levels string ( %s ) is null. It is mandatory to set a non null batch effect levels string. Each level is separated by a comma.", opt$main))
} else 
{
  batchEffLevels <- trimws(unlist(strsplit(opt$batchLevels,split = ",")))
  debug(logger, paste("Batch effect levels list is now set to (", paste(batchEffLevels, collapse = ","), ")", sep=""))
}

if (is.null(opt$genesCount)) {
  error(logger, paste("Genes count cannot be null. It is mandatory to set a non null genes count.", sep=""))
  stop(sprintf("Genes count ( %d ) is null. It is mandatory to set a non null genes count.", opt$genesCount))
} else 
{
  genesCount <- opt$genesCount
  debug(logger, paste("Genes count is now set to ", genesCount, sep=""))
}

if (is.null(opt$sortDEG)) {
  error(logger, paste("The DEG sorting mode string cannot be null. It is mandatory to set a non null type DEG sorting mode string.", sep=""))
  stop(sprintf("The DEG sorting mode  string ( %s ) is null. It is mandatory to set a non null DEG sorting mode string [log2FC_padj or padj].", opt$sortDEG))
} else 
{
  sortDEG <- opt$sortDEG
  debug(logger, paste("The DEG sorting mode string string is now set to ", sortDEG, sep=""))
}

# filter and plot DE genes 
## read dataset
info(logger, paste("Reading DE genes list file into data frame", sep=""))
df <- read.delim(deg)

## order genes: should be ordered by log2FoldChange then by adjusted pvalue
## if DEG list is a down expressed list (down pattern in filename), increasing ordering is applied on log2FoldChange 
## else decreasing ordering is applied
if (sortDEG == "log2FC_padj") {
  info(logger, paste0("Ordering DE genes list file on log2FoldChange and adjusted p-value colums"))
  if (grepl("down", degfn, ignore.case = T)) {
    debug(logger, paste0("Sorting down-regulated genes: log2FoldChange (increasing) and adjusted p-value (increasing)"))
    df <- df[with(df, order(log2FoldChange, padj)), ]
  } else {
    debug(logger, paste0("Sorting up-regulated genes: log2FoldChange (decreasing) and adjusted p-value (increasing)"))
    df <- df[with(df, order(-log2FoldChange, padj)), ]
  }
} else if (sortDEG == "padj") {
  # no log2FoldChange ordering
  info(logger, paste0("Ordering DE genes list file on adjusted p-value colum"))
  df <- df[with(df, order(padj)), ]
}

## extract normalized count: columns should start with "norm.*"
info(logger, paste("Selecting normalized counts into a new matrix", sep=""))
mat.norm <- df[grepl("^norm.*", colnames(df))]
info(logger, paste("Applying log2 + pseudocount(1) transformation on normalized counts into a new matrix", sep=""))
mat.log2p1 <- log2(mat.norm + 1)
debug(logger, paste("Setting matrix row names with genes Id", sep=""))
rownames(mat.log2p1) <- df$Id

# subset
info(logger, paste("Subsetting DE genes, taking n=", genesCount, sep=""))
if (nrow(mat.log2p1) < genesCount || genesCount == -1) {
  genesCount <- nrow(mat.log2p1) 
}
mat.log2p1 <- mat.log2p1[1:genesCount,]

## Centering
info(logger, paste("Centering matrix counts by row", sep=""))
mat.log2p1 <- mat.log2p1 - rowMeans(mat.log2p1)

## annot
info(logger, paste("Setting columns annotation", sep=""))
annotation_col <- data.frame(patient=factor(rep(batchEffLevels, length(mainEffLevels))), group=factor(unlist(lapply(mainEffLevels, function(g){rep(g, length(batchEffLevels))}))))
colnames <- paste0(annotation_col$patient, annotation_col$group)
info(logger, paste("Setting row names to columns annotation", sep=""))
rownames(annotation_col) <- colnames

# update mat colnames 
info(logger, paste("Setting matrix column names to match columns annotation row names", sep=""))
info(logger, paste("before: ", paste(colnames(mat.log2p1), collapse = ","), sep=""))
colnames(mat.log2p1) <- colnames
info(logger, paste("after: ", paste(colnames(mat.log2p1), collapse = ","), sep=""))

## clustering
fout <- paste0(basename(file_path_sans_ext(deg)), "_heatmap_n", genesCount, "_sort", sortDEG, "_log2pc1", ".pdf")
info(logger, paste("Setting plot output file name: ", fout, sep=""))
info(logger, paste("Running the clustering", sep=""))
# row and column clustering
# pheatmap(mat.log2p1, cellheight = 10, filename = fout, annotation_col = annotation_col)
# no row clustering
pheatmap(mat.log2p1, cellheight = 10, filename = fout, annotation_col = annotation_col, cluster_rows = FALSE)
# change colors palette
pheatmap(mat.log2p1, cellheight = 10, filename = fout, annotation_col = annotation_col, cluster_rows = FALSE,
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(200))

unlink('Rplots.pdf')

# R session
info(logger, sessionInfo())

## execution time: end
T2<-Sys.time()
Tdiff= difftime(T2, T1)
info(logger, paste("Execution time : ", Tdiff,"seconds\n", sep=""))

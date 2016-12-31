#!/usr/bin/env Rscript

## Copyright (c) 2016 by Joseph Tran. This work is made available under the terms of the Creative Commons Attribution-ShareAlike 3.0 license, http://creativecommons.org/licenses/by-sa/3.0/.

################################################
##
## @script: filter_gtf_by_deg.r
##
## @author: Joseph Tran (Joseph.Tran@ips2.universite-paris-saclay.fr)
##
## @version: 0.0.1.0
##
## @date: 2016-09-05
##
## @description: This R script helps filtering gff/gtf with differentially expressed genes list. 
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

### rtracklayer
if (is.element('rtracklayer', installed.packages()[,1]))
{
    suppressPackageStartupMessages(require('rtracklayer'));
} else
{
    install.packages('rtracklayer', repos=MIRROR);
    suppressPackageStartupMessages(library('rtracklayer'));
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

## OPTIONS
option_list <- list(
    # general
    make_option(c('-v', '--verbosity-level'), type = "integer", dest = "verbosity", default = 1, help = "Verbosity threshold (5=DEBUG, 4=INFO 3=WARN, 2=ERROR, 1=FATAL)"),
    make_option(c("-g", "--gff"), type="character", dest="gff", help="The reference gff/gtf file path [mandatory]", metavar="gff_or_gtf"),
    make_option(c("-p", "--glob-pattern"), type="character", dest="globPattern", help="The glob pattern to find file(s) for [mandatory]", metavar="pattern"),
    make_option("--working-directory", type="character", dest="workDir", default=getwd(), help="The working directory containing the file(s) matching the pattern [default %default]", metavar="work_dir"),
    make_option("--recursive", type="logical", dest="recursive", default=FALSE, help="Search recursively for file matching the pattern from the working directory [default %default]"),
    make_option("--out-prefix", type="character", dest="outPrefix", help="Prefix to attach to output file(s)", metavar="out_prefix"),
    make_option("--exclude", type="character", dest="exclude", help="TODO: The glob pattern to exclude file(s) from gff/gtf filtering", metavar="exclude")
) 

parser <- OptionParser(usage = "usage:  %prog mandatory [options]", option_list=option_list)
opt <- parse_args(parser, positional_arguments = FALSE, print_help_and_exit=FALSE)

if (opt$help) {
    cat(paste("Filter GFF/GTF with differentially expressed genes list", sep=""))
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

if (file.access(opt$workDir, mode=2) == -1) {
    error(logger, paste("Working directory ", opt$workDir, " does not exist or is not writeable.", sep=""))
    stop(sprintf("Given working directory ( %s ) does not exist or is not writeable.", opt$workDir))      
} else
{
    workDir <- opt$workDir
    debug(logger, paste("OK Working directory ", workDir, " exists.", sep=""))
}

if (file.access(opt$gff, mode=4) == -1) {
    error(logger, paste("gff/gtf file ", opt$gff, " does not exist or is not readable.", sep=""))
    stop(sprintf("Given gff/gtf file ( %s ) does not exist or is not readable.", opt$gff))
} else
{
    debug(logger, paste("OK gff/gtf file ", opt$gff, " does exist and is readable", sep=""))
    fn_gff <- opt$gff
    gff_ext <- file_ext(fn_gff)
    debug(logger, paste("Reading gff/gtf file into R", sep=""))
    gff <- readGFF(fn_gff)
    debug(logger, paste("gff/gtf file imported into R.", sep=""))
}

if (is.null(opt$globPattern)) {
    error(logger, paste("Glob pattern cannot be null. It is mandatory to set a glob pattern to search for differentially expressed genes file(s) to filter gff/gtf.", sep=""))
    stop(sprintf("Glob pattern ( %s ) is null. It is mandatory to set a glob pattern to search for differentially expressed genes file(s) to filter gff/gtf.", opt$globPattern))
} else 
{
    globPattern <- opt$globPattern
    debug(logger, paste("Glob pattern is now set to ", globPattern, sep=""))
}

if (!is.null(opt$exclude)) {
    exclude <- opt$exclude
    debug(logger, paste("Exclude pattern is now set to ", exclude, sep=""))
} else 
{
    exclude <- NULL
    debug(logger, paste("Exclude pattern is now set to NULL", sep=""))
}

if (!is.null(opt$outPrefix)) {
    outPrefix <- opt$outPrefix
    debug(logger, paste("Output file(s) prefix is now set to ", outPrefix, sep=""))
} else
{
    outPrefix <- NULL
    debug(logger, paste("Output file(s) prefix is now set to NULL", sep=""))
}


if (!is.null(opt$recursive)) {
    recursive <- opt$recursive
    debug(logger, paste("Recursive search is now set to ", recursive, sep=""))
}


## list files in workDir using globPattern with or without recursive search
if (recursive) {
    fls <- list.files(workDir, pattern = globPattern, full.names = TRUE, recursive = recursive)
} else 
{
    fls <- list.files(workDir, pattern = globPattern, full.names =TRUE)
}

## filter gtf with deg file(s) 
lapply(fls, function(f) {
    deg <- read.delim(f)
    fgtf <- gff[which(gff$gene_id%in%deg$Id),]
    if (is.null(outPrefix) || length(outPrefix) == 0) {
        export(fgtf, paste(paste0(file_path_sans_ext(f), ".", gff_ext), sep=.Platform$file.sep))
    } else
    {
        export(fgtf, paste(paste0(outPrefix, file_path_sans_ext(f), ".", gff_ext), sep=.Platform$file.sep))
    }
})

# R session
info(logger, sessionInfo())

## execution time: end
T2<-Sys.time()
Tdiff= difftime(T2, T1)
info(logger, paste("Execution time : ", Tdiff,"seconds\n", sep=""))

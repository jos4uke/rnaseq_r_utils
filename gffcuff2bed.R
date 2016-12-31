# convert gff/gtf to bed for fasta sequence extraction using bedtools getfasta
# custom fragment length extraction using the 2 variables: flankUpstreamLenMax and flankDownstreamLenMax relative to the cotranscript coordinates
library(data.table)
library(tools)

# read gff
gff <- "cuffL1gtfTE.final.gtf"
gff_no_ext <- file_path_sans_ext(gff) 
gff.dt <- fread(gff, sep = "\t", header = T, stringsAsFactors = T, strip.white = T, quote = "")

# parse attribute field in gff and get transcript_id value
gff.attributes <- gff.dt$attribute
gff.attributes.nq <- gsub('"', '', gff.attributes)

getAttributeField <- function (x, field, attrsep = ";") {
  s = strsplit(x, split = attrsep, fixed = TRUE)
  sapply(s, function(atts) {
    a = strsplit(trimws(atts), split = " ", fixed = TRUE)
    m = match(field, sapply(a, "[", 1))
    if (!is.na(m)) {
      rv = a[[m]][2]
    }
    else {
      rv = as.character(NA)
    }
    return(rv)
  })
}

transcript_ids <- getAttributeField(gff.attributes.nq, "transcript_id")

# create a bed
## if cotranscript is strand(+) 
## then chromStart equals to s=startco-1-flankUpstreamLenMax or 0 if s<0
## and chromEnd equals to e=endco-1+flankDownstreamLenMax or endall if e>endall because cannot access chromosome max size
## if cotranscript is strand(-)
## then chromStart equals to s=startco-1-flankDownstreamLenMax or 0 if s<0
## and chromEnd equals to e=endco-1+flankUpstreamLenMax or endall if e>endall because cannot access chromosome max size
flankUpstreamLenMax <- 1000
flankDownstreamLenMax <- 2000
bed.dt <- data.table(chrom=gff.dt$`#chrm`,
                     chromStart=ifelse(gff.dt$strand=="+", 
                                       ifelse((s=gff.dt$startco-1-flankUpstreamLenMax)>=0,s,0), 
                                       ifelse((s=gff.dt$startco-1-flankDownstreamLenMax)>=0,s,0)),
                     chromEnd=ifelse(gff.dt$strand=="+", 
                                     ifelse((e=gff.dt$endco-1+flankDownstreamLenMax)<=gff.dt$endall, e, gff.dt$endall), 
                                     ifelse((e=gff.dt$endco-1+flankUpstreamLenMax)<=gff.dt$endall, e, gff.dt$endall)) , # how to get the max end of the chromosome? chrom size file? use endall as a max limit
                     name=paste(gff.dt$TE, transcript_ids, sep = "__"),
                     score=gff.dt$score,
                     strand=gff.dt$strand)

# export bed
fwrite(bed.dt, file=paste0(gff_no_ext, "_UpStrMax",flankUpstreamLenMax, "_DownStrMax", flankDownstreamLenMax,".bed"), sep="\t", quote = FALSE)

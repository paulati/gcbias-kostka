# generate a single file from these three dowloads
# sno/miRNA: http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/wgRna.txt.gz
# tRNAs: http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/tRNAs.txt.gz
# GENCODE Pseudogenes: http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/wgEncodeGencodePseudoGeneV19.txt.gz

script.directory.base <- "/paula/2018/gcbias/scripts"
output.directory.base <- "/paula/2018/gcbias/data/preparation/hg19/rnaGene/split_chr"
input.directory.base <- "/paula/2018/gcbias/data/download/rnaGene"

combine.input.directory.base <- "/paula/2018/gcbias/data/preparation/hg19/rnaGene/split_chr"
combine.output.directory.base <- "/paula/2018/gcbias/data/preparation/hg19/rnaGene"

setwd(script.directory.base)
source("0_common.R")

split.data.by.chr <- function(data.trna.bed, data.pseudogenes.bed, data.miRNA.bed, output.directory.base)
{
  
  out.file.name <- "tRNAs.bed"
  col.name <- "chrom"
  data.bed <- data.trna.bed
  split.data.bed(data.bed, col.name, out.file.name, output.directory.base)  
  
  out.file.name <- "wgEncodeGencodePseudoGeneV19.bed"
  col.name <- "chrom"
  data.bed <- data.pseudogenes.bed
  split.data.bed(data.bed, col.name, out.file.name, output.directory.base)  
  
  out.file.name <- "wgRna.bed"
  col.name <- "chrom"
  data.bed <- data.miRNA.bed
  split.data.bed(data.bed, col.name, out.file.name, output.directory.base)  
  
}

collapse.rna.by.chr <- function()
{
  
  setwd(input.directory.base)
  
  data.miRNA <- read.table("wgRna.txt", sep="\t")
  
  col.names.miRNA <- c("bin", "chrom", "chromStart", "chromEnd", "name", "score", "strand", 
                 "thickStart", "thickEnd", "type")
  
  colnames(data.miRNA) <- col.names.miRNA
  
  
  data.trna <- read.table("tRNAs.txt", sep="\t")
  
  col.names.trna <- c("bin", "chrom", "chromStart", "chromEnd", "name", "score", "strand", 
    "aa", "ac", "intron", "trnaScore", "genomeUrl", "trnaUrl")
  
  colnames(data.trna) <- col.names.trna
  
  data.pseudogenes <- read.table("wgEncodeGencodePseudoGeneV19.txt", sep="\t")
  
  col.names.pseudogenes <- c("bin", "name", "chrom", "strand", "txStart", "txEnd", 
                             "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds", 
                             "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames")
  
  colnames(data.pseudogenes) <- col.names.pseudogenes
  
  data.miRNA.bed <- data.miRNA[, c("chrom", "chromStart", "chromEnd")]
  data.trna.bed <- data.trna[, c("chrom", "chromStart", "chromEnd")]
  data.pseudogenes.bed <- data.pseudogenes[, c("chrom", "txStart", "txEnd")]
  
  split.data.by.chr(data.trna.bed, data.pseudogenes.bed, data.miRNA.bed, output.directory.base)
  

}

#combine 3 files in 1 rna gene bed file
combine.rna.data <- function()
{
  
  for(chr in c(1:22))
  {
    setwd(combine.input.directory.base)
    
    data.tRNA.file.name <- paste("chr", chr, "_tRNAs.bed", sep="")
    data.tRNA <- read.table(data.tRNA.file.name, sep="\t")
    
    data.pseudogene.file.name <- paste("chr", chr, "_wgEncodeGencodePseudoGeneV19.bed", sep="")
    data.pseudogene <- read.table(data.pseudogene.file.name, sep="\t")
    
    data.wgrna.data.wgrna  <- paste("chr", chr, "_wgRna.bed", sep="")
    data.wgrna <- read.table(data.wgrna.data.wgrna, sep="\t")
    
    data.all <- rbind(data.tRNA, data.pseudogene, data.wgrna)
    col.names <- c("chr", "start", "end")
    colnames(data.all) <- col.names
    
    data.all.bed <- bed.object(data.all)
    
    out.file.name <- paste("chr", chr, "_rnaGene.bed", sep="")
    setwd(combine.output.directory.base)
    write.table(data.all.bed, out.file.name, sep="\t", row.names=FALSE, col.names= FALSE, quote = FALSE) 
    
  }
  
}



# 1
# collapse.rna.by.chr()

# 2
# combine.rna.data()

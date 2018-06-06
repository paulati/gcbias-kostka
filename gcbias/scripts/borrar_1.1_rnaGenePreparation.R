

input.directory.base <- "/paula/2018/gcbias/data/download/hg19/rnaGene"
output.directory.base <- "/paula/2018/gcbias/data/preparation/hg19/rnaGene"
script.directory.base <- "/paula/2018/gcbias/scripts"

setwd(script.directory.base)
source("0_common.R")

#combine 3 files in 1 rna gene bed file
combine.rna.data <- function()
{

  for(chr in c(1:22))
  {
    setwd(input.directory.base)
    
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
    setwd(output.directory.base)
    write.table(data.all.bed, out.file.name, sep="\t", row.names=FALSE, col.names= FALSE, quote = FALSE) 
  
  }

}

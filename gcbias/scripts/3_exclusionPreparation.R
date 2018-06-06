
conf.download.directory.base <- "/paula/2018/gcbias/data/download"

conf.input.directory.base <- "/paula/2018/gcbias/data/preparation/hg19"
conf.output.directory.base <- "/paula/2018/gcbias/data/preparation/hg19"
conf.script.directory.base <- "/paula/2018/gcbias/scripts"

source("/paula/2018/gcbias/scripts/0_common.R")

set.input.directory <- function(input.directory.base, folder.name)
{
  input.directory <- paste(input.directory.base, folder.name, sep="/")  
  return(input.directory)
}

read.bed.data <- function(input.directory.base, folder.name, file.name.tail, chr)
{
  input.directory <- set.input.directory(input.directory.base, folder.name)
  print(input.directory)
  setwd(input.directory)
  data.file.name <- paste("chr", chr, file.name.tail, sep="")
  data <- read.table(data.file.name, sep="\t")
  return(data)
}

split.by.chr.exclusion.sets <- function(download.directory.base, input.directory.base, output.directory.base)
{
  base.path <- download.directory.base
  
  #all_mrna
  file.name <- "all_mrna.bed" 
  data <- bed.data(base.path, file.name)
  col.name <- "chr"
  out.file.name <- file.name
  out.directory <- paste(output.directory.base, "/all_mrna", sep="")
  split.data.bed(data, col.name, out.file.name, out.directory)

  #ensGene
  file.name <- "hg19_ensGene.bed" 
  data <- bed.data(base.path, file.name)
  col.name <- "chr"
  out.file.name <- file.name
  out.directory <- paste(output.directory.base, "/ensGene", sep="")
  split.data.bed(data, col.name, out.file.name, out.directory)
  
  #knownGene
  file.name <- "hg19.knownGene.bed" 
  data <- bed.data(base.path, file.name)
  col.name <- "chr"
  out.file.name <- file.name
  out.directory <- paste(output.directory.base, "/knownGene", sep="")
  split.data.bed(data, col.name, out.file.name, out.directory)
  
  #refGene
  file.name <- "hg19.refGene.bed" 
  data <- bed.data(base.path, file.name)
  col.name <- "chr"
  out.file.name <- file.name
  out.directory <- paste(output.directory.base, "/refGene", sep="")
  split.data.bed(data, col.name, out.file.name, out.directory)
  
  #transMapEnsembl
  file.name <- "hg19.transMapEnsembl.bed" 
  data <- bed.data(base.path, file.name)
  col.name <- "chr"
  out.file.name <- file.name
  out.directory <- paste(output.directory.base, "/transMapEnsembl", sep="")
  split.data.bed(data, col.name, out.file.name, out.directory)
  
  #transMapRefSeq
  file.name <- "hg19.transMaprefSeqV4.bed" 
  data <- bed.data(base.path, file.name)
  col.name <- "chr"
  out.file.name <- file.name
  out.directory <- paste(output.directory.base, "/transMapRefSeq", sep="")
  split.data.bed(data, col.name, out.file.name, out.directory)
  
  #transMapRNA
  file.name <- "hg19.transMapRnaV4.bed" 
  data <- bed.data(base.path, file.name)
  col.name <- "chr"
  out.file.name <- file.name
  out.directory <- paste(output.directory.base, "/transMapRNA", sep="")
  split.data.bed(data, col.name, out.file.name, out.directory)
  
  #_rnaGene.bed
  #split with script 1_joinRNAGenePreparation.R
  
  #pseudogeneHuman
  base.path <- paste(input.directory.base, "/pseudogeneHuman", sep="")
  setwd(base.path)
  file.name <- "pseudogeneHuman74.bed" 
  data <- bed.data(base.path, file.name)
  col.name <- "chr"
  out.file.name <- file.name
  out.directory <- paste(output.directory.base, "/pseudogeneHuman", sep="")
  split.data.bed(data, col.name, out.file.name, out.directory)
  
}

merge.exclusion.sets <- function(input.directory.base)
{
  for(chr in c(1:22))
  {
    folder.name <- "all_mrna"
    file.name.tail <- "_all_mrna.bed"
    data.all_mrna <- read.bed.data(input.directory.base, folder.name, file.name.tail, chr)
    data.all_mrna <- data.all_mrna[, c(1:3)]
    
    folder.name <- "ensGene"
    file.name.tail <- "_hg19_ensGene.bed"
    data.ensGene <- read.bed.data(input.directory.base, folder.name, file.name.tail, chr)
    data.ensGene <- data.ensGene[, c(1:3)]
    
    folder.name <- "knownGene"
    file.name.tail <- "_hg19.knownGene.bed"
    data.knownGene <- read.bed.data(input.directory.base, folder.name, file.name.tail, chr)  
    data.knownGene <- data.knownGene[, c(1:3)]
    
    folder.name <- "pseudogeneHuman"
    file.name.tail <- "_pseudogeneHuman74.bed"
    data.pseudogeneHuman <- read.bed.data(input.directory.base, folder.name, file.name.tail, chr)  
    data.pseudogeneHuman <- data.pseudogeneHuman[, c(1:3)]
    
    folder.name <- "refGene"
    file.name.tail <- "_hg19.refGene.bed"
    data.refGene <- read.bed.data(input.directory.base, folder.name, file.name.tail, chr)  
    data.refGene <- data.refGene[, c(1:3)]
    
    folder.name <- "rnaGene"
    file.name.tail <- "_rnaGene.bed"
    data.rnaGene <- read.bed.data(input.directory.base, folder.name, file.name.tail, chr)  
    data.rnaGene <- data.rnaGene[, c(1:3)]
    
    folder.name <- "transMapEnsembl"
    file.name.tail <- "_hg19.transMapEnsembl.bed"
    data.transMapEnsembl <- read.bed.data(input.directory.base, folder.name, file.name.tail, chr)  
    data.transMapEnsembl <- data.transMapEnsembl[, c(1:3)]
    
    folder.name <- "transMapRefSeq"
    file.name.tail <- "_hg19.transMaprefSeqV4.bed"
    data.transMapRefSeq <- read.bed.data(input.directory.base, folder.name, file.name.tail, chr)  
    data.transMapRefSeq <- data.transMapRefSeq[, c(1:3)]
    
    folder.name <- "transMapRNA"
    file.name.tail <- "_hg19.transMapRnaV4.bed"
    data.transMapRNA <- read.bed.data(input.directory.base, folder.name, file.name.tail, chr)  
    data.transMapRNA <- data.transMapRNA[, c(1:3)]
    
    data.all <- rbind(data.all_mrna, data.ensGene, data.knownGene,
                      data.pseudogeneHuman, data.refGene, data.rnaGene, 
                      data.transMapEnsembl, data.transMapRefSeq, data.transMapRNA)
    col.names <- c("chr", "start", "end")
    colnames(data.all) <- col.names
    
    data.all.bed <- bed.object(data.all)
    
    out.file.name <- paste("chr", chr, "_exclusion.bed", sep="")
    output.directory.base <- paste(conf.output.directory.base, "/exclusion_filter", sep="")
    setwd(output.directory.base)
    write.table(data.all.bed, out.file.name, sep="\t", row.names=FALSE, col.names= FALSE, quote = FALSE) 
    
  }  
}


#hay que splitear cada uno de los conjuntos de exclusion 
#split.by.chr.exclusion.sets(conf.download.directory.base, conf.input.directory.base, conf.output.directory.base)

#y despues mergearlos con las lineas de abajo
merge.exclusion.sets(conf.input.directory.base)

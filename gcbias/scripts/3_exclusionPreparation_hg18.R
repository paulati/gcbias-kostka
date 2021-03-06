
conf.download.directory.base <- "/paula/2018/gcbias/data/download"
conf.input.directory.base <- "/paula/2018/gcbias/data/preparation/hg18"
conf.output.directory.base <- "/paula/2018/gcbias/data/preparation/hg18"


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

split.by.chr.exclusion.set <- function(file.name, folder.name, base.path, output.directory.base)
{
  setwd(base.path)
  data <- bed.data(base.path, file.name)
  col.name <- "chr"
  out.file.name <- file.name
  out.directory <- paste(output.directory.base, folder.name, sep="")
  split.data.bed(data, col.name, out.file.name, out.directory)
}

clean.by.chr.exclusion.set <- function(file.name.base, folder.name, base.path, output.directory.base)
{
  setwd(base.path)
  for(i in c(1:22))
  {
    file.name <- paste("chr", i, file.name.base, sep="")
    data <- bed.data(base.path, file.name)
    col.name <- "chr"
    col.value <- paste("chr", i, sep="")
    out.file.name <- file.name
    out.directory <- paste(output.directory.base, folder.name, sep="")
    print(out.directory)
    clean.data.bed(data, col.name, col.value, out.file.name, out.directory)
    #split.data.bed(data, col.name, out.file.name, out.directory)
  }
}

split.by.chr.exclusion.sets <- function(download.directory.base, input.directory.base, output.directory.base)
{
  base.path <- download.directory.base
  
  #all_mrna
  file.name <- "all_mrna.bed" 
  folder.name <- "/all_mrna"
  split.by.chr.exclusion.set(file.name, folder.name, base.path, output.directory.base)
  
  #ensGene
  file.name <- "hg18_ensGene.bed" 
  folder.name <- "/ensGene"
  split.by.chr.exclusion.set(file.name, folder.name, base.path, output.directory.base)
  
  #knownGene
  file.name <- "hg18_knownGene.bed" 
  folder.name <- "/knownGene"
  split.by.chr.exclusion.set(file.name, folder.name, base.path, output.directory.base)

  #refGene
  file.name <- "hg18_refGene.bed" 
  folder.name <- "/refGene"
  split.by.chr.exclusion.set(file.name, folder.name, base.path, output.directory.base)
  

  #los trans fueron generados con liftover
  
  #transMapEnsembl
  #quito los elementos con "random" en el nombre del chr
  base.path.transMapEnsembl <- paste(output.directory.base, "/transMapEnsembl/lift", sep="")
  file.name.base <- "_hg18.transMapEnsembl.bed"
  folder.name <- "/transMapEnsembl"
  clean.by.chr.exclusion.set(file.name.base, folder.name, base.path.transMapEnsembl, output.directory.base)
  
  #transMapRefSeq
  #quito los elementos con "random" en el nombre del chr
  base.path.transMapRefSeq <- paste(output.directory.base, "/transMapRefSeq/lift", sep="")
  file.name.base <- "_hg18.transMaprefSeqV4.bed"
  folder.name <- "/transMapRefSeq"
  clean.by.chr.exclusion.set(file.name.base, folder.name, base.path.transMapRefSeq, output.directory.base)

  #transMapRNA
  #quito los elementos con "random" en el nombre del chr
  base.path.transMapRNA <- paste(output.directory.base, "/transMapRNA/lift", sep="")
  file.name.base <- "_hg18.transMapRnaV4.bed"
  folder.name <- "/transMapRNA"
  clean.by.chr.exclusion.set(file.name.base, folder.name, base.path.transMapRNA, output.directory.base)

  #_rnaGene.bed
  file.name <- "hg18_rnaGene.bed"
  folder.name <- "/rnaGene"
  split.by.chr.exclusion.set(file.name, folder.name, base.path, output.directory.base)
  
  #pseudogeneHuman
  file.name <- "pseudogeneHuman60.bed" 
  folder.name <- "/pseudogeneHuman"
  base.path.pseudogene <- paste(input.directory.base, "/pseudogeneHuman", sep="")
  split.by.chr.exclusion.set(file.name, folder.name, base.path.pseudogene, output.directory.base)
  
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
    file.name.tail <- "_hg18_ensGene.bed"
    data.ensGene <- read.bed.data(input.directory.base, folder.name, file.name.tail, chr)
    data.ensGene <- data.ensGene[, c(1:3)]
    
    folder.name <- "knownGene"
    file.name.tail <- "_hg18_knownGene.bed"
    data.knownGene <- read.bed.data(input.directory.base, folder.name, file.name.tail, chr)  
    data.knownGene <- data.knownGene[, c(1:3)]
    
    folder.name <- "pseudogeneHuman"
    file.name.tail <- "_pseudogeneHuman60.bed"
    data.pseudogeneHuman <- read.bed.data(input.directory.base, folder.name, file.name.tail, chr)  
    data.pseudogeneHuman <- data.pseudogeneHuman[, c(1:3)]
    
    folder.name <- "refGene"
    file.name.tail <- "_hg18_refGene.bed"
    data.refGene <- read.bed.data(input.directory.base, folder.name, file.name.tail, chr)  
    data.refGene <- data.refGene[, c(1:3)]
    
    folder.name <- "rnaGene"
    file.name.tail <- "_hg18_rnaGene.bed"
    data.rnaGene <- read.bed.data(input.directory.base, folder.name, file.name.tail, chr)  
    data.rnaGene <- data.rnaGene[, c(1:3)]
    
    folder.name <- "transMapEnsembl"
    file.name.tail <- "_hg18.transMapEnsembl.bed"
    data.transMapEnsembl <- read.bed.data(input.directory.base, folder.name, file.name.tail, chr)  
    data.transMapEnsembl <- data.transMapEnsembl[, c(1:3)]
    
    folder.name <- "transMapRefSeq"
    file.name.tail <- "_hg18.transMaprefSeqV4.bed"
    data.transMapRefSeq <- read.bed.data(input.directory.base, folder.name, file.name.tail, chr)  
    data.transMapRefSeq <- data.transMapRefSeq[, c(1:3)]
    
    folder.name <- "transMapRNA"
    file.name.tail <- "_hg18.transMapRnaV4.bed"
    data.transMapRNA <- read.bed.data(input.directory.base, folder.name, file.name.tail, chr)  
    data.transMapRNA <- data.transMapRNA[, c(1:3)]
    
    data.all <- rbind(data.all_mrna, data.ensGene, data.knownGene,
                      data.pseudogeneHuman, data.refGene, data.rnaGene, 
                      data.transMapEnsembl, data.transMapRefSeq, data.transMapRNA)
    col.names <- c("chr", "start", "end")
    colnames(data.all) <- col.names
    
    #hay algunos que dicen chr21_random
    # tmp <- is.valid.region(
    #   data.all,
    #   check.zero.based = TRUE,
    #   check.chr = TRUE,
    #   throw.error = FALSE,
    #   verbose = TRUE
    # )
    # 
    # indexes <- which(! tmp)
    # 
    # not.valid <- data.all[indexes, ]
    # 
  
    data.all.bed <- bed.object(data.all)
    
    out.file.name <- paste("chr", chr, "_exclusion.bed", sep="")
    output.directory.base <- paste(conf.output.directory.base, "/exclusion_filter", sep="")
    setwd(output.directory.base)
    write.table(data.all.bed, out.file.name, sep="\t", row.names=FALSE, col.names= FALSE, quote = FALSE) 
    
  }  
}


#hay que splitear cada uno de los conjuntos de exclusion 
# split.by.chr.exclusion.sets(conf.download.directory.base, conf.input.directory.base, conf.output.directory.base)

#y despues mergearlos con las lineas de abajo
# merge.exclusion.sets(conf.input.directory.base)
# 
# 
# download.directory.base <- conf.download.directory.base
# input.directory.base <- conf.input.directory.base
# output.directory.base <- conf.output.directory.base





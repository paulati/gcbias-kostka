input.directory.base <- "/paula/2018/gcbias/data/download"
output.directory.base <- "/paula/2018/gcbias/data/preparation"
script.directory.base <- "/paula/2018/gcbias/scripts"

setwd(script.directory.base)
source("0_common.R")
source("1.1_rnaGenePreparation.R")

#hg19_all_mrna.bed
setwd(input.directory.base)
data.bed <- read.table("hg19_all_mrna.bed", sep="\t")
col.names <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", 
               "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts")
colnames(data.bed) <- col.names
col.name <- "chrom"
out.file.name <- "hg19_all_mrna.bed"
out.directory <- paste(output.directory.base, "/all_mrna", sep="")
split.data.bed(data.bed, col.name, out.file.name, out.directory)

#hg19_ensGene.bed
setwd(input.directory.base)
data.bed <- read.table("hg19_ensGene.bed", sep="\t")
col.names <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", 
               "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts")
colnames(data.bed) <- col.names
col.name <- "chrom"
out.file.name <- "hg19_ensGene.bed"
out.directory <- paste(output.directory.base, "/ensGene", sep="")
split.data.bed(data.bed, col.name, out.file.name, out.directory)


#hg19_knownGene.bed
setwd(input.directory.base)
data.bed <- read.table("hg19_knownGene.bed", sep="\t")
col.names <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", 
               "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts")
colnames(data.bed) <- col.names
col.name <- "chrom"
out.file.name <- "hg19_knownGene.bed"
out.directory <- paste(output.directory.base, "/knownGene", sep="")
split.data.bed(data.bed, col.name, out.file.name, out.directory)


#hg19_refGene.bed
setwd(input.directory.base)
data.bed <- read.table("hg19_refGene.bed", sep="\t")
col.names <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", 
               "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts")
colnames(data.bed) <- col.names
col.name <- "chrom"
out.file.name <- "hg19_refGene.bed"
out.directory <- paste(output.directory.base, "/refGene", sep="")
split.data.bed(data.bed, col.name, out.file.name, out.directory)

#hg19_transMapEnsembl.bed
setwd(input.directory.base)
data.bed <- read.table("hg19_transMapEnsemblV4.bed", sep="\t")
col.names <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", 
               "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts")
colnames(data.bed) <- col.names
col.name <- "chrom"
out.file.name <- "hg19_transMapEnsemblV4.bed"
out.directory <- paste(output.directory.base, "/transMapEnsembl", sep="")
split.data.bed(data.bed, col.name, out.file.name, out.directory)


#hg19_transMapRefSeqV4.bed
setwd(input.directory.base)
data.bed <- read.table("hg19_transMapRefSeqV4.bed", sep="\t")
col.names <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", 
               "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts")
colnames(data.bed) <- col.names
col.name <- "chrom"
out.file.name <- "hg19_transMapRefSeqV4.bed"
out.directory <- paste(output.directory.base, "/transMapRefSeq", sep="")
split.data.bed(data.bed, col.name, out.file.name, out.directory)

#hg19_transMapRNA.bed
setwd(input.directory.base)
data.bed <- read.table("hg19_transMapRnaV4.bed", sep="\t")
col.names <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", 
               "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts")
colnames(data.bed) <- col.names
col.name <- "chrom"
out.file.name <- "hg19_transMapRNA.bed"
out.directory <- paste(output.directory.base, "/transMapRNA", sep="")
split.data.bed(data.bed, col.name, out.file.name, out.directory)

#pseudogene
setwd(input.directory.base)
data.bed <- read.table("pseudogeneHuman74.txt", sep="\t", header = TRUE)
col.name <- "Chromosome"
data.bed <- data.bed[, c("Chromosome", "Start.Coordinate", "Stop.Coordinate" )]
data.bed[, "Chromosome"] <- sapply(data.bed$Chromosome, function(x) paste("chr", x, sep=""))
out.file.name <- "pseudogeneHuman74.bed"
out.directory <- paste(output.directory.base, "/pseudogeneHuman", sep="")
split.data.bed(data.bed, col.name, out.file.name, out.directory)



conf.input.directory.base <- "/paula/2018/gcbias/data/download/neutral_model"
conf.output.directory.base <- "/paula/2018/gcbias/data/preparation/hg18/neutral_model/gencode.v3c.annotation.NCBI36"
conf.input.file.name <- "gencode.v3c.annotation.NCBI36.level_3.gtf"
conf.version <- "hg18"

source("/paula/2018/gcbias/scripts/0_common.R")

split.gencode.features <- function(input.directory.base,input.file.name, output.directory.base, version)
{
  setwd(input.directory.base)
  data <- read.table(input.file.name, sep="\t", skip = 5)
  data$V1 <- paste(version, ".", data$V1, sep="")
  out.file.name <- input.file.name
  split.data.frame(data, "V1", out.file.name, output.directory.base, version)
    
  
  # col.name <- "Chromosome"
  # data.bed <- data.bed[, c("Chromosome", "Start.Coordinate", "Stop.Coordinate" )]
  # data.bed[, "Chromosome"] <- sapply(data.bed$Chromosome, function(x) paste("chr", x, sep=""))
  # out.file.name <- "pseudogeneHuman74.bed"
  # out.directory <- paste(output.directory.base, "/pseudogeneHuman", sep="")
  # split.data.bed(data.bed, col.name, out.file.name, out.directory)



}


split.gencode.features(conf.input.directory.base,
                       conf.input.file.name, 
                       conf.output.directory.base,
                       conf.version)




# input.directory.base <- conf.input.directory.base
# input.file.name <- conf.input.file.name
# output.directory.base <- conf.output.directory.base
# version <- conf.version


source("/paula/2018/gcbias/scripts/0_common.R")

library(stringr)

# 
# elements.base.path <- config.elements.base.path
# elements.file.name <- config.elements.file.name
# elements.chr.flank.base.path <- config.elements.chr.flank.base.path
# elements.chr.flank.base.file.name <- config.elements.chr.flank.base.file.name
# elements.chr.base.path <- config.elements.chr.base.path
# elements.chr.base.file.name <- config.elements.chr.base.file.name


split.elements.and.flanking.beds <- function(elements.base.path, elements.file.name,
                                             elements.chr.flank.base.path, elements.chr.flank.base.file.name,
                                             elements.chr.base.path, elements.chr.base.file.name)
{
  
  setwd(elements.base.path)
  elements <- read.csv(elements.file.name, sep="\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(elements) <- c("chr", "start", "end", "name")
  
  elements[, "start.kostka"] <- elements$start - rep(500000, nrow(elements))
  elements[, "end.kostka"] <- elements$end + rep(500000, nrow(elements))
  
  kostka.regions <- elements[, c("chr", "start.kostka", "end.kostka", "name")]
  setwd(elements.chr.flank.base.path)
  for (i in c(1:nrow(kostka.regions)))
  {
    chr <- kostka.regions[i, "chr"]
    name <- kostka.regions[i, "name"]
    file.name <- paste(chr, "_", name, elements.chr.flank.base.file.name, sep="")
    data <- kostka.regions[i, ]
    write.table(data, file.name, quote=FALSE, col.names = FALSE, row.names = FALSE, sep="\t")
  }
  
  elements.regions <- elements[, c("chr", "start", "end", "name")]
  setwd(elements.chr.base.path)
  for(i in c(1:nrow(kostka.regions)))
  {
    chr <- elements.regions[i, "chr"]
    name <- elements.regions[i, "name"]
    file.name <- paste(chr, "_", name, elements.chr.base.file.name, sep="")
    data <- elements.regions[i, ]
    write.table(data, file.name, quote=FALSE, col.names = FALSE, row.names = FALSE, sep="\t")
  }
  
  
  #separo las regiones de kostka por cromosoma
  setwd(elements.chr.flank.base.path)
  for(chr in 1:22)
  {
    chr.name <- paste("chr", chr, sep="")
    file.name <- paste(chr.name, elements.chr.flank.base.file.name, sep="")
    data <- kostka.regions[kostka.regions$chr == chr.name, ]
    write.table(data, file.name, quote=FALSE, col.names = FALSE, row.names = FALSE, sep="\t")
  }

  elements.regions <- elements[, c("chr", "start", "end", "name")]
  setwd(elements.chr.base.path)
  for(chr in 1:22)
  {
    chr.name <- paste("chr", chr, sep="")
    file.name <- paste(chr.name, elements.chr.base.file.name, sep="")
    data <- elements.regions[elements.regions$chr == chr.name, ]
    write.table(data, file.name, quote=FALSE, col.names = FALSE, row.names = FALSE, sep="\t")
  }
}

split.elements.and.flanking.beds.by.chr <- function(elements.base.path, elements.file.name,
                                             elements.chr.flank.base.path, elements.chr.flank.base.file.name,
                                             elements.chr.base.path, elements.chr.base.file.name)
{
  
  setwd(elements.base.path)
  elements <- read.csv(elements.file.name, sep="\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(elements) <- c("chr", "start", "end", "name")
  
  elements[, "start.kostka"] <- elements$start - rep(500000, nrow(elements))
  elements[, "end.kostka"] <- elements$end + rep(500000, nrow(elements))
  
  kostka.regions <- elements[, c("chr", "start.kostka", "end.kostka", "name")]
 
  setwd(elements.chr.flank.base.path)
  for(chr in 1:22)
  {
    chr.name <- paste("chr", chr, sep="")
    file.name <- paste(chr.name, elements.chr.flank.base.file.name, sep="")
    data <- kostka.regions[kostka.regions$chr == chr.name, ]
    write.table(data, file.name, quote=FALSE, col.names = FALSE, row.names = FALSE, sep="\t")
  }

  elements.regions <- elements[, c("chr", "start", "end", "name")]
  setwd(elements.chr.base.path)
  for(chr in 1:22)
  {
    chr.name <- paste("chr", chr, sep="")
    file.name <- paste(chr.name, elements.chr.base.file.name, sep="")
    data <- elements.regions[elements.regions$chr == chr.name, ]
    write.table(data, file.name, quote=FALSE, col.names = FALSE, row.names = FALSE, sep="\t")
  }
}


split.flanking.beds.by.elements <- function(input.base.path, file.name.tail, output.base.path)
{
  
  for(chr in c(1:22))
  {
  
    file.name <- paste("chr", chr, file.name.tail, sep="")
    data.chr <- bed.data.element(input.base.path, file.name)
    
    if(! is.null(data.chr))
    {
      q <- nrow(data.chr)
      
      setwd(output.base.path)
      for(k in c(1:q))
      {
        element.name <- data.chr$name[k]
        out.file.name <- paste("chr", chr, "_", element.name, "_flnk.bed", sep="")
        write.table(data.chr[k, ], out.file.name, quote=FALSE, col.names = FALSE, row.names = FALSE, sep="\t")
      }
      
    }
  
  }
}



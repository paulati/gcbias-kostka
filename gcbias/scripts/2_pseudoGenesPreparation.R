
library(stringr)

convert.to.bed <- function(input.directory.base, input.file.name, 
                           output.directory.base, output.file.name)
{
  setwd(input.directory.base)
  
  data <-read.csv(input.file.name, sep="\t", header=TRUE)
  
  result <- data.frame(matrix(ncol=3, nrow=0))
  colnames(result) <- c("chr", "start", "end")
  
  pattern.exons <- "\\[[:digit:]+\\,[:blank:][:digit:]+\\]"
  pattern.exon <- "[:digit:]+"
  
  result.index <- 1
  
  for(i in 1:nrow(data))
  {
    chr <- as.character(data[i, "Chromosome"])
    
    exons <- as.character(data[i, "Exons"])
    
    #exons <- "[[67580013, 67580233], [67587393, 67587962], [67640672, 67640842]]"
    print(exons)
  
    exons.split <- unlist(str_extract_all(exons, pattern.exons, simplify = FALSE))
  
    for(exon in exons.split)
    {
      #exon <- "[5165490, 5165633]"
      print(exon)
      
      start.end <- unlist(str_extract_all(exon, pattern.exon, simplify = FALSE))
      start <- as.numeric(start.end[1])
      end <- as.numeric(start.end[2])
      
      print(start)
      print(end)
      
      result[result.index, "chr"] <- paste("chr", chr, sep="")
      result[result.index, "start"] <- start
      result[result.index, "end"] <- end
      result.index <- result.index + 1
      
    }
    
  }

  setwd(output.directory.base)
  write.table(result, output.file.name, sep="\t", col.names = FALSE, row.names = FALSE, quote=FALSE)
  
}

#armo los archivos de pseudogene por chr

#setwd("~/Documentos/r/2017/lara/data/exclusion/split_hg18_pseudogeneHuman60/")
#for(j in 1:22) 
#{
#  chr.name <- paste("chr", j, sep="")
#  data.file <- result[result$chr == chr.name, ]
#  file.name <- paste("chr", j,  "_pseudogeneHuman60.bed", sep="")
#  write.table(x = data.file, file = file.name, sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#}


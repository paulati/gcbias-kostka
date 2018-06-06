library(bedr)

read.data <- function(base.path, file.name)
{
  setwd(base.path)
  data <- tryCatch({
    if (file.size(file.name) > 0){
      read.csv(file.name, sep="\t", header=FALSE)
    }
  }, error = function(err) {
    # error handler picks up where error was generated
    print(paste("Read.table didn't work!:  ",err))
  })
  
  if(! is.null(data))
  {
    colnames(data) <- c("chr", "start", "end")  
  }
  return(data)
}

bed.object <- function(data)
{
  data[, "chr"] <- as.character(data[, "chr"])
  col.format <- c("chr", "start", "end")
  data.format <- data[,col.format] 
  
  data.sort <- bedr.sort.region(data.format)
  #data[order("start"), ]
  
  data.merge <- bedr.merge.region(data.sort)
  
  data.bed <- convert2bed(data.merge,
                          set.type = TRUE,
                          check.zero.based = TRUE,
                          check.chr = TRUE,
                          check.valid = TRUE,
                          check.sort = TRUE,
                          check.merge = TRUE,
                          verbose = TRUE)
  return(data.bed)
  
}

subtract.exclusion.from.bed <- function(chr, data.A.base.path, data.A.file.name, 
                                        data.B.base.path, data.B.base.file.name,
                                        intersect.data.base.path)
{
  
  exclusion.file.name <- paste(chr, data.B.base.file.name, sep="")
  chr.exclusion.data <- read.data(data.B.base.path, exclusion.file.name)
  chr.exclusion.data.bed <- bed.object(chr.exclusion.data)
  
  elements.data <- read.data(data.A.base.path, data.A.file.name)
  
  if(! is.null(elements.data) )
  {
    
    elements.data.bed <- bed.object(elements.data)
    
    chr.alignment.data.util.bed <- bedr.subtract.region(
      #chr.alignment.data.util.bed <- bedr.join.region(
      elements.data.bed,
      chr.exclusion.data.bed,
      fraction.overlap = 1/1e9,
      remove.whole.feature = FALSE,
      check.zero.based = TRUE,
      check.chr = TRUE,
      check.valid = TRUE,
      check.sort = TRUE,
      check.merge = TRUE,
      verbose = TRUE
    )
    
    element.start <- elements.data.bed$start[1]
    element.end <- elements.data.bed$end[1]
    
    # dummy <- chr.exclusion.data.bed[(chr.exclusion.data.bed$start >= element.start
    #                                 & chr.exclusion.data.bed$start <= element.end) |
    #                                   (chr.exclusion.data.bed$end >= element.start & 
    #                                   chr.exclusion.data.bed$end <= element.end) |
    #                                   (chr.exclusion.data.bed$start <= element.start
    #                                 & chr.exclusion.data.bed$start >= element.end), ]
    
    
    setwd(intersect.data.base.path)
    out.file.name <- data.A.file.name
    write.table(x=chr.alignment.data.util.bed, file = out.file.name, col.names = FALSE, row.names = FALSE, quote=FALSE, sep = "\t")    
    
  }
}
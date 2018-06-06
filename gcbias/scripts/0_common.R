library(bedr)

split.data.bed <- function(data.bed, col.name, out.file.name, out.directory)
{
  version <- ""
  split.data.frame(data.bed, col.name, out.file.name, out.directory, version)
}

clean.data.bed <- function(data.bed, col.name, col.value, out.file.name, out.directory)
{
  setwd(out.directory)
  data.column <- data.bed[, col.name]
  indexes <- which(data.column == col.value)
  chr.data <- data.bed[indexes, ]
  write.table(chr.data, out.file.name, sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

split.data.frame <- function(data.bed, col.name, out.file.name, out.directory, version)
{
  setwd(out.directory)
  chr.name.prefix <- "chr"
  
  for(chr in c(1:22))
  {
    file.name <- paste("chr", chr, "_", out.file.name, sep="")  
    chr.name <- paste(chr.name.prefix, chr, sep="")
    if(version != "")
    {
      chr.name <- paste(version, ".", chr.name, sep="")  
    } 
    data.column <- data.bed[, col.name]
    indexes <- which(data.column == chr.name)
    chr.data <- data.bed[indexes, ]
    write.table(chr.data, file.name, sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }  
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


bed.data <- function(base.path, file.name)
{
  setwd(base.path)
  data <- read.table(file.name, sep="\t", header = FALSE, stringsAsFactors = FALSE)
  result <- data[, c(1:3)]
  colnames(result) <- c("chr", "start", "end")
  return(result)
}


bed.data.element <- function(base.path, file.name)
{
  setwd(base.path)
  data <- tryCatch({
    size <- file.size(file.name)
    if ( !is.na(size) & size > 0){
      read.table(file.name, sep="\t", header = FALSE, stringsAsFactors = FALSE)
    }
  }, error = function(err) {
    # error handler picks up where error was generated
    print(paste("Read.table didn't work!:  ",err))
  })
  
  if(! is.null(data))
  {
    result <- data[, c(1:4)]
    colnames(result) <- c("chr", "start", "end", "name")
  }
  else
  {
    result <- data
  }
  return(result)
  
  
}



get.species.lst <- function(ali.orgs.base.path, ali.orgs.file.name)
{
  setwd(ali.orgs.base.path)
  orgs <- scan(ali.orgs.file.name,what="character")
  return(orgs)
}




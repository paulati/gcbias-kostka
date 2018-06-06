source("/paula/2018/gcbias/scripts/0_common.R")

split.element.beds.by.elements <- function(input.base.path, file.name.tail, output.base.path)
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
        out.file.name <- paste("chr", chr, "_", element.name, "_elem.bed", sep="")
        write.table(data.chr[k, ], out.file.name, quote=FALSE, col.names = FALSE, row.names = FALSE, sep="\t")
      }
      
    }
    
  }
}

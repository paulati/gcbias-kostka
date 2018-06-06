
read.data <- function(base.path, file.name)
{
  
  setwd(base.path)
  
  
  col.classes <- c("character", "character", "character", "character", "character", 
                   "character", "character", "character", "character", "character", 
                   "character", "character", "character", "character", 
                   "logical", "logical", "logical")
  
  colnames(data)
  
  data <- read.table(file.name, sep="\t", header=TRUE, stringsAsFactors = FALSE, colClasses = col.classes)
  count <- nrow(data)
  data$kostka.class.label <- unlist(lapply(c(1: count), function(i) { 
    str <- data$class.x[i]
    pttr <- as.character(data$kostka.class[i])
    rpl <- ""
    result.replace <- str_replace(string = str, pattern = pttr, replacement = rpl)
    #quito |
    result <- substr(result.replace, 1, nchar(result.replace)-1)
    return(result)  }))
  
  return(data)
    
}



base.path <- "/paula/2018/kostka/data/results"
file.name <- "to_compare_prequel_testMa4_sneg_.csv"

data <- read.data(base.path, file.name)

obs <- data$kostka.class.label
pred  <- data$class.y

confusion.matrix.testMa4_sneg <- table(pred, obs)

######################################################################

base.path <- "/paula/2018/kostka/data/results"
file.name <- "to_compare_prequel_testMamail_.csv"

data <- read.data(base.path, file.name)

obs <- data$kostka.class.label
pred  <- data$class.y

confusion.matrix.testMamail <- table(pred, obs)

######################################################################

base.path <- "/paula/2018/kostka/data/results"
file.name <- "to_compare_test3_.csv"

data <- read.data(base.path, file.name)

obs <- data$kostka.class.label
pred  <- data$class.y

confusion.matrix.test3 <- table(pred, obs)

######################################################################

base.path <- "/paula/2018/kostka/data/results"
file.name <- "to_compare_test11_.csv"

data <- read.data(base.path, file.name)

obs <- data$kostka.class.label
pred  <- data$class.y

confusion.matrix.test11 <- table(pred, obs)


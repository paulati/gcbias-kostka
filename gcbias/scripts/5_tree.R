library(ape)
library(rphast)

create.tree <- function(input.tree.base.path, tree.file.name,
                        species.base.path, species.file.name,
                        output.tree.base.path) 
{
  setwd(input.tree.base.path)
  big.tree <- read.newick.tree(tree.file.name)
  
  setwd(species.base.path)
  species <- scan(species.file.name, what = 'character')
  
  usage.tree <- prune.tree(big.tree, species, all.but=TRUE)
  
  tree<-read.tree(text=usage.tree)  
  plot(tree)
  
  #borro las longitudes de las ramas
  tree$edge.length<-NULL
  #plot(tree)
  
  result <- write.tree(tree)
  
  setwd(output.tree.base.path)
  write(result, tree.file.name)
  
  #result <- sub(pattern = "[(][(][(]turTru2[,]orcOrc1[)][,]phyCat1[)][,]balAcu1[)]", 
  #              x = interest.node, 
  #              replacement = "\\2(((turTru2,orcOrc1),phyCat1),balAcu1)cetaceos\\3")
  
  #tmp <- read.tree(text=result)  
  #plot(tmp)
  
  return(result)
  
}

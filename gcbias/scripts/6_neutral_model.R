
library(rphast)

create.neutral.mod <- function(align.4d.base.folder, align.4d.file.name,
                               tree.file.base.folder, tree.file.name, 
                               neutral.mod.base.folder, neutral.mod.file.name) {
  
  
  setwd(align.4d.base.folder)
  align4d <- read.msa(align.4d.file.name, pointer.only=TRUE, format = "SS")
  
  setwd(tree.file.base.folder)
  tree <- read.newick.tree(tree.file.name)
  
  #neutralMod <- phyloFit(align4d, tree=tree, subst.mod="REV")
  neutralMod <- phyloFit(align4d, tree=tree, subst.mod="SSREV")    
  
  setwd(neutral.mod.base.folder)
  write.tm(neutralMod, neutral.mod.file.name)
  
}

concat <- function(in.align.4d.base.folder, align4d.file.name.base, out.align.4d.base.folder, out.align.4d.file.name)
{
  
  msas <- list()
  
  for (i in c(1:22) ) {
    
    align4d.file.name <- paste("chr", i, align4d.file.name.base, sep="")
    
    align4d.folder.path = paste(in.align.4d.base.folder, "/chr", i, "/", sep="")
    
    setwd(in.align.4d.base.folder)
    chr.4d  <- read.msa(align4d.file.name, pointer.only=FALSE, format = "SS")
    
    msas[[i]] <- chr.4d
    
  }
  
  result <- concat.msa(msas, ordered = FALSE, pointer.only = FALSE)
  
  setwd(out.align.4d.base.folder)
  format <- guess.format.msa(out.align.4d.file.name, method="extension")
  write.msa(result, out.align.4d.file.name, format)  
  
  
  
}



main.global <- function(in.align.4d.base.folder, in.align4d.file.name.base,
                        out.align.4d.base.folder, out.align.4d.file.name,
                        tree.file.base.folder, tree.file.name, 
                        neutral.mod.base.folder, neutral.mod.file.name)
{
  
  concat(in.align.4d.base.folder, in.align4d.file.name.base, out.align.4d.base.folder, out.align.4d.file.name)
  
  align.4d.base.folder <- out.align.4d.base.folder
  align.4d.file.name <- out.align.4d.file.name
  
  create.neutral.mod (align.4d.base.folder, align.4d.file.name,
                      tree.file.base.folder, tree.file.name, 
                      neutral.mod.base.folder, neutral.mod.file.name) 
  
  
  return(0)
}

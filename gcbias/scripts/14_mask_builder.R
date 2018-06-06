#este depende de los ancestrales generados con prequel
#probar como funciona con los hars de kostka

source("/paula/2018/gcbias/scripts/0_common.R")


maskCpG <- function(ali, names, anc.fa){
  #============================================================
  
  #maskCpG(element.msa, names, anc$seq)
  
  #ali <- har.msa
  #names=c("hg18","panTro2","rheMac2")
  
  species.count <- length(names)
  
  inds = sapply(names, function(x) which(ali$names == x))
  seqs = sapply(inds,  function(x) toupper(ali[[1]][x]))	 #- keep order, all UC
  SEQS = sapply(seqs,  function(x) strsplit(x,split=""))
  SEQS = matrix(unlist(SEQS), byrow=TRUE, nrow=species.count)
  
  #anc  = apply(SEQS,2,parsimony.func)	#recorre por columnas
  anc <- unlist(strsplit(anc.fa$seq, ""))
  
  len  = dim(SEQS)[2]	
  
  #- informative sites
  inf  = SEQS[1,] != SEQS[2,]
  
  #analizo el primero
  cls  = rep(NA, len)
  if( any(SEQS[,2]=="G")){
    if(anc[1] == "C" & anc[2] == "G") {
      cls[1] = 2
    } else {
      cls[1] = 3
    }
  } else {
    cls[1] = 1
  }
  
  #analizo del segundo la anteultimo
  for(i in 2:(len-1)){
    
    #- potential CpG
    if(any(SEQS[,i-1] == "C") | any(SEQS[,i+1]=="G")){
      if(anc[i] == "C" & anc[i+1] == "G"){
        cls[i] = 2
        next
      } 
      if(anc[i] == "G" & anc[i-1] == "C"){
        cls[i] = 2
        next
      }
      cls[i] = 3
      
    } else{
      cls[i] = 1
    }
  }
  
  #analizo el ultimo
  if( any(SEQS[,len-1]=="C")){
    if(anc[len-1] == "C" & anc[len] == "G") {
      cls[len] = 2
    } else {
      cls[len] = 3
    }
  } else {
    cls[len] = 1
  }
  
  msk.ind  = (cls != 1) & (inf)
  
  LSEQS    = sapply(ali[[1]],function(x)strsplit(x,split=""))
  #enmascaro:
  LMSEQS   = lapply(LSEQS,function(x) {x[msk.ind] = "N"; x }) 
  ali[[1]] = sapply(1:length(LSEQS),function(i) paste(LMSEQS[[i]],collapse=""))
  return(ali)
}



main <- function(names, exclusion, anc.file.suffix, alt.anc.file.suffix)
{
  files <- list.files(element.msa.base.path)
  
  for(element.msa.file.name in files)
  {
    if(! element.msa.file.name %in% exclusion)
    {
      anc.file.name <- gsub(".ss", anc.file.suffix, element.msa.file.name)
    }
    else
    {
      anc.file.name <- gsub(".ss", alt.anc.file.suffix, element.msa.file.name)
    }
      
    #"chr11_TSAR.1577_anc.hg19-monDom5.fa"
    
    orgs <- get.species.lst(ali.orgs.base.path, ali.orgs.file.name)
    
    setwd(element.msa.base.path)
    element.msa <- read.msa(element.msa.file.name,seqnames=orgs, format = "SS")
    
    setwd(anc.base.path)
    anc <- read.msa(anc.file.name,format = "FASTA")
    
    ali.msk <- maskCpG(element.msa, names, anc)
    
    setwd(msk.element.msa.base.path)
   
    write.msa(x = ali.msk, format = "SS", file = element.msa.file.name)
    
  } 
}

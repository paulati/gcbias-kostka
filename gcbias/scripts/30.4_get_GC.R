

get.GCfreqs <- function(msa.ss)
{
  freqs <- base.freq.msa(msa.ss)
  bind <- which(as.character(freqs[,"states"]) %in% c("A","C","G","T"))
  freqs <- freqs[bind,]
  freqs[,"freq"] <- freqs[,"freq"] / sum(freqs[,"freq"])
  cind <- which(as.character(freqs[,"states"]) == "C")
  gind <- which(as.character(freqs[,"states"]) == "G")
  gcfreq <- sum(freqs[c(cind,gind),"freq"])
  gcfreq <- round(gcfreq, 3)
  
  #free mem
  rm(cind,gind,freqs,bind)
  gc()
  
  return (gcfreq)
}


#PAULA TODO: ver que esto solo sirve para kostka, generalizarlo
get.anc.GCfreqs <- function(msa)
{
  
  #TODO: esto hay que rehacerlo para el analisis de anabella
  
  seqs <- rbind(	strsplit(msa["hg18"]$seq, split="")[[1]],
                 strsplit(msa["panTro2"]$seq, split="")[[1]],
                 strsplit(msa["rheMac2"]$seq, split="")[[1]])
  
  anc.seq <- apply(seqs, 2, parsimony.func)
  
  inds <- anc.seq %in% c("A","C","G","T")
  
  anc.gcfreq <- sum((table(anc.seq[inds]) / sum(table(anc.seq[inds])))[c("C","G")])
  
  anc.gcfreq <- round(anc.gcfreq, 3)
  
  return(anc.gcfreq)
}

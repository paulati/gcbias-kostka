element.models <- function(branch.name, init.mod, ali, SEL.NEUTRAL,smax=NULL,bmax=NULL){
  #======================================================================
  
  #print(SEL.NEUTRAL)
  #print(smax)
  
  #init.mod <- init.mod.tra
  #ali <- har.msa
  #SEL.NEUTRAL
  #smax=NULL
  #bmax=NULL
  
  M0.kostka = init.mod
  
  
  #PAULA nuevo
  #Ma = add.ls.mod(init.mod,branch="hg19",separate.params=paste("sel[",SEL.NEUTRAL,",", smax, "]",sep=""))
  #Ma.1 = add.ls.mod(init.mod,branch="hg19",separate.params=paste("sel[",0,",", smax, "]",sep=""), selection=SEL.NEUTRAL)
  #Ma.2 = add.ls.mod(init.mod,branch="hg19",separate.params=paste("sel[",0,",", smax, "]",sep=""), selection=SEL.NEUTRAL, bgc=0)
  #Ma.3 = add.ls.mod(init.mod,branch="hg19",separate.params=paste("sel[",0,",", smax, "]",sep=""), selection=0, bgc=0)
  Ma.4 = add.ls.mod(init.mod, branch=branch.name, separate.params=paste("sel[0,", smax, "]",sep=""), selection=SEL.NEUTRAL, bgc=0)
  #Ma.mail <- add.ls.mod(init.mod,branch=branch.name, separate.params=paste("sel[",SEL.NEUTRAL,",", smax, "]",sep=""), bgc = 0)
  #Ma.mail <- add.ls.mod(init.mod,branch=branch.name, separate.params=paste("sel[",0,",", smax, "]",sep=""), bgc = 0)
  
  
  
  Ma.rel.kostka <- add.ls.mod(init.mod,branch=branch.name,separate.params=paste("sel[0,",SEL.NEUTRAL,"]",sep=""))
  #Ma.rel.1 = add.ls.mod(init.mod,branch="hg19",separate.params=paste("sel[0,",SEL.NEUTRAL,"]",sep=""), bgc=0)
  #Ma.rel es identico a Ma.rel.1 y Ma.rel.4, Ma.rel.2 es distinto
  #Ma.rel.2 = add.ls.mod(init.mod,branch="hg19",separate.params=paste("sel[0,",SEL.NEUTRAL,"]",sep=""), selection=SEL.NEUTRAL, bgc=0)
  #Ma.rel.3 = add.ls.mod(init.mod,branch="hg19",separate.params=paste("sel[0,",SEL.NEUTRAL,"]",sep=""), selection=SEL.NEUTRAL, bgc=0)
  Ma.rel.4 = add.ls.mod(init.mod,branch=branch.name,separate.params=paste("sel[0,",SEL.NEUTRAL,"]",sep=""), bgc=0, selection=0)
  

  
  Ma.pos.kostka = add.ls.mod(init.mod,branch=branch.name,separate.params=paste("sel[",SEL.NEUTRAL,",",smax,"]",sep=""),selection=SEL.NEUTRAL)
  #Ma.pos.1 = add.ls.mod(init.mod,branch="hg19",separate.params=paste("sel[",SEL.NEUTRAL,",",smax,"]",sep=""),selection=SEL.NEUTRAL, bgc=0)
  # Ma.pos y Ma.pos.1 son iguales
  #Ma.pos.2 = add.ls.mod(init.mod,branch="hg19",separate.params=paste("sel[",SEL.NEUTRAL,",",smax,"]",sep=""))
  

  Mb.kostka = add.ls.mod(init.mod,branch=branch.name,separate.params=paste("bgc[0,",bmax,"]",sep=""))
  #Mb.1 = add.ls.mod(init.mod,branch="hg19",separate.params=paste("bgc[0,",bmax,"]",sep=""), selection=0)

  
  #original:
  Mab.kostka = add.ls.mod(init.mod,branch=branch.name,separate.params=paste("bgc[0,",bmax,"],sel[0,",smax,"]",sep=""))
  Mab.1 = add.ls.mod(init.mod,branch=branch.name,separate.params=paste("bgc[0,",bmax,"],sel[", SEL.NEUTRAL, ",", smax,"]", sep=""), selection=SEL.NEUTRAL)
  #Mab.2 = add.ls.mod(init.mod,branch="hg19",separate.params=paste("bgc[0,",bmax,"],sel[0,",smax,"]",sep=""), selection=SEL.NEUTRAL)
  #Mab y Mab.1 son iguales, Mab.2 es distinto

  
  result <- list(	null  = M0.kostka,
        sel   = Ma.4,
        sel.r = Ma.rel.4,
        sel.p = Ma.pos.kostka,
        bgc   = Mb.kostka,
        both  = Mab.kostka )
  
  return(result)
  
}








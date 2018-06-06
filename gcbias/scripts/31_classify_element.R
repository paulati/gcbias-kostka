
likes.from.fit <- function(fit){
  #-------------------------------
  #- fix fitting issues; needs to go in SGE code...	
  
  #fit <- RES$SIM$S[[1]]$fit[[1]]
  
  #likelihood de cada modelo - likelihood del modelo nulo:
  tmp <- fit[-1,1] - fit[1,1]
  
  if(any(is.na(tmp)))     {
    return(rep(NA,5))
  }
  if(tmp[1] < min(tmp[2],tmp[3])) { 
    tmp[1] <-  min(tmp[2],tmp[3]) 
  }
  #- throw out messed up fits
  if(tmp[1] < 0) { return(rep(NA,5))}
  if(tmp[4] < 0) { return(rep(NA,5))}
  if(tmp[5] < 0) { return(rep(NA,5))}
  
  
  return(tmp)
  
}


apfu_2 <- function(x, l.sel){
  #-------------------
  #- slightly different, need index here 
  #original:
  result <- quantile( l.sel[[x]]$likes[,3] -  l.sel[[x]]$likes[,2], .95, na.rm=TRUE)
  #z.H0 <- l.sel[[x]]$likes[,2] -  l.sel[[x]]$likes[,3]
  #result <- quantile(z.H0 , .95, na.rm=TRUE)
  return(result)
}	


#- extract likelihood scores for simulations under gBGC and selection
#--------------------------------------------------------------------
apfu <- function(x) {
  #--------------------
  likes = matrix(unlist(lapply(x$fit, likes.from.fit)),ncol=5,byrow=TRUE)
  colnames(likes) = c("s","s.relax","s.pos","bgc","s.bgc")
  return(list(likes=likes,val=x$val))
}

load.RES <- function(analysis.base.path, analysis.file.name)
{
  data.file.path <- paste(analysis.base.path, analysis.file.name, sep="/")
  print(data.file.path)
  
  tst <- try(load(data.file.path))
  if(tst != "RES") { 
    result <- NULL
  } else {
    result <- RES
  }
  
  return(result)
}


#PAULA TODO: ver si rescribo este codigo para abstraer que se analiza msa o msa.msk
# en los arreglos que mando kostka lo tiene mejor

#- Assign the class of the ELEMENT
#=================================
#getElementClass <- function(RES, analysis.base.path, analysis.file.name)
getElementClass <- function(RES)
{
  
  # data.file.path <- paste(analysis.base.path, analysis.file.name, sep="/")
  # print(data.file.path)
  # 
  # tst <- try(load(data.file.path))
  # if(tst != "RES") { return(NULL)}
  
  
  #- need to get classification boundaries
  #========================================
  
  #- extract likelihood scores for simulations under the NULL 
  #----------------------------------------------------------
  
  #PAULA revisar aca!!
  #PAULA revisar: este 1 me parece que no va y tengo que iterar entre todos los simulados
  #armar la matriz que ponga como filas cada uno de los valores de S (l.null.s) y de B (l.null.b) y asi construir el l.null
  #PAULA revisar aca!!
  #tmp <- lapply(RES$SIM$S[[1]]$fit, likes.from.fit)
  
  #original:
  #esto mira solo las simulaciones NSIM para el primer valor de s (s=0) 
  #simulado, en lugar de los 102 valores de s
  #l.null.s = matrix(unlist(lapply(RES$SIM$S[[1]]$fit, likes.from.fit)),ncol=5,byrow=TRUE)
  #colnames(l.null.s) = c("s","s.relax","s.pos","bgc","s.bgc")
  
  #original:
  #esto mira solo las simulaciones NSIM para el primer valor de b (b=0) 
  #l.null.b = matrix(unlist(lapply(RES$SIM$B[[1]]$fit, likes.from.fit)),ncol=5,byrow=TRUE)
  #colnames(l.null.b) = c("s","s.relax","s.pos","bgc","s.bgc")
  
  #reemplazo lo de arriba por esto para reusar la funcion que calcula los otros l
  x.s <- RES$SIM$S[[1]]
  tmp.s <- apfu(x.s)
  l.null.s <- tmp.s$likes
  
  x.b <- RES$SIM$B[[1]]
  tmp.b <- apfu(x.b)
  l.null.b <- tmp.b$likes
  
  l.null   = rbind(l.null.s, l.null.b)
  
  #l.null   = rbind(l.null.ss, l.null.bs)
  na.null  = sum(is.na(rowSums(l.null)))
  
  #tmp <- RES$SIM$S[[102]]$fit
  
  alpha = 0.05
  
  s.col.index <- 1
  bgc.col.index <- 4
  s.bgc.col.index <- 5
  
  da.0  = quantile(l.null[,s.col.index], 1-alpha, na.rm=TRUE) 
  db.0  = quantile(l.null[,bgc.col.index], 1-alpha, na.rm=TRUE)
  dab.0 = quantile(l.null[,s.bgc.col.index], 1-alpha, na.rm=TRUE)
  
  l.gbgc = lapply(RES$SIM$B[-1],apfu )
  l.sel  = lapply(RES$SIM$S[-1],apfu )
  na.s = sum(unlist(lapply(l.sel, function(x) sum(is.na(rowSums(x$likes))))))
  na.b = sum(unlist(lapply(l.gbgc,function(x) sum(is.na(rowSums(x$likes))))))
  
  da.bs  = unlist(lapply(l.gbgc, function(x) quantile(x$likes[,s.col.index]-x$likes[,bgc.col.index],1-alpha, na.rm=TRUE)))
  dab.bs = unlist(lapply(l.gbgc, function(x) quantile(x$likes[,s.bgc.col.index]-x$likes[,bgc.col.index],1-alpha, na.rm=TRUE)))
  
  #PAULA: ver que la funcion max esta devolviendo NA, agrego na.rm=TRUE :
  da.b  = max(c(da.bs,0), na.rm=TRUE)
  dab.b = max(c(dab.bs,0), na.rm=TRUE)
  
  db.as  = unlist(lapply(l.sel, function(x) quantile(x$likes[,bgc.col.index]-x$likes[,s.col.index],1-alpha, na.rm=TRUE)))
  dab.as = unlist(lapply(l.sel, function(x) quantile(x$likes[,s.bgc.col.index]-x$likes[,s.col.index],1-alpha, na.rm=TRUE))) 
  
  #PAULA: ver que la funcion max esta devolviendo NA, agrego na.rm=TRUE :
  db.a = max(c(db.as,0), na.rm = TRUE)
  dab.a = max(c(dab.as,0), na.rm = TRUE)
  
  likelihood.col.index <- 1
  m0.row.index <- 1
  ms.row.index <- 2
  msrelax.row.index <- 3
  mspos.row.index <- 4
  mbgc.row.index <- 5
  msbgc.row.index <- 6
  
  #l0    = RES$HAR$fit$summary[1,1]
  l0    = RES$ELEM$fit[m0.row.index,likelihood.col.index]
  #la    = RES$HAR$fit$summary[2,1]- l0
  la    = RES$ELEM$fit[ms.row.index,likelihood.col.index]- l0
  #la.r  = RES$HAR$fit[3,1]- l0
  la.r  = RES$ELEM$fit[msrelax.row.index,likelihood.col.index]- l0
  #la.p  = RES$HAR$fit$summary[4,1]- l0
  la.p  = RES$ELEM$fit[mspos.row.index,likelihood.col.index]- l0
  #lb    = RES$HAR$fit$summary[5,1]- l0
  lb    = RES$ELEM$fit[mbgc.row.index,likelihood.col.index]- l0
  #lab   = RES$HAR$fit$summary[6,1]- l0
  lab   = RES$ELEM$fit[msbgc.row.index,likelihood.col.index]- l0
  
  class = assign.class(la,lb,lab,db.0,da.0,dab.0,da.b,db.a,dab.b,dab.a)
  
  #l0.msk    = RES$HAR$fit.msk$summary[1,1]
  l0.msk    = RES$ELEM$fit.msk[m0.row.index,likelihood.col.index]
  #la.msk    = RES$HAR$fit.msk$summary[2,1]- l0.msk
  la.msk    = RES$ELEM$fit.msk[ms.row.index,likelihood.col.index]- l0.msk
  #la.r.msk  = RES$HAR$fit.msk$summary[3,1]- l0.msk
  la.r.msk  = RES$ELEM$fit.msk[msrelax.row.index,likelihood.col.index]- l0.msk
  #la.p.msk  = RES$HAR$fit.msk$summary[4,1]- l0.msk
  la.p.msk  = RES$ELEM$fit.msk[mspos.row.index,likelihood.col.index]- l0.msk
  #lb.msk    = RES$HAR$fit.msk$summary[5,1]- l0.msk
  lb.msk    = RES$ELEM$fit.msk[mbgc.row.index,likelihood.col.index]- l0.msk
  #lab.msk   = RES$HAR$fit.msk$summary[6,1]- l0.msk
  lab.msk   = RES$ELEM$fit.msk[msbgc.row.index,likelihood.col.index]- l0.msk
  
  class.msk = assign.class(la.msk,lb.msk,lab.msk,db.0,da.0,dab.0,da.b,db.a,dab.b,dab.a)
  
  
  
  #return(list(class=class, class.msk=class.msk))
  
  
  #esto distingue seleccion de relajacion:
  #- refine the selection class
  #============================
  
  #SG = abs(RES$HAR$fit$summary[1,3])
  SG = abs(RES$ELEM$fit[1,3])
  
  #sel.neu.inds = which(unlist(lapply(RES$SIM$S , function(x) x$val <= SG)))
  
  #asumo H0: Ls >= Lr, entonces (Ls - Lr >= 0)
  #en seleccion s > SG, H0.inds son los elementos que verifican H0
  #quito el primer elemento porque tambien lo quitan en el calculo de l.sel
  #H0.inds <- which(unlist(lapply(RES$SIM$S[-1], function(x) x$val >= SG)))
  #original:
  sel.neu.inds = which(unlist(lapply(RES$SIM$S , function(x) x$val <= SG)))
  
  #busco los valores criticos con significancia 0.05
  dsps = sapply(sel.neu.inds,apfu_2, l.sel)
  #dsps = sapply(H0.inds, apfu_2, l.sel)
  
  # less.zero.index <- which(as.numeric(dsps) < 0)
  # less.zero <- dsps[less.zero.index]
  # prob <- length(less.zero.index) / length(dsps)
  # #as.numeric(dsps)
  
  
  #PAULA : agrego na.rm
  dsp  = max(dsps, na.rm = TRUE)
  #dsp  = mean(dsps, na.rm = TRUE)
  #probar min max

  # if(class[1] ==1){
  #   if(la.r - la.p > dsp  ){
  #     class[1] = 10 #a-
  #   } else{
  #     class[1] = 11 #a+
  #   }
  # }
  
  if(class[1] ==1){
    if(la.p - la.r > dsp){
      class[1] = 11
    } else{
      class[1] = 10
    }
  }
  if(class.msk[1] ==1){
    if(la.p.msk - la.r.msk > dsp){
      class.msk[1] = 11
    } else{
      class.msk[1] = 10
    }
  }

  
  #- parameter estimates
  #=====================
  if(class[1]==0){
    #est=RES$HAR$fit$summary[1,]
    est=RES$ELEM$fit[1,]
  }
  if(class[1]==22 | class[1] == 20){
    #est=RES$HAR$fit$summary[5,]
    est=RES$ELEM$fit[5,]
  }
  if(class[1]==3){
    #est=RES$HAR$fit$summary[6,]
    est=RES$ELEM$fit[6,]	
  }
  if(class[1]==10){
    #est=RES$HAR$fit$summary[3,]
    est=RES$ELEM$fit[3,]
  }
  if(class[1]==11){
    #est=RES$HAR$fit$summary[4,]
    est=RES$ELEM$fit[4,]
  }
  
  if(class.msk[1]==0){
    #est.msk=RES$HAR$fit.msk$summary[1,]
    est.msk=RES$ELEM$fit.msk[1,]
  }
  if(class.msk[1]==22 | class.msk[1] == 20){
    #est.msk=RES$HAR$fit.msk$summary[5,]
    est.msk=RES$ELEM$fit.msk[5,]
  }
  if(class.msk[1]==3){
    #est.msk=RES$HAR$fit.msk$summary[6,]
    est.msk=RES$ELEM$fit.msk[6,]
  }
  if(class.msk[1]==10){
    #est.msk=RES$HAR$fit.msk$summary[3,]
    est.msk=RES$ELEM$fit.msk[3,]
  }
  if(class.msk[1]==11){
    #est.msk=RES$HAR$fit.msk$summary[4,]
    est.msk=RES$ELEM$fit.msk[4,]
  }
  
#  ali.orgs.base <- "~/Documentos/r/2017/lara/data/har.msa/chr16/"
#  setwd(ali.orgs.base)
# ali     = read.msa(alifile)
# ali.m   = RES$HAR$ali.msk
# sbs     = get.subs(ali)
# sbs.m   = get.subs(ali.m)
# anc     = getAnc(ali)
# hg      = toupper(ali[[1]][which(ali$names=="hg18")])
# hg      = strsplit(hg, split="")[[1]]
# ind1    = anc %in% c("A","C","G","T")
# ind2    = hg %in% c("A","C","G","T")
# gc.hg   = sum(table(hg[ind1])[ c("C","G")])/sum(table(hg[ind1]))
# gc.anc  = sum(table(anc[ind2])[ c("C","G")])/sum(table(anc[ind2]))
# dgc     = gc.hg - gc.anc
# dgc    	= gc.hg - gc.anc
# 
# anc       = getAnc(ali.m)
# hg        = toupper(ali.m[[1]][which(ali.m$names=="hg18")])
# hg        = strsplit(hg, split="")[[1]]
# gc.hg     = sum(table(hg[ind1])[ c("C","G")])/sum(table(hg[ind1]))
# gc.anc    = sum(table(anc[ind2])[ c("C","G")])/sum(table(anc[ind2]))
# dgc.m     = gc.hg - gc.anc
  
#  print(class)
#  cat("\n")
  
  #print("!!!!!")
  #print(RES$ELEM$fit[1,3])
  
  result <- list(class=class, class.msk=class.msk, 
                 db.0=db.0, da.0 = da.0,dab.0=dab.0,da.b = da.bs,db.a = db.as, dab.b = dab.bs, dab.a = dab.as,
                 #est=est,est.msk=est.msk,
                 #subs=sbs,subs.msk=sbs.m, 
                 #dgc.msk=dgc.m,dgc=dgc, 
                 rate.MN=RES$MOD$rate.MN, 
                 rate.M0 = RES$ELEM$fit[1,2], 
                 rate.M0.msk = RES$ELEM$fit.msk[1,2],
                 #SG=RES$ELEM$fit[1,3],SG.msk = RES$ELEM$fit.msk[1,3],
                 SG=RES$ELEM$fit[1,3], SG.msk = RES$ELEM$fit.msk[1,3],
                 nas=c(na.null,na.s,na.b),
                 fit=RES$ELEM$fit, fit.msk=RES$ELEM$fit.msk)
  
  return(result)  

}


assign.class <- function(la,lb,lab,db.0, da.0, dab.0, da.b, db.a, dab.b, dab.a){
  #-------------------------------------------------------------------------------
  
  #la <- la.msk
  #lb <- lb.msk
  #lab <- lab.msk
  
  
  cuts  = c(db.0, da.0, dab.0,  da.b, db.a , dab.b , dab.a )
  vals  = c(lb  ,   la,   lab, la-lb, lb-la, lab-lb, lab-la)
  
  #PAULA cocina: agrego esta condicion para que no me devuelva como distintos los que son muy cercanos (vals - cuts > 1e-5)
  #vals - cuts
  
  comps = (vals > cuts)
  names(comps) = c("b>0","a>0","ab>0","a>b","b>a","ab>b","ab>a")
  
  #- the easy cases, only one model rejects null
  #---------------------------------------------
  if(all(comps[1:3] == c(FALSE,FALSE,FALSE))) return(c(0,comps))   #- null
  if(all(comps[1:3] == c(TRUE ,FALSE,FALSE))) {
    if(comps[5] == TRUE  ) 		    return(c(22,comps))  #- bgc + high-conf
    if(comps[5] == FALSE )		    return(c(20,comps))  #- bgc
  }
  if(all(comps[1:3] == c(FALSE,TRUE ,FALSE))) return(c(1,comps))   #- sel
  if(all(comps[1:3] == c(FALSE,FALSE,TRUE ))) return(c(3,comps))   #- bgc & sel
  
  #- two models reject null
  #------------------------
  
  #- bgc and sel reject null
  if(all(comps[1:3] == c(TRUE ,TRUE ,FALSE))){
    if(all(comps[4:5] == c(TRUE ,FALSE))) return(c(1,comps))   #- sel
    if(all(comps[4:5] == c(FALSE,TRUE ))) return(c(22,comps))  #- bgc + high-conf
    if(all(comps[4:5] == c(FALSE,FALSE))) return(c(20,comps))  #- bgc
    #if(comps[4:5] == c(TRUE ,FALSE)) cannot be
  }
  
  #- sel and sel&bgc reject null
  if(all(comps[1:3] == c(FALSE ,TRUE ,TRUE ))){
    if(comps[7] == TRUE ) return(c(3,comps)) #- bgc & sel
    if(comps[7] == FALSE) return(c(1,comps)) #- sel
  }
  
  #- bgc and sel&bgc reject null
  if(all(comps[1:3] == c(TRUE ,FALSE,TRUE ))){
    if(comps[6] == TRUE ) return(c(3,comps)) #- bgc & sel
    if(comps[6] == FALSE) {
      if(comps[5] == TRUE  ) return(c(22,comps)) 
      if(comps[5] == FALSE ) return(c(20,comps))
    }
  }
  
  #- all three models reject null
  #------------------------------
  if(all(comps[1:3] == c(TRUE ,TRUE ,TRUE ))){
    #- sel rejects bgc
    if(comps[4] == TRUE){
      if(comps[7] == TRUE ) return(c(3,comps)) #- bgc & sel
      if(comps[7] == FALSE) return(c(1,comps)) #- sel
    }
    #- bgc rejects sel
    if(comps[5] == TRUE){
      if(comps[6] == TRUE ) return(c(3,comps))  #- bgc & sel
      if(comps[6] == FALSE) return(c(22,comps)) #- bgc + high-conf
    }
    #- no bgc vs sel decision
    if(all(comps[4:5] == c(FALSE,FALSE))){
      if(all(comps[6:7] == c(TRUE ,TRUE ))) return(c(3,comps))
      if(all(comps[6:7] == c(TRUE ,FALSE))) return(c(3,comps))  #- cons
      if(all(comps[6:7] == c(FALSE,TRUE ))) return(c(20,comps)) #- cons
      if(all(comps[6:7] == c(FALSE,FALSE))) return(c(20,comps)) #- cons
    }
  }
  
  return(NA)
}

# res.har <- tst$ELEM$fit
# res.har.msk <- tst$ELEM$fit.msk
# sim <- tst$SIM
# sel.neutral <- tst$MOD$snu
# SEL.NEUTRAL <- sel.neutral
# result.refined <- refine.selection.class(res.har, res.har.msk, sim, SEL.NEUTRAL)

#- refine the selection class
#============================
#refine.selection.class <- function(res.har, res.har.msk, sim, SEL.NEUTRAL)
# refine.selection.class <- function(res.har, res.har.msk, sim, sg)
# {
#   
#   # smax = ceiling(max(c(res.har[,3],res.har.msk[,3])))
#   # 
#   # #- 15 fold speedup (compared to neutral) minimum S_max
#   # if(smax < ceiling(SEL.NEUTRAL)*15){
#   #   smax = ceiling(SEL.NEUTRAL)*15
#   # }
#   
#   #ses = sort( c(seq(0,smax,len=101),SEL.NEUTRAL));
#   
#   ##################################################
#   #PAULA TODO: esto es medio raro, pero probe con SG <- sel.neutral y no mejora
#   
#   #SG = abs(res.har[1,3])  #original!
#   
#   # res.har[1,3] esta asignado a res$SG en el resultado de classify.elem, 
#   #cambie eso y ahora asigno a res$SG el valor de res.har[2,3].
#   # por lo tanto, está bien usar el valor de sg que vienen como parametro pero 
#   # hay que revisar como se asigna en el classify.elem
#   SG <- sg
#   ##################################################
#   
#   print(sg)
#   
#   #SG <- sel.neutral
#   #SG = abs(res.har[1,3])
#   #SG <- sg
#   
#   sel.neu.inds = which(unlist(lapply(sim$S , function(x) x$val <= SG)))
#   #sel.neu.inds = which(unlist(lapply(sim$S , function(x) x$val >= SG))) #doy vuelta la hipotesis, H0 es a+ y H1 es a-
#   
#   print(length(sel.neu.inds))
#   
#   # if(102 %in% sel.neu.inds)
#   # {
#   #   sel.neu.inds <- sel.neu.inds[1:length(sel.neu.inds)-1]
#   # }
#   
#   # x$val es el valor del s simulado
#   
#   l.sel  <- lapply(sim$S[-1], apfu)
#   
#   dsps <- sapply(sel.neu.inds, apfu_2, l.sel)
#   
#   dsp  <- max(dsps)
#   
#   # es lo mismo restar res.har[1,1] que no hacerlo porque en la.p - la.r se anula
#   la.p <- res.har[4,1] - res.har[1,1]
#   la.r <- res.har[3,1] - res.har[1,1]
#   
#   
#   if(la.p - la.r > dsp){
#   #if(la.r - la.p > dsp){
#     result = "rel"
#   } else{
#     result = "pos"
#   }
#   
#   la.p.msk <- res.har.msk[4,1]
#   la.r.msk <- res.har.msk[3,1]
#   
#   if(la.p.msk - la.r.msk > dsp){
#     result.msk = "pos"
#   }  else  {
#     result.msk = "rel"
#   }
#   
#   result <- list(refine=result, refine.msk=result.msk)
#   
#   return(result) 
#   
# }


decode.class <- function(bits.class){
  
  #col.names <- c("b>0", "a>0", "ab>0", "a>b", "b>a", "ab>b", "ab>a" )
  
  if(paste(bits.class[c("b>0", "a>0", "ab>0")], collapse = "") == "000"){
    result <- "0"
  } else if(paste(bits.class[c("b>0", "a>0", "ab>0")], collapse="") == "100" & bits.class["b>a"] == 1) {
    result <- "b+"
  } else if(paste(bits.class[c("b>0", "a>0", "ab>0")], collapse="") == "100" & bits.class["b>a"] == 0) {
    result <- "b-"
  } else if(paste(bits.class[c("b>0", "a>0", "ab>0")], collapse="") == "010"){
    result <- "a"
  }  else if(paste(bits.class[c("b>0", "a>0", "ab>0")], collapse="") == "001"){
    result <- "ab"
  } else if(paste(bits.class[c("b>0", "a>0", "ab>0", "a>b", "b>a")], collapse="") == "11010"){
    result <- "a"
  } else if(paste(bits.class[c("b>0", "a>0", "ab>0", "a>b", "b>a")], collapse="") == "11001"){
    result <- "b+"
  } else if(paste(bits.class[c("b>0", "a>0", "ab>0", "a>b", "b>a")], collapse="") == "11000"){
    result <- "b-"
  } else if(paste(bits.class[c("b>0", "a>0", "ab>0", "ab>a")], collapse="") == "0111"){
    result <- "ab"
  } else if(paste(bits.class[c("b>0", "a>0", "ab>0", "ab>a")], collapse="") == "0110"){
    result <- "a"
  } else if(paste(bits.class[c("b>0", "a>0", "ab>0", "ab>b")], collapse="") == "1011"){
    result <- "ab"
  } else if(paste(bits.class[c("b>0", "a>0", "ab>0", "b>a", "ab>b")], collapse="") == "10110"){
    result <- "b+"
  } else if(paste(bits.class[c("b>0", "a>0", "ab>0", "b>a", "ab>b")], collapse="") == "10100"){
    result <- "b-"
  } else if(paste(bits.class[c("b>0", "a>0", "ab>0", "a>b", "b>a", "ab>a")], collapse="") == "111101"){
    result <- "ab"
  } else if(paste(bits.class[c("b>0", "a>0", "ab>0", "a>b", "b>a", "ab>a")], collapse="") == "111100"){
    result <- "a"
  } else if(paste(bits.class[c("b>0", "a>0", "ab>0", "a>b", "b>a", "ab>b")], collapse="") == "111011"){
    result <- "ab"
  } else if(paste(bits.class[c("b>0", "a>0", "ab>0", "a>b", "b>a", "ab>b")], collapse="") == "111010"){
    result <- "b+"
  } else if(paste(bits.class[c("b>0", "a>0", "ab>0", "a>b", "b>a", "ab>b", "ab>a")], collapse="") == "1110011"){
    result <- "ab"
  } else if(paste(bits.class[c("b>0", "a>0", "ab>0", "a>b", "b>a", "ab>b", "ab>a")], collapse="") == "1110010"){
    result <- "ab"
  } else if(paste(bits.class[c("b>0", "a>0", "ab>0", "a>b", "b>a", "ab>b", "ab>a")], collapse="") == "1110001"){
    result <- "b-"
  } else if(paste(bits.class[c("b>0", "a>0", "ab>0", "a>b", "b>a", "ab>b", "ab>a")], collapse="") == "1110000"){
    result <- "b-"
  }
  
  if(bits.class[1] == 11)
  {
    result <- "a+"
  } else if(bits.class[1] == 10)
  {
    result <- "a-"
  }
  
  return(result)
  
}
  


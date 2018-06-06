
do.simulations <- function(element.msa, element.fit, element.fit.msk, sel.neutral, init.mod.tra, branch.name)
{
  #- the masking pattern for the HAR
  element.map <- make.map(element.msa)
  
  #- max encountered estimates for S and B
  shat.col.index <- 3
  bhat.col.index <- 4
  smax <- ceiling(max(c(element.fit[,shat.col.index],element.fit.msk[,shat.col.index])))
  bmax <- ceiling(max(c(element.fit[,bhat.col.index],element.fit.msk[,bhat.col.index])))
  
  #- 15 fold speedup (compared to neutral) minimum S_max
  fold.speedup <- 15 #15 5 1
  if(smax < ceiling(sel.neutral)*fold.speedup){
    smax <- ceiling(sel.neutral)*fold.speedup
  }
  if(bmax < round(ceiling(sel.neutral)*fold.speedup*1.5)){
    bmax <- round(ceiling(sel.neutral)*fold.speedup*1.5)
  }
  
  smin <- -smax
  #ses <- sort( c(seq(0,smax,len=101),sel.neutral));
  ses <- sort( c(seq(smin,smax,len=101),sel.neutral));
  
  #  ses <- sort( c(seq(0,smax,len=101),sel.neutral));
  bes <- seq(0,bmax,len=101)
  
  SIM <- list(S = vector(len=102,mode="list"), B = vector(len=101,mode="list"))
  
  
  for(i in 1:102){
    
    print(ses[i])
    
    #- transform generating model
    #============================
    cur.s = ses[i]
    #init.mod.tra.sim = add.alt.mod(init.mod.tra,"hg18",selection=cur.s)
    init.mod.tra.sim = add.ls.mod(init.mod.tra,branch=branch.name,selection=cur.s)
    SIM$S[[i]]$mod = init.mod.tra.sim
    SIM$S[[i]]$val = cur.s
    msas = vector(len=NSIM,mode="list")
    rslt = vector(len=NSIM,mode="list")
    for(j in 1:NSIM){
      # msas[[j]] <- simulate.msa(init.mod.tra.sim, nsim=dim(har.msa)[2])
      msas[[j]] <- simulate.msa(init.mod.tra.sim, nsim=dim(element.msa)[2] ) 
      # este nro 2 no tiene nada que ver con el numero NSIM, tiene que ver con indexar el elem que tiene la cantidad de bases en el alineamiento
      # dim devuelve un array de dos elementos donde el primero es la cantidad de especies y el segundo la cantidad de bases en el alineamiento
      msas[[j]] <- gappify(msas[[j]],element.map)
      rslt[[j]] <- fit.models(branch.name, init.mod.tra,msas[[j]],sel.neutral,smax=smax)
    }
    SIM$S[[i]]$dat = msas
    SIM$S[[i]]$fit = rslt
    
    #clean
    rm(msas)
    rm(rslt)
    gc() 
  }
  
  for(i in 1:101){
    
    print(bes[i])
    
    #- transform generating model
    #============================
    cur.b = bes[i]
    #init.mod.tra.sim = add.alt.mod(init.mod.tra,"hg18",bgc=cur.b)
    init.mod.tra.sim = add.ls.mod(init.mod.tra, branch=branch.name, bgc=cur.b)
    SIM$B[[i]]$mod = init.mod.tra.sim
    SIM$B[[i]]$val = cur.b
    msas = vector(len=NSIM,mode="list")
    rslt = vector(len=NSIM,mode="list")
    for(j in 1:NSIM){
      #msas[[j]] <- simulate.msa(init.mod.tra.sim, nsim=dim(har.msa)[2])
      msas[[j]] <- simulate.msa(init.mod.tra.sim, nsim=dim(element.msa)[2])
      # este nro 2 no tiene nada que ver con el numero NSIM, tiene que ver con indexar el elem que tiene la cantidad de bases en el alineamiento
      # dim devuelve un array de dos elementos donde el primero es la cantidad de especies y el segundo la cantidad de bases en el alineamiento
      msas[[j]] <- gappify(msas[[j]],element.map)
      rslt[[j]] <- fit.models(branch.name, init.mod.tra,msas[[j]],sel.neutral,bmax=bmax)
    }
    SIM$B[[i]]$dat = msas
    SIM$B[[i]]$fit = rslt
    
    #clean
    rm(msas)
    rm(rslt)
    gc() 
  }
  
  
  
  return(SIM)
  
}

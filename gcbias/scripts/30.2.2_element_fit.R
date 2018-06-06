element.fit <- function(mods, 
                        branch.name, 
                        ali, SEL.NEUTRAL,smax=NULL,bmax=NULL){
  #======================================================================
  

  res.null                =  phyloFit(    ali,
                                          init.mod=mods$null,
                                          scale.only=TRUE,
                                          no.opt=c("backgd", "ratematrix"),
                                          quiet=TRUE)

  res.s             = phyloFit(     ali,
                                    #init.mod=Ma.mail,  #use Ma.4 en todos los tests
                                    #init.mod=Ma.mail,  #use Ma.4 en todos los tests
                                    init.mod=mods$sel,
                                    scale.only=TRUE,
                                    no.opt=c("backgd", "ratematrix"),
                                    quiet=TRUE)
  
  res.s.relax             = phyloFit(     ali,
                                          init.mod=mods$sel.r,
                                          #init.mod=add.alt.mod(init.mod,"hg19",sep=paste("sel[0,",SEL.NEUTRAL,"]",sep="")),
                                          scale.only=TRUE,
                                          no.opt=c("backgd", "ratematrix"),
                                          quiet=TRUE)
  
  res.s.pos               = phyloFit(     ali,
                                          init.mod=mods$sel.p,
                                          #init.mod=add.alt.mod(init.mod,"hg19",sep=paste("sel[",SEL.NEUTRAL,",",smax,"]",sep=""),selection = SEL.NEUTRAL),
                                          scale.only=TRUE,
                                          no.opt=c("backgd", "ratematrix"),
                                          quiet=TRUE)
  
  
  res.b                   = phyloFit(     ali,
                                          init.mod=mods$bgc,
                                          #init.mod=add.alt.mod(init.mod,"hg19",sep=paste("bgc[0,",bmax,"]",sep="")),
                                          scale.only=TRUE,
                                          no.opt=c("backgd", "ratematrix"),
                                          quiet=TRUE)
  
  res.sb                  = phyloFit(     ali,
                                          init.mod=mods$both,
                                          #init.mod=add.alt.mod(init.mod,"hg19",sep=paste("bgc[0,",bmax,"],sel[0,",smax,"]",sep="")),
                                          scale.only=TRUE,
                                          no.opt=c("backgd", "ratematrix"),
                                          quiet=TRUE)
  
  res = matrix(0,nrow=6,ncol=4)
  
  colnames(res) = c(      "likelihood"    ,
                          "hg-branchlen"  ,
                          "shat"          ,
                          "bhat"          )
  
  rownames(res) = c(      "null"          ,
                          "s"             ,
                          "s.relax"       ,
                          "s.pos"         ,
                          "bgc"           ,
                          "s.bgc"         )
  
  res[1,1] = res.null$likelihood
  res[1,2] = summary.tree(res.null$tree)[ which(summary.tree(res.null$tree)[,"name"]==branch.name)  ,"tparent"]
  #res[1,3] = res.null$ls.model$selection
  #res[1,4] = res.null$ls.model$bgc
  
  #PAULA aca falta La:
  res[2,1] = res.s$likelihood
  tmp      = summary.tree(res.s$tree)[ which(summary.tree(res.s$tree)[,"name"]==branch.name)  ,"tparent"]
  tmp      = sum(diag(res.s$ls.model$rate.matrix)  * res.s$backgd)*tmp*(-1)
  res[2,2] = tmp
  res[2,3] = res.s$ls.model$selection        
  
  #esta es la fila 3 La.r
  res[3,1] = res.s.relax$likelihood
  tmp      = summary.tree(res.s.relax$tree)[ which(summary.tree(res.s.relax$tree)[,"name"]==branch.name)  ,"tparent"]
  tmp      = sum(diag(res.s.relax$ls.model$rate.matrix)  * res.s.relax$backgd)*tmp*(-1)
  res[3,2] = tmp
  #res[2,3] = res.s.relax$ls.model$selection
  res[3,3] = res.s.relax$ls.model$selection
  
  #esta es la fila 4 La.p
  res[4,1] = res.s.pos$likelihood
  tmp      = summary.tree(res.s.pos$tree)[ which(summary.tree(res.s.pos$tree)[,"name"]==branch.name)  ,"tparent"]
  #VER
  tmp      = sum(diag(res.s.pos$ls.model$rate.matrix)  * res.s.pos$backgd)*tmp*(-1)
  res[4,2] = tmp
  #VER
  res[4,3] = res.s.pos$ls.model$selection
  
  #esta es la fila 5 Lb
  res[5,1] = res.b$likelihood
  tmp      = summary.tree(res.b$tree)[ which(summary.tree(res.b$tree)[,"name"]==branch.name)  ,"tparent"]
  #VER
  tmp      = sum(diag(res.b$ls.model$rate.matrix)  * res.b$backgd)*tmp*(-1)
  res[5,2] = tmp
  #VER
  res[5,4] = res.b$ls.model$bgc
  
  #esta es la fila 6 Lab
  res[6,1] = res.sb$likelihood
  tmp      = summary.tree(res.sb$tree)[ which(summary.tree(res.sb$tree)[,"name"]==branch.name)  ,"tparent"]
  tmp      = sum(diag(res.sb$ls.model$rate.matrix)  * res.sb$backgd)*tmp*(-1)
  res[6,2] = tmp
  res[6,3] = res.sb$ls.model$selection
  res[6,4] = res.sb$ls.model$bgc
  
  return(res)
  
}








fit.models <- function(branch.name, init.mod, ali, SEL.NEUTRAL,smax=NULL,bmax=NULL){
#======================================================================
  
  #print(SEL.NEUTRAL)
  #print(smax)
  
      #init.mod <- init.mod.tra
      #ali <- har.msa
      #SEL.NEUTRAL
      #smax=NULL
      #bmax=NULL
  
        M0.kostka = init.mod
        res.null                =  phyloFit(    ali,
                                                init.mod=M0.kostka,
                                                scale.only=TRUE,
                                                no.opt=c("backgd", "ratematrix"),
						quiet=TRUE)

        #PAULA nuevo
        #Ma = add.ls.mod(init.mod,branch="hg19",separate.params=paste("sel[",SEL.NEUTRAL,",", smax, "]",sep=""))
        #Ma.1 = add.ls.mod(init.mod,branch="hg19",separate.params=paste("sel[",0,",", smax, "]",sep=""), selection=SEL.NEUTRAL)
        #Ma.2 = add.ls.mod(init.mod,branch="hg19",separate.params=paste("sel[",0,",", smax, "]",sep=""), selection=SEL.NEUTRAL, bgc=0)
        #Ma.3 = add.ls.mod(init.mod,branch="hg19",separate.params=paste("sel[",0,",", smax, "]",sep=""), selection=0, bgc=0)
        Ma.4 = add.ls.mod(init.mod, branch=branch.name, separate.params=paste("sel[0,", smax, "]",sep=""), selection=SEL.NEUTRAL, bgc=0)
        #Ma.mail <- add.ls.mod(init.mod,branch=branch.name, separate.params=paste("sel[",SEL.NEUTRAL,",", smax, "]",sep=""), bgc = 0)
        #Ma.mail <- add.ls.mod(init.mod,branch=branch.name, separate.params=paste("sel[",0,",", smax, "]",sep=""), bgc = 0)

        #Ma.1 y Ma.2 son iguales
        #Ma, Ma.3, Ma.4 son iguales
        res.s             = phyloFit(     ali,
                                          #init.mod=Ma.mail,  #use Ma.4 en todos los tests
                                          #init.mod=Ma.mail,  #use Ma.4 en todos los tests
                                          init.mod=Ma.4,
                                          scale.only=TRUE,
                                          no.opt=c("backgd", "ratematrix"),
                                          quiet=TRUE)
        
        Ma.rel.kostka <- add.ls.mod(init.mod,branch=branch.name,separate.params=paste("sel[0,",SEL.NEUTRAL,"]",sep=""))
        #Ma.rel.1 = add.ls.mod(init.mod,branch="hg19",separate.params=paste("sel[0,",SEL.NEUTRAL,"]",sep=""), bgc=0)
        #Ma.rel es identico a Ma.rel.1 y Ma.rel.4, Ma.rel.2 es distinto
        #Ma.rel.2 = add.ls.mod(init.mod,branch="hg19",separate.params=paste("sel[0,",SEL.NEUTRAL,"]",sep=""), selection=SEL.NEUTRAL, bgc=0)
        #Ma.rel.3 = add.ls.mod(init.mod,branch="hg19",separate.params=paste("sel[0,",SEL.NEUTRAL,"]",sep=""), selection=SEL.NEUTRAL, bgc=0)
        Ma.rel.4 = add.ls.mod(init.mod,branch=branch.name,separate.params=paste("sel[0,",SEL.NEUTRAL,"]",sep=""), bgc=0, selection=0)
        
        res.s.relax             = phyloFit(     ali,
                                                init.mod=Ma.rel.kostka,
                                                #init.mod=add.alt.mod(init.mod,"hg19",sep=paste("sel[0,",SEL.NEUTRAL,"]",sep="")),
                                                scale.only=TRUE,
                                                no.opt=c("backgd", "ratematrix"),
						                                    quiet=TRUE)

        
        Ma.pos.kostka = add.ls.mod(init.mod,branch=branch.name,separate.params=paste("sel[",SEL.NEUTRAL,",",smax,"]",sep=""),selection=SEL.NEUTRAL)
        #Ma.pos.1 = add.ls.mod(init.mod,branch="hg19",separate.params=paste("sel[",SEL.NEUTRAL,",",smax,"]",sep=""),selection=SEL.NEUTRAL, bgc=0)
        # Ma.pos y Ma.pos.1 son iguales
        #Ma.pos.2 = add.ls.mod(init.mod,branch="hg19",separate.params=paste("sel[",SEL.NEUTRAL,",",smax,"]",sep=""))
        
        # Ma.pos.2 es el original del mail
        res.s.pos               = phyloFit(     ali,
                                                init.mod=Ma.pos.kostka,
                                                #init.mod=add.alt.mod(init.mod,"hg19",sep=paste("sel[",SEL.NEUTRAL,",",smax,"]",sep=""),selection = SEL.NEUTRAL),
                                                scale.only=TRUE,
                                                no.opt=c("backgd", "ratematrix"),
						                                    quiet=TRUE)

        
        Mb.kostka = add.ls.mod(init.mod,branch=branch.name,separate.params=paste("bgc[0,",bmax,"]",sep=""))
        #Mb.1 = add.ls.mod(init.mod,branch="hg19",separate.params=paste("bgc[0,",bmax,"]",sep=""), selection=0)
        res.b                   = phyloFit(     ali,
                                                init.mod=Mb.kostka,
                                                #init.mod=add.alt.mod(init.mod,"hg19",sep=paste("bgc[0,",bmax,"]",sep="")),
                                                scale.only=TRUE,
                                                no.opt=c("backgd", "ratematrix"),
						quiet=TRUE)

        #original:
        Mab.kostka = add.ls.mod(init.mod,branch=branch.name,separate.params=paste("bgc[0,",bmax,"],sel[0,",smax,"]",sep=""))
        Mab.1 = add.ls.mod(init.mod,branch=branch.name,separate.params=paste("bgc[0,",bmax,"],sel[", SEL.NEUTRAL, ",", smax,"]", sep=""), selection=SEL.NEUTRAL)
        #Mab.2 = add.ls.mod(init.mod,branch="hg19",separate.params=paste("bgc[0,",bmax,"],sel[0,",smax,"]",sep=""), selection=SEL.NEUTRAL)
        #Mab y Mab.1 son iguales, Mab.2 es distinto
        
        res.sb                  = phyloFit(     ali,
                                                init.mod=Mab.kostka,
                                               #init.mod=add.alt.mod(init.mod,"hg19",sep=paste("bgc[0,",bmax,"],sel[0,",smax,"]",sep="")),
                                                scale.only=TRUE,
                                                no.opt=c("backgd", "ratematrix"),
						quiet=TRUE)

        #modif paula
        #res.sb                  = phyloFit(     ali,
        #                                        init.mod=add.ls.mod(init.mod,branch="hg19",separate.params=paste("bgc[0,",bmax,"],sel[",SEL.NEUTRAL,",",smax,"]",sep="")),
                                                #init.mod=add.alt.mod(init.mod,"hg19",sep=paste("bgc[0,",bmax,"],sel[0,",smax,"]",sep="")),
        #                                        scale.only=TRUE,
        #                                        no.opt=c("backgd", "ratematrix"),
        #                                        quiet=TRUE)
        
                
        #paula ver esto: la condicion para Lab es S >= Sg and B >= 0 y eso se escribiria como sel[SEL.NEUTRAL, smax]       
        #PROBARLO CON EL FIX

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








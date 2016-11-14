## BUG: que pasa cuando no se usa leaves???!!!
## Revisar como resolvi lo de go ids, leaves e info de clustering

library(igraph)
igraph.offset<-0
#igraph.offset<-1


require(hwriter)
require(Rgraphviz)
require(annotate)
require(GOstats)
require(csbl.go)
#set.prob.table(organism=9606,type="similarity")
# 9606 (for Homo sapiens), 
# 4932 (Saccharomyces cerevisiae), 
# 6239 (Caenorhabditis elegans), 
# 7227 (Drosophila melanogaster), 
# 10090 (Mus musculus), 
# 10116 (Rattus norvegicus)

 varnorm<-function(x){
  (x-mean(x))/sd(x)
 }


#funcion para comparar dos particiones
ComparaPart<-function(ppp1,ppp2){
 pp1<-sort(ppp1)
 pp2<-sort(ppp2)
 gp1<-names(pp1)
 gp2<-names(pp2)
 
 g.common<-intersect(gp1,gp2)
 p1<-pp1[g.common]
 p2<-pp2[g.common]
 
 p1=pp1
 p2=pp2
 clust1<-(unique(p1))
 clust2<-(unique(p2))
 m<-matrix(0,nrow=length(clust1)+1,ncol=length(clust2)+1)
 colnames(m)<-c(unique(p2),"nc")
 rownames(m)<-c(unique(p1),"nc")
 for(i in clust1){
  g<-names(p1[p1==i])
  iin<-g%in%names(p2)
  if(sum(iin)>0){
   tt<-table(p2[g[iin]])
   m[as.character(i),names(tt)]=tt
  }
  if(any(!iin)){
   m[as.character(i),"nc"]=sum(!iin)
  }
 }
 for(i in clust2){
  g<-names(p2[p2==i])
  iin<-g%in%names(p1)
  if(any(!iin)){
   m["nc",as.character(i)]=sum(!iin)
  }
 
 }
 io<-order(as.numeric(colnames(m[,-ncol(m)])))
 return(cbind(m[,io],nc=m[,ncol(m)]))
 #ior<-order(as.numeric(rownames(m[-nrow(m),])))
 #return(m[ior,io])
} 

## Detecto categorias:
##     significativas + con  uMin < |mapeos| < uMax + con mas de minHits hits de la lista
## lbioTheme
## mbioTheme
#set.prob.table(organism=9606,type="similarity")
# 9606 (for Homo sapiens), 
# 4932 (Saccharomyces cerevisiae), 
# 6239 (Caenorhabditis elegans), 
# 7227 (Drosophila melanogaster), 
# 10090 (Mus musculus), 
# 10116 (Rattus norvegicus)




## Funcion para generar lbioTheme: lista con info sobre
## leaves, mica y clustersGO de los conceptos GO sobreexpresados
## en hg
## input
##   hg
##   uMin   <-5            #tamanio minimo de una categoria a testear
##   uMax   <-300          #tamanio maximo de una categoria a testear
##   uSmallSize<-15        #tamanio considerado chico
##   minHits<-3            #nro minimo de hits en una categorias (de tamanio mayor uSmallSize) sobreexpresada para tenerla en cuenta
##   k.minHits<-TRUE       #para tamanios chicos (uMin < s < uSmallSize) relajo condicion sobre nro de hits
##                         #usando un factor que va de 0.5 para uMin, hasta minHits/uSmallSize para uSmallSize
##   pv         <-0.01     #pvalue de Fisher
##   pv.adjusted<-0.1      #no se usa
##   adjust.method="fdr"   #no se usa
##   
##   bplot                 #plotimage nro de categorias de cada tipo 
##
##   onlyLeaves           <-TRUE #individualizar y clusterear solo leaves
##   bcutree              <-TRUE #usa cutree at 1/(1+micaBranchMinIC) level, si FALSE usa cutDynamicTree
##   minClusterSize       <-2 #solo si bcutree es false 
##   micaBranchMinIC      <-3 #todo lo que esta por arriba de la dist 1/(1+micaBranchMinIC) es clustereado separado
##   minGOnodes2clusterize<-2 #para mas de minGOnodes2clusterize uso el dendrograma para mergear conceptos
##   gometric             <-"Resnik"
##
## output
##   lbioTheme
biocat.bioThemes<-function(hg=hg,uMin=5,uMax=300,uSmallSize=15,bplot=FALSE,verbose=FALSE,
                           minHits=3,k.minHits=TRUE,pv=0.01,organism="At",
                           onlyLeaves=TRUE,bcutree=TRUE,minClusterSize=2,
                           micaBranchMinIC=3,minGOnodes2clusterize=2,gometric="Resnik"){
                           
  org<-c("Hs","Sc","Dm","Mm","At")
  orgCode <-c(9606,4932,7227,10090,3702)
  names(orgCode)<-org

  if(!organism%in%org){
    warning(paste("Organism",organism,"not in {",paste(org,collapse=","),"} recognized set.\n"))
  }

  set.prob.table(organism=orgCode[organism],type="similarity")

  if(verbose){
   paramdf<-data.frame(var=cbind(
        c("uMin","uMax","uSmallSize","minHits","k.minHits","pv","onlyLeaves","bcutree",
        "micaBranchMinIC","minGOnodes2clusterize","gometric"),
        value=c(uMin,uMax,uSmallSize,minHits,k.minHits,pv,onlyLeaves,bcutree,micaBranchMinIC,minGOnodes2clusterize,gometric),
        description=c("tamanio minimo de una categoria a testear",
          "tamanio maximo de una categoria a testear",
          "tamanio considerado chico",
          "nro minimo de hits en una categorias (de tamanio mayor uSmallSize) sobreexpresada para tenerla en cuenta",
          "para tamanios chicos (uMin < s < uSmallSize) relajo condicion sobre nro de hits",
          "pvalue de Fisher",
          "al exportar a html, reportar solo leaves",
          "usa cutree at 1/(1+micaBranchMinIC) level, si FALSE usa cutDynamicTree",
          "todo lo que esta por arriba de la dist 1/(1+micaBranchMinIC) es clustereado separado",
          "para mas de minGOnodes2clusterize uso el dendrograma para mergear conceptos",
          "metrica para ordenar categorias GO")))
   #write.table(paramdf,file="params.csv",sep="\t",col.name=TRUE,row.names=FALSE,quote=FALSE)
   print(paramdf)
  }


  hh=1/(1+micaBranchMinIC)  
  
  lbioTheme<-list() 
  laux<-lapply(hg,function(a){     
      lres<-list()
      #if(length(a$gset)<
      for(ccat in c("BP","MF","CC","KEGG")){ 
       onto<-ccat
       if(onto=="BP") envGOANC=GOBPANCESTOR
       if(onto=="MF") envGOANC=GOMFANCESTOR
       if(onto=="CC") envGOANC=GOCCANCESTOR
       if(is.null(a[[ccat]])) next
       if(is.na(a[[ccat]])) next
       ucount<-universeCounts(a[[ccat]])
       vMinHits<-rep(minHits,length(ucount))
       #para categorias chicas
       #cambio linealmente el criterio (factor del tamanio del nodo) para nro de hits de acuerdo al tamanio de las categorias
       #(e.g. para la uMin pido que 0.5*uMin)
       if(k.minHits){
          ismall<-which(ucount<uSmallSize)
          xx<-ucount[ismall]
          aux<-floor((((minHits/uSmallSize-0.5)*(xx-uMin)/(uSmallSize-uMin)+0.5))*xx)
          vMinHits[ismall]<-aux
       }
       igo<-which(ucount>=uMin & ucount<uMax & geneCounts(a[[ccat]])>=vMinHits)
       gonames<-intersect(sigCategories(a[[ccat]],pv),names(igo))
       #if(!is.null(drop.evidence)) gonames
       #i<-which(p.adjust(pvalues(a[[ccat]]),method=adjust.method)<pv.adjusted)
       
       #ya tengo las categorias significativas
       #ahora me quedo con las leaves
       isAncestor<-NA
       if(onlyLeaves & ccat!="KEGG"){
        goleaves<-c()
        goAncestors<-mget(gonames,envGOANC)
        lancestors<-lapply(goAncestors,function(x){
                         aux<-gonames%in%x
                         return(gonames[aux])
                         })
        isAncestor<-names(reverseSplit(lancestors))
        isLeave<-gonames[!(gonames%in%isAncestor)]
        #gonames<-gonames[!(gonames%in%isAncestor)]
       #}
       #if(ccat!="KEGG"){
        #lres[[ccat]]<-list(ids=gonames,isAncestor=isAncestor)
        lres[[ccat]]<-list(ids=gonames,leaves=isLeave)
       }else{
        lres[[ccat]]<-list(ids=gonames)
       }
      }
      return(lres)})
  names(laux)<-seq_along(laux)     
  lbioTheme[["clusters"]]<-laux
  lbioTheme[["params"]]<-c(uMin=uMin,uMax=uMax,minHits=minHits,
                           pv=pv,onlyLeaves=onlyLeaves)      

  mbioTheme<-matrix(unlist(lapply(lbioTheme[["clusters"]][seq_along(hg)],function(x){
                  return(unlist(lapply(x[1:4],function(y){length(y[[1]])})))
             })),ncol=4,byrow=TRUE)
  colnames(mbioTheme)<-c("BP","MF","CC","KEGG")        
  
  for(i in 1:ncol(mbioTheme)){
    cat(colnames(mbioTheme)[i]," fraction non-annotated clusters:",
       signif(sum(mbioTheme[,i]==0)/length(mbioTheme[,i]),2),"\n")
  }
  cat("Fraction of clusters non-annotated to any bio cat",
      signif(1-sum(apply(mbioTheme,1,any)/length(mbioTheme[,1])),2),"\n")  
  
  if(bplot){
   image(t(mbioTheme),zlim=c(1,max(mbioTheme)),axes=FALSE,xlab="BioCategory",ylab="Cluster")
   axis(1,c(0,.33,.66,1),c("BP","MF","CC","KEGG"))
  }
 if(FALSE){
  unlist(lapply(hg,function(x){length(x$gset)}))
  
  
 }

 cat("Recognizing GO branches...\n")

 ## Clusterizo categorias GO
 ## y guardo en lbioTheme goids, leavesId clusterizados en grupos "homogeneos"
 aa<-lapply(laux,function(x){
 
# for(jj in seq_along(laux)){
#  x<-laux[[jj]]
  
  if(length(x)==0) return(x)
  if(names(x)[1]=="uMin") return(x)
  
  #if(names(x)[1]=="uMin") { warning(paste(jj))}
  for(ccat in c("BP","MF","CC")){
    #cat(paste(jj,ccat,"\n"))
    if(is.null(x[[ccat]])) next
    a<-x[[ccat]]
    if(onlyLeaves){
     gonames<-a$leaves
    }else{
     gonames<-a$ids
    }
    if(length(gonames)==0)next;
    if(length(gonames)<minGOnodes2clusterize){
      hclusters<-seq_along(gonames)
      names(hclusters)<-gonames
    }else{
     go.dist<-similarity.to.distance(term.sim(gonames,gometric),gometric)
     go.dist<-jitter(go.dist,factor=0.1)
     hc<-hclust(as.dist(go.dist),method="average")

     if(bcutree){
      icheck<-which(diff(hc$height)<0)
      if(length(icheck)>0){
       for(ii in seq_along(icheck)){
	hc$height[icheck[ii]+1]<-hc$height[icheck[ii]]+1e-6
       }
      }
      hclusters<-cutree(hc,h=hh)
     }else{
      if(FALSE){
        cth<-cutreeHybrid(hc,distM=go.dist,minClusterSize=minClusterSize)
	cth$labels
	cth$cores
	cth$smallLabels
	cth$branches
      }else{
       cth<-cutreeDynamic(hc,distM=go.dist,method="hybrid",
                          cutHeight=hh,minClusterSize=minClusterSize)      
       hclusters<-cth
       names(hclusters)<-hc$labels
      # hclusters<-hclusters[hclusters!=0]
      } 
     }
    }
    
    mica<-c()
    for(j in seq_along(unique(hclusters))){
      mica<-c(mica,biocat.goMICA(names(hclusters)[hclusters==j]))
    }
        
    a[["leaves"]]<-hclusters
    a[["mica"]]  <-mica
    x[[ccat]]<-a
 }
 return(x)})


 #identifico clusters que no tienen ninguna categoria
 #ia0<-unlist(lapply(aa,function(x){any(unlist(lapply(x,function(y){length(unlist((y))}))>0)}))

 
 lbioTheme$paramsGOcluster<-c(gometric=gometric,bcutree=bcutree,minClusterSize=minClusterSize,
             micaBranchMinIC=micaBranchMinIC,minGOnodes2clusterize=minGOnodes2clusterize,onlyLeaves=onlyLeaves)
 lbioTheme[["clusters"]]<-aa     
 return(lbioTheme)
}

## Genera una tabla de reporte del analisys de overexpression
## para la categoria seleccionada
## lbioTheme y hg DEBEN estar definidos correctamente
## input
##  iclus:nro cluster
##  cat:  categoria
## out:
##   data.frame con info sobre categorias seleccionadas
##   y temas biologicos
biocat.Report<-function(iclus,cat=c("BP","MF","CC","KEGG")[1]){
   a<-hg[[iclus]]
   if(lbioTheme[["params"]]["onlyLeaves"]==1){
    b<-lbioTheme[["clusters"]][[iclus]][[cat]]$leaves
   }else{
    b<-lbioTheme[["clusters"]][[iclus]][[cat]]$ids
   }
   if(is.null(a[[cat]])) next
   if(is.na(a[[cat]])) next
   cat("Cluster:",iclus," # genes:",length(a$gset),"\n")
   if(cat=="KEGG"){
    anot<-annotation(hg[[1]]$KEGG)
    eenv<-get(paste(orgdb,"SYMBOL",sep=""))                            
    require(paste(anot,"db",sep="."),character.only=TRUE)
#     ddf<-data.frame(cbind(
#                    c("TestName:","GeneSetSize","GeneUniverseSize",
#                      "#Tested:",paste("#Signif (pv<",pvalueCutoff(a[[cat]]),"):",sep=""),
#                      "Annotation"),      
#                    c(paste(paste(a[[cat]]@testName,collapse="-"),pvalueCutoff(a[[cat]]),sep="-"),
#                      geneMappedCount(a[[cat]]),
#                      universeMappedCount(a[[cat]]),
#                      length(nodes(a[[cat]]@goDag)),
#                      a[[cat]]@annotation)))
#      colnames(ddf)<-c(" "," ")    
    gbc<-geneIdsByCategory(a$KEGG)
    laux<-lapply(gbc,function(x){
         #aux<-unlist(lapply(mget(x,eenv,ifnotfound=NA),function(y){return(paste(y,collapse=":"))}))
         aux<-unlist(lapply(mget(x,eenv,ifnotfound=NA),function(y){return(y[1])}))                 
         aux[is.na(aux)]<-names(aux)[is.na(aux)]
         return(aux)
        })
    genes<-unlist(lapply(laux[sigCategories(a$KEGG)],function(x){paste(x,collapse=", ")}))
    df<-cbind(summary(a$KEGG),genes)
    rownames(df)<-NULL
   }else{
    df<-data.frame(cbind(
                   c("TestName:","GeneSetSize","GeneUniverseSize",
                     "#Tested:",paste("#Signif (pv<",a[[cat]]@pvalueCutoff,"):",sep=""),
                     "#Leave biotheme concepts","#MICA branches",
                     
                     "Annotation"),      
                   c(paste(paste(a[[cat]]@testName,collapse="-"),a[[cat]]@testDirection,sep="-"),
                     geneMappedCount(a[[cat]]),
                     universeMappedCount(a[[cat]]),
                     length(nodes(a[[cat]]@goDag)),
                     paste(length(sigCategories(a[[cat]],a[[cat]]@pvalueCutoff))),                     
                     length(b),length(unique(b)),
                     a[[cat]]@annotation)))
     colnames(df)<-c(" "," ")    
    
     df1<-data.frame(names(lbioTheme$params),lbioTheme$params)
     colnames(df1)<-c(" "," ")
    
     df2<-data.frame(names(lbioTheme$paramsGOcluster),lbioTheme$paramsGOcluster)
     colnames(df2)<-c(" "," ")
     ddf<-rbind(df,rbind(df1,df2))
     rownames(ddf)<-NULL
    } 
    return(ddf)
}


## Genera plot de GOgraph coloreado
## lbioTheme DEBE estar definidos correctamente
## input:
##   iclus: ncluster
##   cat: categoria a plotear
## output:
##   grafico GO DAG coloreado
biocat.goGraph<-function(iclus,cat=c("BP","MF","CC")[1]){
     gonames<-lbioTheme[["clusters"]][[iclus]][[cat]]$leaves
     if(length(gonames)<1){
       cat("No leaves categories to plot","\n")
       return(NULL)
     }  
     if(cat=="BP") eenv=GOBPPARENTS
     if(cat=="MF") eenv=GOMFPARENTS
     if(cat=="CC") eenv=GOCCPARENTS
     gogr<-GOGraph(names(gonames),eenv)
#      ccol<-heat.colors(length(unique(gonames)))[gonames]
#      ccol<-terrain.colors(length(unique(gonames)))[gonames]
#      ccol<-topo.colors(length(unique(gonames)))[gonames]
     ccol<-rainbow(length(unique(gonames)))[gonames]
     names(ccol)<-names(gonames)
     lnattrs<-list(fillcolor=ccol)  
     plot(gogr,nodeAttrs=lnattrs)
}

## Dado un grupo de genes, encuentra el concepto MICA asociado
## input:
##   gos: vector con nombre de nodos
## output: vector nombrado con el valo de IC del mica y el nombre respsectivo
biocat.goMICA<-function(gos){
  cat<-getGOOntology(gos[1])
  if(cat=="BP") eenv=GOBPANCESTOR
  if(cat=="MF") eenv=GOMFANCESTOR
  if(cat=="CC") eenv=GOCCANCESTOR
  gosAnc<-mget(gos,eenv)
  gosAllAnc<-unique(unlist(gosAnc))
  gosAllAnc<-gosAllAnc[!(gosAllAnc=="all")]
  icAllAnc<-diag(term.sim(gosAllAnc,"Resnik"))
  
  mica<-gosAnc[[gos[1]]]
  for(j in seq(along=gos)[-1]){
    mica<-intersect(mica,gosAnc[[gos[j]]])
  }
  mica<-mica[!(mica=="all")]
  mica<-names(which.max(icAllAnc[mica]))[1]
  return(icAllAnc[mica])
}

if(FALSE){ #aver...
 ggos<-lbioTheme[[1]][[7]]$BP$leaves
 ggos<-names(ggos)[ggos==1]
 biocat.goMICA(ggos)
}

## Dado un grafo completo 'gg', extrae el subgrafo correspondeinte
## al cluster iclus, y lo colorea de acuerdo a los clustersGO
## de conceptos leaves sobreexpresados
## input:
##   gg:  grafo igraph
##   iclust: cluster a extraer
##   cat: categoria a utilizar para colorear
##   includeNoLeaves: si TRUE se colorean ademas los genes mapeados a 
##                   nodos ancestros 
##   kSizeScale: se escalean tamanios deacuerdo al grado?
##   size0,(size1): tamanio fijo nodos genes no incluidos (incluidos) en 
##                  categorias sobreexpresadas (si kSizeScale=TRUE se ignoran)
##   smin,(smax): si kSizeScale=TRUE, el tamanio de c/nodo resulta proporcional
##                al grado del nodo
##   layout: layout previamente calculado para el grafo
##   bOutLayout: devuele layout a la salida?
## output:
##    si bOutLayout es TRUE devuelvce el layout utilizado
biocat.ClusterGraph<-function(gg,iclus,cat=c("BP","MF","CC")[1],
                              inputLayout=NULL,bOutLayout=FALSE,
                              includeNoLeaves=TRUE,
                              size0=5,size1=10,kSizeScale=FALSE,smin=5,smax=10){
 x<-lbioTheme[["clusters"]][[iclus]]
 goleaves<-x[[cat]]$leaves
 #pongo a los no-clustereados en clusters separados
 i0<-which(goleaves==0)
 if(length(i0)>0){
  offset<-length(u)-1
  goleaves[i0]<-offset+seq_along(i0)
 }
 tt<-table(goleaves)
 
 #que genes aparecen en cada goCluster?
 gogenes<-list()
 for(ii in 1:length(tt)){
  clus<-ii
  gos<-names(goleaves[goleaves==clus])
  gogenes[[as.character(ii)]]<-unique(unlist(geneIdsByCategory(hg[[iclus]][[cat]])[gos]))             
 }
 genesgo<-reverseSplit(gogenes)
 if(any(!names(genesgo)%in%hg[[iclus]]$gset)) warning("OhOh")   ## check!!

 ## estaba con un BUG??...lo corrijo
 #ig <- match(names(hg[[iclus]]$gset),V(gg)$name)- igraph.offset
 ## igual esto tampoco era muy cool
 #ig <- match(hg[[iclus]]$gset,V(gg)$name)- igraph.offset
 #sgg<-induced.subgraph(gg,ig)
 ## ahora lo hago asi
 sgg<-induced.subgraph(gg,hg[[iclus]]$gset)
 
 if(!is.null(inputLayout)){
   l <- inputLayout
 }else{
  l <- layout.drl(sgg, options=list(simmer.attraction=0))
 }

 # gSymbol <-names(hg[[iclus]]$gset[match(names(genesgo),hg[[iclus]]$gset)])
 gSymbol <-hg[[iclus]]$gset[match(names(genesgo),hg[[iclus]]$gset)]
 igSymbol<-match(gSymbol,V(sgg)$name) - igraph.offset
 if(includeNoLeaves){
  go0<-x[[cat]]$ids[!x[[cat]]$ids%in%names(x[[cat]]$leave)]
  noLeaveGenes<-geneIdsByCategory(hg[[iclus]][[cat]])[go0]
  noLeaveGenes<-unique(unlist(noLeaveGenes))
  noLeaveGenes<-noLeaveGenes[!noLeaveGenes%in%names(genesgo)]
  #gSymbol0 <-names(hg[[iclus]]$gset[match(noLeaveGenes,hg[[iclus]]$gset)])
  gSymbol0 <-hg[[iclus]]$gset[match(noLeaveGenes,hg[[iclus]]$gset)]
  igSymbol0<-match(gSymbol0,V(sgg)$name)-igraph.offset
 }
 
 if(kSizeScale){
  V(sgg)$size <-smax*(igraph::degree(sgg) - smin)/(smax-smin)+smin
 }else{
  V(sgg)$size <-size0
 } 
 V(sgg)[igSymbol]$size<-size1
 V(sgg)$color          <-rgb(.4,.4,.4)

 #color de genes mapeados a nodos no-leaves
 if(includeNoLeaves){
  ccol<-"black"
  V(sgg)[igSymbol0]$color<-ccol
 }
 
 #color de genes mapeados a nodos leaves
 #BUG: que hago cuando un gen esta anotado a mas de una cat GO????
 # por ahora tomo el 1
 ccol<-as.numeric(unlist(lapply(genesgo,function(x){x[1]})))
 ccol<-rainbow(length(unique(unlist(genesgo))))[ccol]
 V(sgg)[igSymbol]$color<-ccol
 
 #V(sgg)$label<-V(sgg)$name
 E(sgg)$color<-rgb(.1,.1,.1,.25)
 plot(sgg, layout=l,
     edge.color=E(sgg)$color, vertex.label="",
     vertex.size=V(sgg)$size, edge.width=2)

 if(bOutLayout) return(l)
}

if(FALSE){
 biocat.ClusterGraph(gg,10,includeNoLeaves=TRUE)
}

##
## Grafico de Treatment Pattern Activity
 biocat.patternGraph<-function(pat){
  #pat<-c(0,0,0,1,0,1,1,0,0,1)
  #names(pat)<-rep(c("shoots.heat","shoots.cold"),5)
  plot(seq_along(pat),pat,typ="n",xlab="",ylab="",axes=FALSE)
  text(seq_along(pat),rep(0,length(pat)),names(pat),col=pat+1,srt=90,adj=0.1,
      font=pat+1)
 }

## Export resultados html
## Usa lbioTheme, hg, gg
## Falta KEGG
biocat.exportHtml.v0<-function(iclus,fname="clusterMCODE",cat=c("BP","MF","CC","KEGG")[1:3],
                            bplot=FALSE,kscale=1){
                            
 for(icat in seq_along(cat)){ 
  ccat<-cat[icat]
  if(length(lbioTheme[["clusters"]][[iclus]][[ccat]]$leaves)==0) return()
  
  fnameHtml<-paste(fname,iclus,ccat,"html",sep=".")
  phtml  <-openPage(fnameHtml)
   
  #tabla con datos de analisis de overrepresentation
  df.repo<-biocat.Report(iclus,ccat)
  
  hw1<-hwrite(df.repo[1:8,],border=0,col.names=FALSE,row.names=FALSE);
  hw2<-hwrite(df.repo[9:19,],border=0,col.names=FALSE,row.names=FALSE);
  
  hwrite(cbind(hw1,hw2),phtml,border=1,col.names=FALSE);
  
  #DAG coloreado
  ffname1<-paste(fname,"GOdag",iclus,ccat,"png",sep=".")
  if(bplot){
   biocat.goGraph(iclus,ccat)
   dev.copy(png,file=ffname1,width=600*kscale,height=400*kscale)
   dev.off() 
  }else{
   png(file=ffname1,width=600*kscale,height=400*kscale)
   biocat.goGraph(iclus,ccat)
   dev.off() 
  }
  #hwriteImage(ffname1,page=phtml)
 
  #Grafo interacciones
  ffname2<-paste(fname,"subgraph",iclus,ccat,"png",sep=".")
  if(bplot){
   if(icat==1){
    cLayout<-biocat.ClusterGraph(gg,iclus,ccat,inputLayout=NULL,bOutLayout=TRUE)
   }else{
    biocat.ClusterGraph(gg,iclus,ccat,inputLayout=cLayout,bOutLayout=FALSE)
   }
   dev.copy(png,file=ffname2,width=500*kscale,height=500*kscale)
   dev.off() 
  }else{
   png(file=ffname2,width=500*kscale,height=500*kscale)
   if(icat==1){
    cLayout<-biocat.ClusterGraph(gg,iclus,ccat,inputLayout=NULL,bOutLayout=TRUE)
   }else{
    biocat.ClusterGraph(gg,iclus,ccat,inputLayout=cLayout)
   }
   dev.off() 
  }  
  #hwriteImage(ffname2,page=phtml)
  
  himg<-hwriteImage(c(ffname1,ffname2),table=FALSE)
  hwrite(himg,phtml,br=TRUE,center=TRUE)
  
  #Tablas GO's
  xx<-lbioTheme[["clusters"]][[iclus]][[ccat]]
  hclusters<-xx$leaves
  summary.hg<-summary(hg[[iclus]][[ccat]])
  summary.hg[,2:4]<-signif(summary.hg[,2:4],2)
  ccol<-rainbow(length(unique(hclusters)))
  ccol<-substr(ccol,1,7)
  for(j in seq_along(unique(hclusters))){
     ih<-which(hclusters==unique(hclusters)[j])
     x<-names(hclusters)[ih]
     ddf<-summary.hg[match(x,summary.hg[,1]),]
     
     saux<-paste(names(xx$mica[j]),getGOTerm(names(xx$mica[j])),sep="-")
     hwrite(saux,page=phtml,heading=3)
     hwrite(ddf,page=phtml,row.bgcolor=ccol[j])                        
  }

  
  closePage(phtml) 
 } 
}     


## Export resultados html
## Usa lbioTheme, hg, gg
## Falta KEGG
## 
## Report                 | Colored DAG
## colored gg subgraph    | MICAs labels
## GO tables
biocat.exportHtml.v1<-function(iclus,fname="clusterMCODE",cat=c("BP","MF","CC","KEGG")[1:3],
                            bplot=FALSE,kscale=1,orgdb=NULL){
 if(!is.null(orgdb)){
  require(paste(orgdb,"db",sep="."),character.only=TRUE)                            
  eenv<-get(paste(orgdb,"SYMBOL",sep=""))                            
  ncbi="http://www.ncbi.nlm.nih.gov/gene?term="
 } 

 for(icat in seq_along(cat)){ 
  ccat<-cat[icat]
  if(length(lbioTheme[["clusters"]][[iclus]][[ccat]]$leaves)==0) return()
  
  gbycat<-geneIdsByCategory(hg[[iclus]][[ccat]])
  
  fnameHtml<-paste(fname,iclus,ccat,"html",sep=".")
  phtml  <-openPage(fnameHtml)
   
  #tabla con datos de analisis de overrepresentation
  df.repo<-biocat.Report(iclus,ccat)
  
  hw.report<-hwrite(df.repo[1:8,],border=0,col.names=FALSE,row.names=FALSE);
  hw.params<-hwrite(df.repo[9:19,],border=0,col.names=FALSE,row.names=FALSE);
  
  
  #DAG coloreado
  ffname1<-paste(fname,"GOdag",iclus,ccat,"png",sep=".")
  if(bplot){
   biocat.goGraph(iclus,ccat)
   dev.copy(png,file=ffname1,width=600*kscale,height=400*kscale)
   dev.off() 
  }else{
   png(file=ffname1,width=600*kscale,height=400*kscale)
   biocat.goGraph(iclus,ccat)
   dev.off() 
  }
  hw.imgdag<-hwriteImage(ffname1)
  
  #Grafo interacciones
  ffname2<-paste(fname,"subgraph",iclus,ccat,"png",sep=".")
  if(bplot){
   if(icat==1){
    cLayout<-biocat.ClusterGraph(gg,iclus,ccat,inputLayout=NULL,bOutLayout=TRUE)
   }else{
    biocat.ClusterGraph(gg,iclus,ccat,inputLayout=cLayout,bOutLayout=FALSE)
   }
   dev.copy(png,file=ffname2,width=500*kscale,height=500*kscale)
   dev.off() 
  }else{
   png(file=ffname2,width=500*kscale,height=500*kscale)
   if(icat==1){
    cLayout<-biocat.ClusterGraph(gg,iclus,ccat,inputLayout=NULL,bOutLayout=TRUE)
   }else{
    biocat.ClusterGraph(gg,iclus,ccat,inputLayout=cLayout)
   }
   dev.off() 
  }  
  hw.imggraph<-hwriteImage(ffname2)
  
  
  hwrite(cbind(hw.report,hw.imggraph),phtml,border=1,col.names=FALSE);


  #Tablas GO's
  xx<-lbioTheme[["clusters"]][[iclus]][[ccat]]
  hclusters<-xx$leaves
  summary.hg<-summary(hg[[iclus]][[ccat]])
  summary.hg[,2:4]<-signif(summary.hg[,2:4],2)
  ccol<-rainbow(length(unique(hclusters)))
  ccol<-substr(ccol,1,7)
  hw.mica<-list()
  hw.table<-list()
  n<-length(unique(hclusters))
  for(j in 1:n){
     ih<-which(hclusters==unique(hclusters)[j])
     x<-names(hclusters)[ih]
     ddf<-summary.hg[match(x,summary.hg[,1]),]
     
     if(!is.null(orgdb)){
      aa<-lapply(gbycat[ddf[,1]],function(x){
                       xx<-mget(x,eenv)
                       #yy<-paste("link=",ncbi,names(xx),sep="")
                       paste(unlist(xx),collapse=" ")
                    })
      ddf<-cbind(ddf,unlist(aa))                   
     }

     #saux<-paste(names(xx$mica[j]),getGOTerm(names(xx$mica[j])),sep="-")
     saux<-getGOTerm(names(xx$mica[j]))[[1]]
     hw.mica[[j]]<-hwrite(saux)
     hw.table[[j]]<-hwrite(ddf,row.bgcolor=ccol[j],row.names=FALSE,border=0)                        
  }

  #escribo html
  phtml  <-openPage(fnameHtml)
  
  hw.micacolor<-hwrite(rep("*",n),dim=c(n,1),bgcolor=ccol,border=0)
  hw.micatable<-hwrite(unname(unlist(hw.mica)),dim=c(n,1),border=0)
  hw.legend<-hwrite(c(hw.micacolor,hw.micatable),border=0)
  
  hwrite(cbind(hw.report,hw.imgdag),phtml,border=0,col.names=FALSE);
  hwrite(cbind(hw.imggraph,hw.legend),phtml,border=0,col.names=FALSE);
  
  for(i in 1:n){
    hwrite(hw.mica[[i]],phtml,heading=3)
    hwrite(hw.table[[i]],phtml,border=0,br=TRUE)
  }
  
 
  closePage(phtml) 
 } 
}     



## Export resultados html
## Usa lbioTheme, hg, dataX
## Falta KEGG
## input
##   iclus: cluster a considerar
##   pCluster: dataFrame con info sobre id-tissue-treatment-color-size
##   fname:    pattern para files
##   cat: categorias GO que seran exportadas
##   bplot:
##   kscale: scaling de png's (kscale*(600,400))
##   orgdb : nombre libreria de anotacion para el organismo (e.g. org.At.tair)
##
## output: htmlpage
##   cluster expresion profile| correlations
##   Report                   | Colored DAG
##   GO tables
biocat.exportHtml.v2<-function(iclus,pCluster=pCluster,fname="cluster",cat=c("BP","MF","CC","KEGG")[1:3],
                            bplot=FALSE,kscale=1,orgdb=NULL){
 if(!is.null(orgdb)){
  require(paste(orgdb,"db",sep="."),character.only=TRUE)                            
  eenv<-get(paste(orgdb,"SYMBOL",sep=""))                            
  ncbi="http://www.ncbi.nlm.nih.gov/gene?term="
 } 

 for(icat in seq_along(cat)){ 
  ccat<-cat[icat]
  if(length(lbioTheme[["clusters"]][[iclus]][[ccat]]$leaves)==0) return()
  
  gbycat<-geneIdsByCategory(hg[[iclus]][[ccat]])
  
  fnameHtml<-paste(fname,iclus,ccat,"html",sep=".")
   
  ######## 
  ## tabla con datos de analisis de overrepresentation
  df.repo<-biocat.Report(iclus,ccat)
  
  hw.report<-hwrite(df.repo[1:8,],border=0,col.names=FALSE,row.names=FALSE);
  hw.params<-hwrite(df.repo[9:19,],border=0,col.names=FALSE,row.names=FALSE);
  
  ###########
  ## Perfiles de Expresion
    treatment<-pCluster[iclus,"treatment"]
    tissue   <-pCluster[iclus,"tissue"]
    
    if(tissue=="Cell"){
     slides<-as.character(pheno[which(as.character(pheno[,"TYPE"])==treatment &
                                      as.character(pheno[,"SRC"])=="CellCulture"),"ID"])
    }else{
     slides<-as.character(pheno[which(as.character(pheno[,"TYPE"])==treatment &
                                      as.character(pheno[,"SRC"])==tissue),"ID"])
    }
    psets<-names(hg[[iclus]]$gset)
    datax<-apply(dataX[psets,slides],1,varnorm)
  
    MEs<-clusters[[tissue]][[treatment]][["MEs"]]
    #ime<-grep(paste("ME",pCluster[iclus,"color"],sep=""),colnames(MEs))
    ime<-which(paste("ME",pCluster[iclus,"color"],sep="")==colnames(MEs))
    me <-varnorm(MEs[,ime])
    
    #Si no existe...genero la figura de expresion para el cluster
    fname0<-paste(paste(pCluster[iclus,],collapse="."),".png",sep="")
    #if(!(fname0%in%dir())){
    if(TRUE){
     png(file=fname0,width=800,height=400)
     layout(matrix(1:2,1,2))
     matplot(datax,typ="l",col=rgb(0,.4,0,0.2),
            main=paste(pCluster[iclus,],collapse="."))
     lines(me,col="gray",lwd=5)
    
     #correlaciones con ME
     ccor<-cor(me,datax)
     plot(rank(ccor[1,]),ccor[1,],ylim=c(-1,1),
          xlab="rank cor(gen,ME)",ylab="cor(gen,ME)",cex=0.8)
     abline(h=seq(-1,1,.25),lty=3,col="gray")
     m<-mean(abs(ccor))
     abline(h=c(m),col=1,lwd=2)
     text(1,m*0.9,paste("cor",signif(m,2)),cex=0.8,col=1,pos=4)
     if(sum(ccor>0)>0){
      mp<-mean(ccor[ccor>0])
      abline(h=c(mp),col=c(2),lwd=2)
      text(1,mp*1.1,paste("cor+",signif(mp,2)),cex=0.8,col=2,pos=4)
     } 
     if(sum(ccor<0)>0){
      mm<-mean(ccor[ccor<0])
      abline(h=c(mm),col=c(3),lwd=2)
      text(1,mm*0.9,paste("cor-",signif(mm,2)),cex=0.8,col=3,pos=4)
     } 
     dev.off()
    }
    hw.genex<-hwriteImage(fname0)
   


  
  
  #####
  ## DAG coloreado
  ffname1<-paste(fname,"GOdag",iclus,ccat,"png",sep=".")
  if(bplot){
   biocat.goGraph(iclus,ccat)
   dev.copy(png,file=ffname1,width=600*kscale,height=400*kscale)
   dev.off() 
  }else{
   png(file=ffname1,width=600*kscale,height=400*kscale)
   biocat.goGraph(iclus,ccat)
   dev.off() 
  }
  hw.imgdag<-hwriteImage(ffname1)
  
  ########
  ##Tablas GO's
  xx<-lbioTheme[["clusters"]][[iclus]][[ccat]]
  hclusters<-xx$leaves
  summary.hg<-summary(hg[[iclus]][[ccat]])
  summary.hg[,2:4]<-signif(summary.hg[,2:4],2)
  ccol<-rainbow(length(unique(hclusters)))
  ccol<-substr(ccol,1,7)
  hw.mica<-list()
  hw.table<-list()
  n<-length(unique(hclusters))
  for(j in 1:n){
     ih<-which(hclusters==unique(hclusters)[j])
     x<-names(hclusters)[ih]
     ddf<-summary.hg[match(x,summary.hg[,1]),]
     
     if(!is.null(orgdb)){
      aa<-lapply(gbycat[ddf[,1]],function(x){
                       xx<-unlist(lapply(mget(x,eenv),function(y){y[1]}))
                       ina<-which(is.na(xx))
                       xx[ina]<-names(xx)[ina]
                       #yy<-paste("link=",ncbi,names(xx),sep="")
                       paste(unlist(xx),collapse=" ")
                    })
      ddf<-cbind(ddf,unlist(aa))                   
     }

     #saux<-paste(names(xx$mica[j]),getGOTerm(names(xx$mica[j])),sep="-")
     saux<-getGOTerm(names(xx$mica[j]))[[1]]
     hw.mica[[j]]<-hwrite(saux)
     hw.table[[j]]<-hwrite(ddf,row.bgcolor=ccol[j],row.names=FALSE,border=0)                        
  }

  #escribo html
  phtml  <-openPage(fnameHtml)
  
  hw.micacolor<-hwrite(rep("*",n),dim=c(n,1),bgcolor=ccol,border=0)
  hw.micatable<-hwrite(unname(unlist(hw.mica)),dim=c(n,1),border=0)
  hw.legend<-hwrite(c(hw.micacolor,hw.micatable),border=0)
  
  hwrite(hw.genex,phtml,border=0)
  hwrite(cbind(hw.report,hw.imgdag),phtml,border=0,col.names=FALSE);
  
  
  for(i in 1:n){
    hwrite(hw.mica[[i]],phtml,heading=3)
    hwrite(hw.table[[i]],phtml,border=0,br=TRUE)
  }
  
 
  closePage(phtml) 
 } 
}     


## Export resultados html
## Usa lbioTheme, hg,
## Falta KEGG
## input
##   iclus: cluster a considerar
##   fname:    pattern para files
##   cat: categorias GO que seran exportadas
##   bplot:
##   kscale: scaling de png's (kscale*(600,400))
##   orgdb : nombre libreria de anotacion para el organismo (e.g. org.At.tair)
##
## output: htmlpage
##   Report                   | Colored DAG
##   GO tables
biocat.exportHtml.v3<-function(iclus,fname="cluster",cat=c("BP","MF","CC","KEGG")[1:3],
                            bplot=FALSE,kscale=1,orgdb=NULL){
 if(!is.null(orgdb)){
  require(paste(orgdb,"db",sep="."),character.only=TRUE)                            
  eenv<-get(paste(orgdb,"SYMBOL",sep=""))                            
  ncbi="http://www.ncbi.nlm.nih.gov/gene?term="
 } 

 for(icat in seq_along(cat)){ 
  ccat<-cat[icat]
  if(length(lbioTheme[["clusters"]][[iclus]][[ccat]]$leaves)==0) return()
  
  gbycat<-geneIdsByCategory(hg[[iclus]][[ccat]])
  
  fnameHtml<-paste(fname,iclus,ccat,"html",sep=".")
   
  ######## 
  ## tabla con datos de analisis de overrepresentation
  df.repo<-biocat.Report(iclus,ccat)
  
  hw.report<-hwrite(df.repo[1:8,],border=0,col.names=FALSE,row.names=FALSE);
  hw.params<-hwrite(df.repo[9:19,],border=0,col.names=FALSE,row.names=FALSE);
  
  
  #####
  ## DAG coloreado
  ffname1<-paste(fname,"GOdag",iclus,ccat,"png",sep=".")
  if(bplot){
   biocat.goGraph(iclus,ccat)
   dev.copy(png,file=ffname1,width=600*kscale,height=400*kscale)
   dev.off() 
  }else{
   png(file=ffname1,width=600*kscale,height=400*kscale)
   biocat.goGraph(iclus,ccat)
   dev.off() 
  }
  hw.imgdag<-hwriteImage(ffname1)
  
  ########
  ##Tablas GO's
  xx<-lbioTheme[["clusters"]][[iclus]][[ccat]]
  hclusters<-xx$leaves
  summary.hg<-summary(hg[[iclus]][[ccat]])
  summary.hg[,2:4]<-signif(summary.hg[,2:4],2)
  ccol<-rainbow(length(unique(hclusters)))
  ccol<-substr(ccol,1,7)
  hw.mica<-list()
  hw.table<-list()
  n<-length(unique(hclusters))
  for(j in 1:n){
     ih<-which(hclusters==unique(hclusters)[j])
     x<-names(hclusters)[ih]
     ddf<-summary.hg[match(x,summary.hg[,1]),]
     
     if(!is.null(orgdb)){
      aa<-lapply(gbycat[ddf[,1]],function(x){
                       xx<-unlist(lapply(mget(x,eenv),function(y){y[1]}))
                       ina<-which(is.na(xx))
                       xx[ina]<-names(xx)[ina]
                       #yy<-paste("link=",ncbi,names(xx),sep="")
                       paste(unlist(xx),collapse=" ")
                    })
      ddf<-cbind(ddf,unlist(aa))                   
     }

     #saux<-paste(names(xx$mica[j]),getGOTerm(names(xx$mica[j])),sep="-")
     saux<-getGOTerm(names(xx$mica[j]))[[1]]
     hw.mica[[j]]<-hwrite(saux)
     hw.table[[j]]<-hwrite(ddf,row.bgcolor=ccol[j],row.names=FALSE,border=0)                        
  }

  #escribo html
  phtml  <-openPage(fnameHtml)
  
  hw.micacolor<-hwrite(rep("*",n),dim=c(n,1),bgcolor=ccol,border=0)
  hw.micatable<-hwrite(unname(unlist(hw.mica)),dim=c(n,1),border=0)
  hw.legend<-hwrite(c(hw.micacolor,hw.micatable),border=0)
  
  #hwrite(hw.genex,phtml,border=0)
  hwrite(cbind(hw.report,hw.imgdag),phtml,border=0,col.names=FALSE);
  
  
  for(i in 1:n){
    hwrite(hw.mica[[i]],phtml,heading=3)
    hwrite(hw.table[[i]],phtml,border=0,br=TRUE)
  }
  
 
  closePage(phtml) 
 } 
}     


 ## Export resultados html
## Usa lbioTheme, hg,
## Falta KEGG
## input
##   iclus: cluster a considerar
##   fname:    pattern para files
##   cat: categorias GO que seran exportadas
##   bplot:
##   kscale: scaling de png's (kscale*(600,400))
##   orgdb : nombre libreria de anotacion para el organismo (e.g. org.At.tair)
##
## output: htmlpage
##   Treatment Activity Pattern
##   Report                   | Colored DAG
##   GO tables
biocat.exportHtml.v4<-function(iclus,pat,fname="cluster",cat=c("BP","MF","CC","KEGG")[1:3],
                            bplot=FALSE,kscale=1,orgdb=NULL){
 if(!is.null(orgdb)){
  require(paste(orgdb,"db",sep="."),character.only=TRUE)                            
  eenv<-get(paste(orgdb,"SYMBOL",sep=""))                            
  ncbi="http://www.ncbi.nlm.nih.gov/gene?term="
 } 

 for(icat in seq_along(cat)){ 
  ccat<-cat[icat]
  if(length(lbioTheme[["clusters"]][[iclus]][[ccat]]$leaves)==0) return()
  
  gbycat<-geneIdsByCategory(hg[[iclus]][[ccat]])
  
  fnameHtml<-paste(fname,iclus,ccat,"html",sep=".")
   
  ########
  ## Treatment Activity Pattern
  ffname11<-paste(fname,"TAP",iclus,ccat,"png",sep=".") 
  kkscale<-kscale*0.6
  if(bplot){
   biocat.patternGraph(pat)
   dev.copy(png,file=ffname11,width=600*kkscale,height=400*kkscale)
   dev.off() 
  }else{
   png(file=ffname11,width=600*kkscale,height=400*kkscale)
   biocat.patternGraph(pat)
   dev.off() 
  }
  hw.imgpat<-hwriteImage(ffname11)
     
  ######## 
  ## tabla con datos de analisis de overrepresentation
  df.repo<-biocat.Report(iclus,ccat)
  
  hw.report<-hwrite(df.repo[1:8,],border=0,col.names=FALSE,row.names=FALSE);
  hw.params<-hwrite(df.repo[9:19,],border=0,col.names=FALSE,row.names=FALSE);
  
#    ## Tabla KEGG
#    hw.reportKEGG<-NULL
#    if(ccat== "BP"  &  !is.null(hg[[iclus]][["KEGG"]])){
#     if(length(sigCategories(  hg[[iclus]][["KEGG"]]))>0){
#      df.repoKEGG<-biocat.Report(iclus,"KEGG")
#      hw.reportKEGG<-hwrite(df.repoKEGG,border=0,col.names=FALSE,row.names=FALSE);
#     }
#    }

  
  #####
  ## DAG coloreado
  ffname1<-paste(fname,"GOdag",iclus,ccat,"png",sep=".")
  if(bplot){
   biocat.goGraph(iclus,ccat)
   dev.copy(png,file=ffname1,width=600*kscale,height=400*kscale)
   dev.off() 
  }else{
   png(file=ffname1,width=600*kscale,height=400*kscale)
   biocat.goGraph(iclus,ccat)
   dev.off() 
  }
  hw.imgdag<-hwriteImage(ffname1)
  
  ########
  ##Tablas GO's
  xx<-lbioTheme[["clusters"]][[iclus]][[ccat]]
  hclusters<-xx$leaves
  summary.hg<-summary(hg[[iclus]][[ccat]])
  summary.hg[,2:4]<-signif(summary.hg[,2:4],2)
  ccol<-rainbow(length(unique(hclusters)))
  ccol<-substr(ccol,1,7)
  hw.mica<-list()
  hw.table<-list()
  n<-length(unique(hclusters))
  for(j in 1:n){
     ih<-which(hclusters==unique(hclusters)[j])
     x<-names(hclusters)[ih]
     ddf<-summary.hg[match(x,summary.hg[,1]),]
     
     if(!is.null(orgdb)){
      aa<-lapply(gbycat[ddf[,1]],function(x){
                       xx<-unlist(lapply(mget(x,eenv),function(y){y[1]}))
                       ina<-which(is.na(xx))
                       xx[ina]<-names(xx)[ina]
                       #yy<-paste("link=",ncbi,names(xx),sep="")
                       paste(unlist(xx),collapse=" ")
                    })
      ddf<-cbind(ddf,unlist(aa))                   
     }

     #saux<-paste(names(xx$mica[j]),getGOTerm(names(xx$mica[j])),sep="-")
     saux<-getGOTerm(names(xx$mica[j]))[[1]]
     hw.mica[[j]]<-hwrite(saux)
     hw.table[[j]]<-hwrite(ddf,row.bgcolor=ccol[j],row.names=FALSE,border=0)                        
  }

  #escribo html
  phtml  <-openPage(fnameHtml)
  
  hw.micacolor<-hwrite(rep("*",n),dim=c(n,1),bgcolor=ccol,border=0)
  hw.micatable<-hwrite(unname(unlist(hw.mica)),dim=c(n,1),border=0)
  hw.legend<-hwrite(c(hw.micacolor,hw.micatable),border=0)
  
  #hwrite(hw.genex,phtml,border=0)
  hwrite(hw.imgpat,phtml)
  
#   if(!is.null(hw.reportKEGG)) hwrite(hw.reportKEGG,phtml,border=0,col.names=FALSE);
  
  hwrite(cbind(hw.report,hw.imgdag),phtml,border=0,col.names=FALSE);
  
  
  for(i in 1:n){
    hwrite(hw.mica[[i]],phtml,heading=3)
    hwrite(hw.table[[i]],phtml,border=0,br=TRUE)
  }
  
 
  closePage(phtml) 
 } 
}     

     



 #######################
 #Export resultados html
if(FALSE){
 bplot=TRUE
 for(ccat in c("BP","MF","CC","KEGG")[1:3]){
  onto<-ccat
  for(iclus in seq_along(hg)){
   fnameHtml<-paste("cluster",iclus,ccat,"html",sep=".")
   a<-hg[[iclus]]
   if(is.null(a[[ccat]])) next
   if(is.na(a[[ccat]])) next
   cat("Cluster:",iclus," # genes:",length(a$gset),"\n")
   a$BP
   a$MF
   a$CC
   a$KEGG
   df<-data.frame(cbind(
                   c("TestName:","#Tested:",paste("#Signif (pv<",a[[ccat]]@pvalueCutoff,"):",sep=""),"GeneSetSize","GeneUniverseSize","Annotation"),
             
                   c(paste(paste(a[[ccat]]@testName,collapse="-"),a[[ccat]]@testDirection,sep="-"),
                     length(nodes(a[[ccat]]@goDag)),
                     paste(length(sigCategories(a[[ccat]],pv))),
                     geneMappedCount(a[[ccat]]),
                     universeMappedCount(a[[ccat]]),
                     a[[ccat]]@annotation)))
    colnames(df)<-c(" "," ")    
    phtml<-openPage(fnameHtml)
    hwrite(df,phtml,border=0);
    
    #me quedo con las categorias
    #que tengan entre umax y umin 
    ucount<-universeCounts(a[[ccat]])
    igo<-which(ucount>=uMin & ucount<uMax & geneCounts(a[[ccat]])>=minHits)
    gonames<-intersect(sigCategories(a[[ccat]],pv),names(igo))
    #geneIdsByCategory(a$BP)[gonames]
    if(length(gonames)==0)next;
    if(length(gonames)<5){
      hclusters<-seq_along(gonames)
      names(hclusters)<-gonames
    }else{
     go.dist<-similarity.to.distance(term.sim(gonames,gometric),gometric)
     hc<-hclust(as.dist(go.dist),method="average")
     bcutree<-FALSE
     if(bcutree){
      icheck<-which(diff(hc$height)<0)
      if(length(icheck)>0){
       for(ii in seq_along(icheck)){
	hc$height[icheck[ii]+1]<-hc$height[icheck[ii]]+1e-12
       }
      }
      hclusters<-cutree(hc,h=0.25)
     }else{
      if(FALSE){
        cth<-cutreeHybrid(hc,distM=go.dist,minClusterSize=minClusterSize)
	cth$labels
	cth$cores
	cth$smallLabels
	cth$branches
      }else{
       cth<-cutreeDynamic(hc,distM=go.dist,method="hybrid",minClusterSize=minClusterSize)      
       hclusters<-cth
       names(hclusters)<-hc$labels
       hclusters<-hclusters[hclusters!=0]
      } 
     }
    }
    ccol<-rainbow(length(unique(hclusters)))[hclusters]
    names(ccol)<-names(hclusters)
    if(bplot){
     layout(matrix(1:2,1,2))
     plot(hc,main=paste("Cluster",iclus, " #genes:",length(a$gset)," onto:",onto))
     #onto<-ccat
     if(onto=="BP") eenv=GOBPPARENTS
     if(onto=="MF") eenv=GOMFPARENTS
     if(onto=="CC") eenv=GOCCPARENTS
     gg<-GOGraph(gonames,eenv)
 
     lnattrs<-list(fillcolor=ccol)  
     plot(gg,nodeAttrs=lnattrs)
     fname<-paste("cluster",iclus,ccat,"png",sep=".")
     dev.copy(png,file=fname,width=800,height=500)
     dev.off() 
     hwriteImage(fname,page=phtml)
     closePage(phtml)
    }
    
    for(j in seq_along(unique(hclusters))){
     ih<-which(hclusters==unique(hclusters)[j])
     x<-names(hclusters)[ih]
     df<-data.frame(GOTerm=unlist(getGOTerm(x)),
                    nodeCount=universeCounts(a[[ccat]])[x],
                    hits=geneCounts(a[[ccat]])[x])
     other<-df #getGOTerm(x)
     saux<-"+---------------------------------+"
     if(j==1) saux=paste("Cluster:",iclus,"(",length(hg[[iclus]]$gset),"genes) - ",onto)
     htmlpage(list(names(hclusters)[ih]),
              filename  =fnameHtml,
              title     =saux,
              othernames=other,
              table.head=c("GOID","Description","nodeCount","count"),
              repository=list("go"),append=TRUE)#(j!=1))
    }
    ih<- which(!(hc$labels%in%names(hclusters)))
    if(length(ih)>0){
     x<-hc$labels[!(hc$labels%in%names(hclusters))]
     df<-data.frame(GOTerm=unlist(getGOTerm(x)),
                    nodeCount=universeCounts(a[[ccat]])[x],
                    hits=geneCounts(a[[ccat]])[x])
     other<-df #getGOTerm(x)
     saux<-"+---------------------------------+"
     if(j==1) saux=paste("Cluster:",iclus,"(",length(hg[[iclus]]$gset),"genes) - ",onto)
     htmlpage(list(names(hclusters)[ih]),
              filename  =fnameHtml,
              title     =saux,
              othernames=other,
              table.head=c("GOID","Description","nodeCount","count"),
              repository=list("go"),append=TRUE)#(j!=1))
    }

  } 
 }
}




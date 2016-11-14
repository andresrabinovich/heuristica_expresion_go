#Clusterizador de términos GO
library(GO.db)
setwd("~/doctorado/programacion/test_de_sobrerepresentacion/")

#Aplico filtros de evidencia
source("ecFilter.R")
dropEvidence <- c("IEA", "NAS", "IC", "ND", "RCA", "IEP")
initECfilter("org.At.tair",verbose=TRUE)
setECfilter(dropEvidence,verbose=TRUE)

#Cargamos todas las categorias GO con anotaciones que nos interesen
categorias_go <- as.list(org.At.tairGO2ALLTAIRS)

#Cargamos las categorias que estan en BP nada mas
categorias_bp = as.list(GOBPPARENTS) 

#Filtramos las categorias con anotaciones para quedarnos solo con las de bp
anotaciones_por_categoria <- categorias_go[names(categorias_go) %in% names(categorias_bp)]

#Puede haber un mismo gen anotado en varias tipos de evidencia (por ejemplo, AT3G24320 en IMP y en IGI en GO:0000002)
#así que dejamos cada gen solo una vez
anotaciones_por_categoria <- lapply(anotaciones_por_categoria, unique)

#Nos quedamos solo con las categorías que tengan alguno de los genes del tratamiento anotados
anotaciones <- anotaciones_por_categoria[unlist(lapply(anotaciones_por_categoria, function(x) sum(x %in% nombres_de_genes[rownames(genes)])) > 0)]

#Calculo el information content de todos, con GO:0008150 la cantidad de anotaciones que tiene el nodo raiz BP
ic <- unlist(lapply(anotaciones, function(x) -log2(length(unique(x))/length(unique(anotaciones[["GO:0008150"]])))))
ic <- c(ic, 0)
names(ic)[length(ic)] <- "all"


an<-as.list(GOBPANCESTOR)

categorias <- names(anotaciones)

#Calculamos la similaridad semantica de resnik
dyn.load("c/semsim.so")
simRes<-.Call("semsim", an[categorias], ic)

distRes <- 1/(1+simRes)
h<-hclust(as.dist(distRes))
clus<-cutreeDynamic(h, distM = distRes)

#simres<-matrix(0, nrow=length(categorias), ncol=length(categorias))
#rownames(simres) <- colnames(simres) <- categorias

#i = 1
#j = length(categorias)
#ptm <- proc.time()
#for(categoria1 in categorias){
#  for(categoria2 in categorias[i:j]){
    #simres[categoria1, categoria2] <- 1
#    simres[categoria1, categoria2] <- max(ic[unlist(an[categoria1])[unlist(an[categoria1]) %in% unlist(an[categoria2])]])
    #  }
#  i <- i + 1
#}
#simres <- t(simres) + simres
#diag(simres) <- diag(simres)/2
#proc.time() - ptm

#ptm <- proc.time()
#lapply(an[categorias], function(categoria2) { max(ic[unlist(an[categoria1])[unlist(an[categoria1]) %in% unlist(an[categoria2])]])})
#proc.time() - ptm
#Calculo la similaridad de resnik entre las categorías
#an<-as.list(GOBPANCESTOR)
#
#lapply(an[names(anotaciones)[1]], function(x) { x[lapply(an[names(anotaciones)[1:3]], function(y) { x %in% y } } )])
#ic[unlist(an[names(anotaciones)[1]])[unlist(an[names(anotaciones)[1]]) %in% unlist(an[names(anotaciones)[2]])][1]]
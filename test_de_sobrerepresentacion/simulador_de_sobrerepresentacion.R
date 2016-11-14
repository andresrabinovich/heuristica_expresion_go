simulador_de_sobrerepresentacion() <- function(){}
library(GO.db)
library(org.At.tair.db)
library(combinant)
setwd("~/doctorado/programacion/test_de_sobrerepresentacion/")

#Aplico filtros de evidencia
source("ecFilter.R")
dropEvidence <- c("IEA", "NAS", "IC", "ND", "RCA", "IEP")
initECfilter("org.At.tair",verbose=TRUE)
setECfilter(dropEvidence,verbose=TRUE)

#Cargamos todas las categorias GO con anotaciones que nos interesen
categorias_go <- as.list(org.At.tairGO2ALLTAIRS)

#Cargamos las categorias que estan en BP nada mas
categorias_bp <- as.list(GOBPPARENTS) 

#Filtramos las categorias con anotaciones para quedarnos solo con las de bp
anotaciones_por_categoria <- categorias_go[names(categorias_go) %in% names(categorias_bp)]

#Puede haber un mismo gen anotado en varias tipos de evidencia (por ejemplo, AT3G24320 en IMP y en IGI en GO:0000002)
#así que dejamos cada gen solo una vez
anotaciones_por_categoria <- lapply(anotaciones_por_categoria, unique)

#Tomamos el universo como todos los genes anotados
#universo<-unique(unlist(anotaciones_por_categoria))
universo<-sample(unique(unlist(anotaciones_por_categoria)), 2000) #Por ahora nos quedamos solo con estos 2000 simulando filtrado del experiment
I = length(universo)

#Nos quedamos solo con las categorías que tengan alguno de los genes del tratamiento anotados
anotaciones <- anotaciones_por_categoria[unlist(lapply(anotaciones_por_categoria, function(x) sum(x %in% universo)) > 0)]


intervalo <- c(150, 200)
lanotaciones_por_categoria<-lapply(anotaciones_por_categoria, length)
categorias_sobrerepresentadas <- names(anotaciones_por_categoria[lanotaciones_por_categoria <= intervalo[2] & lanotaciones_por_categoria >= intervalo[1]])


lapply()

pvalue_por_categoria(anotaciones_por_categoria[categorias_sobrerepresentadas[1]][[1]], sample(anotaciones_por_categoria[categorias_sobrerepresentadas[1]][[1]], 10), I)


pvalue_por_categoria<-function(genes_a, genes_s, I){
  Ia <- length(genes_a)
  Is <- length(genes_s)
  na <- length(intersect(genes_a, genes_s))
  return(pvalue(Ia, Is, na, I))
}

pvalue<-function(Ia, Is, na, I){
  p <- 0
  for(k in seq(na, min(Ia, Is))){
    p <- p + nCm(Ia, k)*nCm((I-Ia), (Is-k))/nCm(I, Is)
  }
  return(p)
}
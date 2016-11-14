#Directorio por defecto
setwd("~/trabajo/heuristica_expresion_go/clustering/")

#Cargamos librerías
library(genefilter) #Para el filtrado de genes K over A
library(dynamicTreeCut) #Para clusterizar con dynamic tree cut

#Leemos los archivos de datos de stress y el detalle de que es cada columna de stress (lo que llamamos tratamientos)
datos_stress             <- read.table(unz("datos_originales_wygel/AtGE_Abiostress_gcRMA.txt.zip", "AtGE_Abiostress_gcRMA.txt"),header=TRUE, as.is=TRUE)
detalles_de_tratamientos <- read.table("datos_originales_wygel/AtGE-Abiostress-sampleList.cvs",header=TRUE, as.is=TRUE)

#Ponemos como rownames los IDS de las tablas
rownames(datos_stress)             <- datos_stress[, 1]
rownames(detalles_de_tratamientos) <- detalles_de_tratamientos[, 1]

#Traemos la base de datos del chipset que nos dice a que gen corresponde cada probeset
require(ath1121501.db)

#Matchiamos los probesets a sus genes correspondientes
nombres_de_genes    <- unlist(mget(rownames(datos_stress), ath1121501ACCNUM, ifnotfound = NA))

#Descartamos todos los probesets que no matchean a ningun gen porque no sabemos que son
nombres_de_genes <- nombres_de_genes[!is.na(nombres_de_genes)]
datos_stress     <- datos_stress[names(nombres_de_genes), ]

#Clusterizamos en cada tratamiento
correlaciones_medias <- matrix(0, ncol=7, nrow=11)
colnames(correlaciones_medias) <- c("single", "complete", "average", "ward.D", "ward.D2", "mcquitty", "centroid")
rownames(correlaciones_medias) <- unique(detalles_de_tratamientos[, "TYPE"])
for(tratamiento in unique(detalles_de_tratamientos[, "TYPE"])[8]){
  print(tratamiento)
  
  #Traemos el detalle de la serie temporal de este tratamiento
  serie_temporal <- detalles_de_tratamientos[detalles_de_tratamientos[,"TYPE"] == tratamiento & detalles_de_tratamientos[,"SRC"]=="Shoots", "ID"]
  
  #Promediamos las dos réplicas en la serie temporal
  ipar <- seq(2,length(serie_temporal),2)
  inon <- ipar-1
  stress_promedio <- 0.5*(datos_stress[,serie_temporal[ipar]]+datos_stress[,serie_temporal[inon]])
  
  #Aplicamos filtros de tipo desviación standar. Nos quedamos con el quantile > cuantil 
  cuantil            <- 0.9
  sd_stress_promedio <- apply(stress_promedio, 1 ,sd)
  cuantil_elegido    <- quantile(sd_stress_promedio, cuantil)
  genes              <- stress_promedio[sd_stress_promedio > cuantil_elegido, ]
  
  #Aplicamos un filtro de tipo k over a. Pedimos un nivel de señal mayor a 4 en la mitad de los puntos
  koa   <- kOverA(0.5*ncol(genes), 4, na.rm = TRUE)
  genes <- genes[apply(genes, 1, function(x){koa(x)}), ]
  
  #Estandarizamos la serie temporal
  genes <- t(apply(genes,1,function(x){return( (x-mean(x))/sd(x) )}))
  genes_de_tratamiento_actual <- nombres_de_genes[rownames(genes)]
  
  #Vamos a sacar los genes que tienen probesets redundantes. No son splicing alternativo pero
  #de acuerdo con esto http://transvar.org/results/exp_cor/rp.shtml pueden ser probesets promiscuos
  #probe sets that are promiscuous, i.e., map onto multiple gene targets according to the TAIR probe set-to-gene mappings we used. 
  genes_de_tratamiento_actual <- genes_de_tratamiento_actual[!(genes_de_tratamiento_actual %in% genes_de_tratamiento_actual[duplicated(genes_de_tratamiento_actual)])]
  genes                       <- genes[names(genes_de_tratamiento_actual), ]

  for(metodo in colnames(correlaciones_medias)[2]){
  #Clusterizamos jerarquico y hacemos dynamic tree cut usando como distancia 0.5*(1-correlacion), con deepSplit = 1
    distancia             <- as.dist(0.5*(1 - (cor(t(genes)))))
    clusters_jerarquicos  <- hclust(distancia, method = metodo)
    clusters              <- as.table(cutreeDynamic(clusters_jerarquicos, distM = as.matrix(distancia), deepSplit = 1))
  }
}

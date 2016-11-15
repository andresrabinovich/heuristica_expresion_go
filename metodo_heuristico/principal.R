#seteo parametros HyperG
#source("http://bioconductor.org/biocLite.R")
#biocLite(c("GO.db", "GOstats", "ath1121501.db"))
library(GO.db)
require(GOstats)
require(ath1121501.db)
#biocLite(c("igraph", "GOSemSim", "ggplot2"))
library(igraph)
library(GOSemSim)
library(ggplot2)
#biocLite(c("grid", "gridExtra", "dendextend"))
library(grid)
library(gridExtra)
library(dendextend)
#biocLite(c("genefilter", "dynamicTreeCut", "VennDiagram", "STRINGdb"))
library(genefilter) #Para el filtrado de genes K over A
library(dynamicTreeCut) #Para clusterizar con dynamic tree cut
library(VennDiagram)
library(STRINGdb)

setwd("~/doctorado/programacion/heuristica_expresion_go/test_de_sobrerepresentacion/")
source("librerias.R")
source("hypertest.R")
source("similaridad_semantica.R")
source("hipotesis_nula.R")
source("~/doctorado/programacion/heuristica_expresion_go/clustering/clusterizador-paso-1.R")
setwd("~/doctorado/programacion/heuristica_expresion_go/test_de_sobrerepresentacion/")
source("correccion_bh.R")

setwd("~/doctorado/programacion/heuristica_expresion_go/metodo_heuristico/")

chip     <- "ath1121501"


ancestros   <-as.list(GOBPANCESTOR)
hijos       <-as.list(GOBPOFFSPRING)
acc         <-as.list(ath1121501ACCNUM)

#universe <- unique(unlist(as.list(ath1121501ACCNUM)))
universe <- nombres_de_genes[rownames(genes)]

genes_originales<-genes

evidencia_a_filtrar = c("IEA", "NAS", "IC", "ND", "RCA")

#Me quedo únicamente con los genes de mi universo que tengan alguna anotacion y al conjunto de genes y terminos go los llamo mi ontologia
anotaciones     <- g2at("BP", evidencia_a_filtrar)
anotaciones     <- unique(anotaciones[, c("gene_id", "go_id")])
ontologia       <- subset(anotaciones, gene_id %in% universe, select=c(gene_id, go_id))
universo        <- unique(ontologia[, "gene_id"])
#Sacamos los genes del universo que tengan mas de un probe para el mismo gen
universo        <- unlist(acc[which(acc %in% universo)])
universo        <- names(table(universo)[table(universo) == 1])
genes           <- genes_originales[names(unlist(acc[which(acc %in% universo)])), ]
rownames(genes) <- acc[rownames(genes)]
ontologia       <- subset(anotaciones, gene_id %in% universo, select=c(gene_id, go_id))

p                           <- hypertest(rownames(genes), ontologia=ontologia)

dag<-graph_from_graphnel(p$dag, name=TRUE)
ancestros_sin_all <- lapply(ancestros, function(x) x[-length(x)])
caminos<-matrix(0, ncol=vcount(dag), nrow=vcount(dag))
colnames(caminos)<-rownames(caminos)<-names(V(dag))
for(n in names(V(dag))){
  for(a in ancestros_sin_all[[n]]){
    caminos[n, a] <- length(all_simple_paths(dag, n, to = a, mode = "in"))
  }
}

ic                          <- -log2(p$anotaciones/p$anotaciones["GO:0008150"])


pd <- function(c1, c2, a){
  #if (c1 != c2 & c1 != a & c2 != a){
  #cat(c1, c2, a, "\n", sep=" ")
    return(abs(caminos[c1, a] - caminos[c2, a]))
  #}else{return(NULL)}
}

dishin <- function(a, b, ic) {
  common_ancestors <- intersect(c(ancestros_sin_all[[a]], a), c(ancestros_sin_all[[b]], b))#ancestros comunes entre a y b
  if(length(common_ancestors) == 1) return(ic[common_ancestors])
  #common_ancestors <- ic[common_ancestors][order(ic[common_ancestors], decreasing = T)]#los ordeno por ic
  #print(common_ancestors)
  dca <- vector()#aca van a estar la diferencia de caminos entre a y b hasta cada ancestro comun
  for (i in 1:length(common_ancestors)) {
    dca[common_ancestors[i]] <- pd(a, b, common_ancestors[i])
  }
  
  max_ic_diff_path <- c()
  for(n in names(table(dca))){
    max_ic_diff_path <- c(max_ic_diff_path, max(ic[names(which(dca == n))], na.rm = T))
  }
  
  #dca_order <- dca[order(dca)]#los ordeno y como ya estaban ordenados por ic desde antes me queda el de >ic primero (paso polémico)
  #max_ic_diff_path <- dca_order[1]#aca van a estar los dca propiamente dicho
  #if (length(dca_order)>2) {
    #  for (n in 2:length(dca_order)) {
      #if (dca_order[n] != dca_order[n-1]) {
        #max_ic_diff_path <- append(max_ic_diff_path, dca_order[n])
        #}
      #}
    #}else{
    #  return(1)
  #}
  
  
  return(mean(max_ic_diff_path))
}
bla <- p$terminos_go
m <- matrix(0, nrow = length(bla), ncol = length(bla))
colnames(m) <- rownames(m) <- bla
for (x in 1:length(bla)) {
  #print(bla[x])
  for(y in (x+1):length(bla)){
    #print(bla[y])
    m[x, y] <- dishin(bla[x], bla[y])
  }
}
m[lower.tri(m)] <- t(m)[lower.tri(m)]
m

simRes<-similaridad_semantica(ancestros[p$terminos_go][bla], ic)

p1 <- hist(simRes, breaks=100)                     # centered at 4
p2 <- hist(m, breaks=100)                     # centered at 4
plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,3))  # first histogram
plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,3), add=T)  # second


#Trae de string la magia del grafo
genes_detectados<-as.data.frame(nombres_de_genes[rownames(genes)])
colnames(genes_detectados)<-"genes"
n<-get_STRING_species(version="10", species_name="Arabidopsis thaliana")
string_db <- STRINGdb$new(version="10", species=3702, score_threshold=0, input_directory="" )
example1_mapped <- string_db$map(genes_detectados, "genes", removeUnmappedRows = TRUE)

#genes <- genes_originales


terminado <- FALSE
while(terminado == FALSE){
  
  
  distancia             <- as.dist(0.5*(1 - (cor(t(genes)))))
  clusters_jerarquicos  <- hclust(distancia, method = metodo)  
  clusters <- cutree(tree = clusters_jerarquicos, h = (clusters_jerarquicos$height[length(clusters_jerarquicos$height)]+clusters_jerarquicos$height[length(clusters_jerarquicos$height)-1])/2)
  
  #d1=color_branches(dend = clusters_jerarquicos, h=(clusters_jerarquicos$height[length(clusters_jerarquicos$height)]+clusters_jerarquicos$height[length(clusters_jerarquicos$height)-1])/2, col = c("red", "blue"))
  #plot(d1)
  table(clusters)
  p                             <- vector("list", length(unique(clusters)))
  ic                            <- vector("list", length(unique(clusters)))
  terminos_sobrerepresentados   <- vector("list", length(unique(clusters)))
  distRes                       <- vector("list", length(unique(clusters)))
  j = i = 1
  #plot(rowMeans(t(genes[clusters==i, ])), typ="b")
  #matplot(t(genes[clusters==i, ]), typ="b")
  
  for(i in unique(clusters)){
    p[[i]]                           <- hypertest(nombres_de_genes[rownames(genes)][clusters == i], ontologia=ontologia)
    terminos_sobrerepresentados[[i]] <- names(p[[i]]$pvalues[cbh(p[[i]]$pvalues, 0.05)])  
    ic[[i]]                          <- -log2(p[[i]]$anotaciones/p[[i]]$anotaciones["GO:0008150"])
    if(length(terminos_sobrerepresentados[[i]]) > 0){
      distRes[[i]] <- 1/(1+similaridad_semantica(ancestros[p[[i]]$terminos_go][terminos_sobrerepresentados[[i]]], ic[[i]]))
    }else{
      distRes[[i]] <- NA
    }
    
    if(nrow(genes[clusters == i, ]) < 50){
      terminado = TRUE
    }else{
      if(length(terminos_sobrerepresentados[[i]]) < 2){
        if(length(terminos_sobrerepresentados[[i]]) == 1){
          if(ic[[i]][terminos_sobrerepresentados[[i]]] > quantile(ic[[i]], c(0.5))){
            terminado = TRUE
          }
        }
      }else{
        distancia<-distRes[[i]]
        h<-hclust(as.dist(distancia[terminos_sobrerepresentados[[i]], terminos_sobrerepresentados[[i]]]))
        clus<-cutreeDynamic(h, distM = distancia[terminos_sobrerepresentados[[i]], terminos_sobrerepresentados[[i]]], minClusterSize = 3)
        clus<-cutree(tree = h, h = 0.3) 
        a<-c()
        lim <- max(as.integer(names(table(clus))))
      
        for(j in 1:(lim-1)){
          genes_a<- p[[i]]$genes_blancos[p[[i]]$genes_blancos %in% unique(p[[i]]$mi_ontologia[which(p[[i]]$mi_ontologia[, 2] %in% terminos_sobrerepresentados[[i]][clus==j]), 1])]
          for(k in (j+1):lim){
            genes_b<- p[[i]]$genes_blancos[p[[i]]$genes_blancos %in% unique(p[[i]]$mi_ontologia[which(p[[i]]$mi_ontologia[, 2] %in% terminos_sobrerepresentados[[i]][clus==k]), 1])]
  
            #modificación para que la intersección la calcule con intersección de conjuntos          
            #a<-c(a, p_interseccion_conjuntos(genes_a, genes_b, p[[i]]$genes_blancos))
            #modificación para que la intersección la calcule con un test de Fisher          
            a<-c(a, phyper(sum(genes_b%in%genes_a)-1L, length(genes_a), length(clusters==j) - length(genes_a), length(genes_b), lower.tail = F))
          }
        }
      
        a
        sum(a<0.05)/length(a)
        if(sum(a<0.05)/length(a) >= 0.8) {
          terminado = TRUE
        }
      }
    }
  }
  arbol[[length(arbol)]]$hermano   <- (length(arbol) - 1)
  arbol[[length(arbol)-1]]$hermano <- length(arbol)
  
  terminado
  if(!terminado) genes <- genes[clusters == i, ]
  dim(genes)genes
  table(clusters)
}
dim(genes)
clusters_finales[[14]] <- genes[clusters == i, ]

hits <- example1_mapped[sample(1:nrow(example1_mapped), nrow(genes[clusters == i, ])), "STRING_id"]
interacciones <- string_db$get_interactions(hits)

cat(nombres_de_genes[rownames(genes[clusters == 2, ])], sep="\n")

#g<-graph_from_edgelist(as.matrix(interacciones[interacciones$combined_score >=400, 1:2]), directed = FALSE)

interacciones <- as.matrix(interacciones[(apply(interacciones[, c(7:10, 12:13)] > 0, 1, sum) > 0 & interacciones$combined_score >=400), 1:2])
g<-graph_from_edgelist(interacciones, directed = FALSE)
plot(g, vertex.size=4, vertex.label = NA)
#table(components(g)$membership)
#plot(degree_distribution(g), log="xy")
gclus<-cluster_fast_greedy(g)
mod<-modularity(gclus)
mod

modularidad <- c()
for(i in 1:1000){
  hits <- example1_mapped[sample(1:nrow(example1_mapped), nrow(genes[clusters == i, ])), "STRING_id"]
  interacciones <- string_db$get_interactions(hits)
  interacciones <- as.matrix(interacciones[(apply(interacciones[, c(7:10, 12:13)] > 0, 1, sum) > 0 & interacciones$combined_score >=400), 1:2])
  g<-graph_from_edgelist(interacciones, directed = FALSE)
  #plot(g, vertex.size=4, vertex.label = NA)
  #table(components(g)$membership)
  #plot(degree_distribution(g), log="xy")
  gclus<-cluster_fast_greedy(g)
  modularidad <- c(modularidad, modularity(gclus))
}
shapiro.test(modularidad)
sum(modularidad>mod)/5000

cat(unlist(lapply(strsplit(V(g)[gclus$membership == 1]$name, "\\."), function(x) x[2])), sep="\n")

hits <- example1_mapped[example1_mapped$genes %in% nombres_de_genes[rownames(genes[clusters == 2, ])], "STRING_id"]
interacciones <- string_db$get_interactions(hits)

g<-graph_from_edgelist(as.matrix(string_interactions[, 1:2]), directed = F)
g
plot(g, vertex.size=4, vertex.label = NA)
table(components(g)$membership)
plot(degree_distribution(g), log="xy")
gclus<-cluster_walktrap(g)
modularity(gclus)

gr<-graph_from_edgelist(as.matrix(string_interactions_random[, 1:2]), directed = F)
gr
plot(gr, vertex.size=4, vertex.label = NA)
table(components(gr)$membership)
plot(degree_distribution(gr), log="xy")
grclus<-cluster_walktrap(gr)
modularity(grclus)

g<-graph_from_edgelist(as.matrix(string_interactions[apply(string_interactions[, c(7:10, 12:13)] > 0, 1, sum) > 0, 1:2]), directed = F)
g
plot(g, vertex.size=4, vertex.label = NA)
table(components(g)$membership)


cat(names(components(g)$membership[components(g)$membership == 1]), sep="\n")
cat(unlist(lapply(a, function(x) strsplit(x, "\\.")[[1]][2])), sep="\n")

b<-sample(nombres_de_genes[rownames(genes)], 311)
cat(b, sep="\n")

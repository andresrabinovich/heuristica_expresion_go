#seteo parametros HyperG
library(GO.db)
require(GOstats)
require(ath1121501.db)
library(igraph)
library(GOSemSim)
library(ggplot2)
library(grid)
library(gridExtra)

setwd("~/doctorado/programacion/test_de_sobrerepresentacion/")
source("librerias.R")
source("hypertest.R")
source("similaridad_semantica.R")
source("hipotesis_nula.R")
source("~/doctorado/programacion/clustering/clusterizador-paso-1.R")

chip     <- "ath1121501"
#universe <- unique(unlist(as.list(ath1121501ACCNUM)))
universe <- nombres_de_genes[rownames(genes)]

evidencia_a_filtrar = c("IEA", "NAS", "IC", "ND", "RCA")

p<-hypertest(nombres_de_genes[rownames(genes)][clusters == 1], universe, ontologia_a_filtrar = "BP", evidencia_a_filtrar)

ic <- -log2(p$anotaciones/p$anotaciones["GO:0008150"])

ancestros<-as.list(GOBPANCESTOR)
hijos<-as.list(GOBPOFFSPRING)

#terminos<-toTable(GOTERM)
#terminos<-lapply(terminos[p$dag@nodes], function(x) x@Term)

#Sin corrección
terminos_sobrerepresentados <- names(p$pvalues[p$pvalues < 0.05])
#Corrección BH
#categorias <- names(pvalues(hgt)[cbh(pvalues(hgt), 0.05)])
#Corrección B
#categorias <- names(pvalues(hgt)[pvalues(hgt)< 0.05/length(pvalues(hgt))])

#dyn.load("c/semsim.so")
#simRes2<-.Call("semsim", ancestros[p$terminos_go][1:5], ic)
simRes<-similaridad_semantica(ancestros[p$terminos_go][terminos_sobrerepresentados], ic)
#simRes2<-mgoSim(terminos_sobrerepresentados, terminos_sobrerepresentados, ont="BP", organism="arabidopsis", measure = "Resnik", combine=NULL)
distRes <- 1/(1+simRes)
#k <- kappa(p$mi_ontologia, names(which(p$pvalues < 0.05)))
#dk<- 0.5*(1-k)
#l<-list(semsim=distRes, dkappa=dk)
ldist<-list(semsim=distRes)

distancia <- ldist[[i]]
h<-hclust(as.dist(distancia[terminos_sobrerepresentados, terminos_sobrerepresentados]))
clus<-cutreeDynamic(h, distM = distancia[terminos_sobrerepresentados, terminos_sobrerepresentados], minClusterSize = 5)


mica<-names(which.max(ic[Reduce(intersect, ancestros[terminos_sobrerepresentados])]))
l<-list(unique(unlist(ancestros[terminos_sobrerepresentados])), hijos[mica][[1]])
clus_minimo<-unique(c(Reduce(intersect, l), mica, terminos_sobrerepresentados))

nAttrs<-list()
nAttrs$label[clus_minimo] <- ""
nAttrs$label[terminos_sobrerepresentados] <- clus
colores <- rainbow(max(clus)+1)#topo.colors(max(clus)+1)#
nAttrs$fillcolor[clus_minimo] <- "#FFFFFFFF"
names(nAttrs$fillcolor) <- clus_minimo
nAttrs$fillcolor[terminos_sobrerepresentados] <- colores[clus+1]
nAttrs$fontsize[clus_minimo]<-20
names(nAttrs$fontsize) <- clus_minimo

plot(subGraph(clus_minimo, p$dag), nodeAttrs = nAttrs)


#Primero agregamos todos los micas así calculamos una sola vez las dos matrices de distancia
mica <-vector(mode="list", length=2)
for(i in 1:length(l)){
  distancia <- ldist[[i]]
  h<-hclust(as.dist(distancia))
  clus<-cutreeDynamic(h, distM = distancia, minClusterSize = 5)
  for(j in sort(unique(clus))){
    mica[[i]]<-c(mica[[i]], names(which.max(ic[Reduce(intersect, ancestros[terminos_sobrerepresentados[clus==j]])])))
  }
}

simRes<-similaridad_semantica(ancestros[p$terminos_go][c(mica[[1]], terminos_sobrerepresentados)], ic)
#simRes<-mgoSim(terminos_sobrerepresentados, terminos_sobrerepresentados, ont="BP", organism="arabidopsis", measure = "Resnik", combine=NULL)
distRes <- 1/(1+simRes)
#k <- kappa(p$mi_ontologia, c(mica[[2]], names(which(p$pvalues < 0.05))))
#dk<- 0.5*(1-k)
ldist<-list(semsim=distRes)#list(semsim=distRes, dkappa=dk)

td<-matrix(ncol=11, nrow=0)
colnames(td) <- c("cluster", "terminos", "tipo_distancia", "distancia_media", "distancia_maxima", "ic_medio", "ic_mica", "media_armonica", "media_armonica_h0", "pv", "pvh0")
for(i in (1:length(ldist))[1]){
  distancia <- ldist[[i]]
  h<-hclust(as.dist(distancia[terminos_sobrerepresentados, terminos_sobrerepresentados]))
  clus<-cutreeDynamic(h, distM = distancia[terminos_sobrerepresentados, terminos_sobrerepresentados], minClusterSize = 5)
  g<-vector("list", length(unique(clus))) 
  for(j in sort(unique(clus))){
    mica<-names(which.max(ic[Reduce(intersect, ancestros[terminos_sobrerepresentados[clus==j]])]))
    l<-list(unique(unlist(ancestros[terminos_sobrerepresentados[clus==j]])), hijos[mica][[1]])
    clus_minimo<-unique(c(Reduce(intersect, l), mica, terminos_sobrerepresentados[clus==j]))
    clus_extendido<-unique(c(unlist(ancestros[terminos_sobrerepresentados[clus==j]]), terminos_sobrerepresentados[clus==j]))
    clus_extendido<-clus_extendido[-which(clus_extendido == "all")]
    #terminos_del_cluster <- unique(c(terminos_sobrerepresentados[clus==j], mica))
    terminos_del_cluster <- clus_minimo
    #terminos_del_cluster <- clus_extendido
    #terminos_del_cluster <- terminos_sobrerepresentados[clus==j]
    
    simRes<-similaridad_semantica(ancestros[p$terminos_go][c(mica[[1]], terminos_del_cluster)], ic)
    #simRes<-mgoSim(terminos_del_cluster, terminos_del_cluster, ont="BP", organism="arabidopsis", measure = "Resnik", combine=NULL)
    distancia <- 1/(1+simRes)
    
    dj<-distancia[terminos_del_cluster, terminos_del_cluster]
    dm<-mean(dj[upper.tri(dj)])
    dmax<-dj[which.max(dj)]
    icm<-mean(ic[terminos_del_cluster])
    icmin<-ic[mica]
    icmmin<-mean(ic[terminos_del_cluster])/ic[mica]
    media.harmonica<-exp(mean(log(-log(p$pvalues[terminos_del_cluster]))))
    media.harmonica_h0<-media_armonica_h0(p, terminos_del_cluster, "CC")
    s<-sort(ic[terminos_del_cluster])
    df<-as.data.frame(cbind(names(s), s, rbind(cbind(p$pvalues[names(s)], rep("pvalue", length(p$pvalues[names(s)]))), cbind(media.harmonica_h0$pvalues[names(s)], rep("Control Nulo", length(media.harmonica_h0$pvalues[names(s)]))))))
    rownames(df)<-1:nrow(df)
    colnames(df)<-c("Termino", "IC", "pv", "Tipo")
    df$IC<-as.numeric(as.character(df$IC))
    df$pv<-as.numeric(as.character(df$pv))
    g[[j]]<-ggplot(data=df, aes(x=IC, y = pv, group=Tipo, shape=Tipo, colour=Tipo)) +
      geom_line() +
      geom_point(size=5) +
      xlab("Contenido de información") +
      ylab("P-Value") +
      ggtitle(paste("Rama", j))
    
    print(g[[j]])
    td<-rbind(td, c(j, paste(as.character(terminos_del_cluster), collapse=", "), names(ldist[i]), dm, dmax, icm, icmin, media.harmonica, media.harmonica_h0$mh, paste(p$pvalues[terminos_del_cluster], collapse = ', '), paste(media.harmonica_h0$pvalues[terminos_del_cluster], collapse = ', ')))
  }
}
td<-data.frame(td)

grid.arrange(g[[1]], g[[2]], g[[3]], g[[4]], g[[5]], g[[6]], g[[7]], g[[8]], g[[9]], g[[10]], g[[11]], ncol = 3, top = "P-Value de los términos de las ramas en función del contenido de información (Resnik + CC)")
#----------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------




nAttrs<-list()
nAttrs$label <- ""#p$dag@nodes#terminos#
colores <- c("white", "black", "green", "red", "blue", "yellow", "orange")#rainbow(max(clus)+1)#topo.colors(max(clus)+1)#
#colores[1] <- "#000000FF"
nAttrs$fillcolor <- rep("#FFFFFFFF", length(p$dag@nodes))
names(nAttrs$fillcolor) <- nodes(p$dag)
nAttrs$fillcolor[terminos_sobrerepresentados] <- colores[clus+1]
nAttrs$label[terminos_sobrerepresentados] <- clus
nAttrs$fontsize<-rep(60, length(p$dag@nodes))
names(nAttrs$fontsize) <- names(nAttrs$fillcolor)
plot(p$dag, nodeAttrs = nAttrs)

mica<-names(which.max(ic[Reduce(intersect, ancestros[terminos_sobrerepresentados])]))
l<-list(unique(unlist(ancestros[terminos_sobrerepresentados])), hijos[mica][[1]])
clus_minimo<-unique(c(Reduce(intersect, l), mica, terminos_sobrerepresentados))
plot(subGraph(clus_minimo, p$dag), nodeAttrs = nAttrs)
clus_extendido<-unique(c(unlist(ancestros[terminos_sobrerepresentados]), terminos_sobrerepresentados))
clus_extendido<-clus_extendido[-which(clus_extendido == "all")]

#mica<-names(which.max(ic[Reduce(intersect, ancestros[terminos_sobrerepresentados])]))
#clus_extendido<-unique(c(unlist(ancestros[terminos_sobrerepresentados]), terminos_sobrerepresentados))
#clus_extendido<-clus_extendido[-which(clus_extendido == "all")]
#plot(subGraph(clus_extendido, hgt@goDag), nodeAttrs = nAttrs)

for(i in 1:length(ancestros)){
  ancestros[[i]]<-c(ancestros[[i]], names(ancestros[i]))
}

for(j in unique(clus)){
  j = 4
  mica<-names(which.max(ic[Reduce(intersect, ancestros[terminos_sobrerepresentados[clus==j]])]))
  l<-list(unique(unlist(ancestros[terminos_sobrerepresentados[clus==j]])), hijos[mica][[1]])
  clus_minimo<-unique(c(Reduce(intersect, l), mica, terminos_sobrerepresentados[clus==j]))
  plot(subGraph(clus_minimo, p$dag), nodeAttrs = nAttrs)
  clus_extendido<-unique(c(unlist(ancestros[terminos_sobrerepresentados[clus==j]]), terminos_sobrerepresentados[clus==j]))
  clus_extendido<-clus_extendido[-which(clus_extendido == "all")]
  #plot(subGraph(clus_extendido, p$dag), nodeAttrs = nAttrs)
}  

library(ggplot2)
s<-sort(ic[terminos_del_cluster])
df<-as.data.frame(cbind(s, t(rbind((p$pvalues[names(s)]), (media.harmonica_h0$pvalues[names(s)])))))
colnames(df)<-c("IC", "pv", "pvh0")
ggplot(data=df, aes(x=IC)) + geom_line(aes(y=pv), colour="red") + geom_point(aes(y=pv), colour="red") + geom_line(aes(y=pvh0), colour="green")  + geom_point(aes(y=pvh0), colour="green") + xlab("Contenido de información") +  ylab("P-Value") + ggtitle("P-Value de los términos de la rama 4 y control nulo en función del contenido de información")

  genes_cluster_go<-unique(ans[which(ans[, "go_id"] == "GO:0006950"), "gene_id"])#unique(subset(g2at, go_id %in% c(categorias[clus==i], mica), select=c(gene_id)))[, "gene_id"]
  genes_cluster<-nombres_de_genes[rownames(genes)][clusters == i]
  numW <- length(genes_cluster_go)
  numWdrawn <- length(intersect(genes_cluster_go, genes_cluster))
  numB <- length(universe) - numW
  numDrawn <- length(genes_cluster)
  p<-phyper(numWdrawn - 1L, numW, numB, numDrawn, lower.tail=FALSE)
  
  genes_cluster_go<-unique(subset(g2at, go_id %in% terminos_sobrerepresentados[clus==j], select=c(gene_id)))[, "gene_id"]
  numW <- length(genes_cluster_go)
  numWdrawn <- length(intersect(genes_cluster_go, genes_cluster))  
  psin<-(1-phyper(numWdrawn - 1L, numW, numB, numDrawn, lower.tail=TRUE))
  nAttrs$fillcolor[mica]<-"violet"
  plot(subGraph(clus_extendido, hgt@goDag), nodeAttrs = nAttrs, main=paste("p con mica=", p, " p sin mica=", psin, sep=""))
  nAttrs$fillcolor[mica]<-"#FFFFFFFF"
}


#llll<-sample(names(ancestros), 20)
#lllll<-llll

# Start the clock!
#ptm <- proc.time()

# Loop through the vector, adding one
#dyn.load("c/semsim_multithreading.so")
#simRes<-.Call("semsim", ancestros[llll], ic)

# Stop the clock
#proc.time() - ptm

# Start the clock!
ptm <- proc.time()

#Calculamos la similaridad semantica de resnik
dyn.load("c/semsim.so")
simRes2<-.Call("semsim", ancestros[llll], ic)
# Stop the clock
proc.time() - ptm
sum(simRes-simRes2)

#org.At.tairGO2ALLTAIRS@L2Rchain[[2]]@filter <- "evidence NOT IN ('IEA','NAS','IC','ND','RCA')"
#org.At.tairGO2ALLTAIRS@L2Rchain[[1]]@filter <- "Ontology = 'BP'"

#g2at<-toTable(org.At.tairGO2ALLTAIRS)

ontologias <- list()
ontologias[[1]] <- list(evidencia=c("IEA", "NAS", "IC", "ND", "RCA", "IEP"), nombre="bpa", tipo="BP")
ontologias[[2]] <- list(evidencia=c("IEA", "NAS", "IC", "ND", "RCA"), nombre="bpb", tipo="BP")
ontologias[[3]] <- list(evidencia=c("IEA", "NAS", "IC", "ND", "RCA", "IEP"), nombre="cc", tipo="CC")

cbh <- function(pvalues, alfa){
  pvalues<-sort(pvalues)
  m <- length(pvalues)
  cbh<-rep(FALSE, m)
  for(i in 1:m){
    cbh[i] <- (pvalues[i] <= i*alfa/m)
  }
  return(cbh)
}


for(ontologia in ontologias[2]){
  #evidencia <- ontologia$evidencia
  evidencia <- NULL

  params <- new("GOHyperGParams",
                geneIds=nombres_de_genes[rownames(genes)][clusters == 2], 
                universeGeneIds=universe,
                annotation=chip,
                ontology=ontologia$tipo,
                #ontology="BP",
                pvalueCutoff=1,
                conditional=FALSE,
                testDirection="over")

  hgt = tryCatch(
    {
      hyperGTest(params)
    }, 
    error = function(e) 
    {
      print(paste("MY_ERROR:  ",e))
      return(NA)
    }
  )
}

g2at<-toTable(org.At.tairGO2ALLTAIRS)
g2at<-unique(g2at[, c("gene_id", "go_id")])
g2at<-subset(g2at, gene_id %in% universe, select=c(gene_id, go_id))
anotaciones<-table(g2at[, "go_id"])

#Calculo el information content de todos, con GO:0008150 la cantidad de anotaciones que tiene el nodo raiz BP
ic <- -log2(anotaciones/anotaciones["GO:0008150"])
ic <- c(ic, 0)
names(ic)[length(ic)] <- "all"
library(GO.db)
an<-as.list(GOBPANCESTOR)


#Sin corrección
#categorias <- names(pvalues(hgt)[pvalues(hgt)< 0.05])
#Corrección BH
#categorias <- names(pvalues(hgt)[cbh(pvalues(hgt), 0.05)])
#Corrección B
#categorias <- names(pvalues(hgt)[pvalues(hgt)< 0.05/length(pvalues(hgt))])

categorias <- names(p)[p < 0.05]

#Calculamos la similaridad semantica de resnik
dyn.load("c/semsim.so")
simRes<-.Call("semsim", an[categorias], ic)

distRes <- 1/(1+simRes)
distRes <- 1/(1+m)
h<-hclust(as.dist(distRes))
clus<-cutreeDynamic(h, distM = distRes)
clus<-cutree(tree=h, h = 0.991)
#require(fpc)
#ch<-rep(0, 500)
#for(i in 2:500){
#  km<-kmeans(distRes, 5)
#  ch[i]<-calinhara(distRes, km$cluster)
#}
#plot(ch)
#clus<-km$cluster
#clus<-cutree(h, h = 0.27)

#library("Rgraphviz")
#dag<-igraph.from.graphNEL(hgt@goDag)
#dag<-delete.vertices(dag, names(pvalues(hgt)[pvalues(hgt) >= 0.05]))
nAttrs<-list()
nAttrs$label <- hgt@goDag@nodes
colores <- rainbow(max(clus)+1)
colores[1] <- "#000000FF"
nAttrs$fillcolor <- rep("#FFFFFFFF", length(hgt@goDag@nodes))
names(nAttrs$fillcolor) <- nodes(hgt@goDag)
nAttrs$fillcolor[terminos_sobrerepresentados] <- colores[clus+1]

plot(hgt@goDag, nodeAttrs = nAttrs)

for(j in unique(clus)){
  j = 2
  mica<-names(which.max(ic[Reduce(intersect, an[categorias[clus==j]])]))
  clus_extendido<-unique(c(unlist(an[categorias[clus==j]]), categorias[clus==j]))
  clus_extendido<-clus_extendido[-which(clus_extendido == "all")]
  
  genes_cluster_go<-unique(ans[which(ans[, "go_id"] == "GO:0006950"), "gene_id"])#unique(subset(g2at, go_id %in% c(categorias[clus==i], mica), select=c(gene_id)))[, "gene_id"]
  genes_cluster<-nombres_de_genes[rownames(genes)][clusters == i]
  numW <- length(genes_cluster_go)
  numWdrawn <- length(intersect(genes_cluster_go, genes_cluster))
  numB <- length(universe) - numW
  numDrawn <- length(genes_cluster)
  p<-phyper(numWdrawn - 1L, numW, numB, numDrawn, lower.tail=FALSE)
  
  genes_cluster_go<-unique(subset(g2at, go_id %in% categorias[clus==j], select=c(gene_id)))[, "gene_id"]
  numW <- length(genes_cluster_go)
  numWdrawn <- length(intersect(genes_cluster_go, genes_cluster))  
  psin<-(1-phyper(numWdrawn - 1L, numW, numB, numDrawn, lower.tail=TRUE))
  nAttrs$fillcolor[mica]<-"violet"
  plot(subGraph(clus_extendido, hgt@goDag), nodeAttrs = nAttrs, main=paste("p con mica=", p, " p sin mica=", psin, sep=""))
  nAttrs$fillcolor[mica]<-"#FFFFFFFF"
}


#-------------------------------------------0--------------------------------------------
m<-matrix(0, ncol=68, nrow=68)
rownames(m)<-colnames(m)<-categorias
for(a in 1:length(categorias)){
  for(b in 1:length(categorias)){
    m[a,b] <- goSim(categorias[a], categorias[b],  ont = "BP", organism = "arabidopsis", measure = "Resnik")
  }
}
#------------------------------------------0----------------------------------------
TABLENAME = "genes"
GENEIDS="gene_id"
SQL <- "SELECT DISTINCT go_id FROM %s INNER JOIN go_%s_all USING (_id) WHERE %s IN (%s) AND evidence NOT IN (%s)"
inClause2 <- toSQLStringSet(evidencia) # may get reused below
keep.all <- switch(testDirection(p),
                   over=FALSE,
                   under=TRUE,
                   stop("Bad testDirection slot"))
inClause1 <- if (!keep.all){
  geneIds(p)
}else{
  universe
}
inClause1 <- toSQLStringSet(inClause1) # may get reused below
SQL <- sprintf(SQL, TABLENAME, ontology(p), GENEIDS, inClause1, inClause2)
wantedGO <- dbGetQuery(db, SQL)[[1]]

SQL <- "SELECT DISTINCT %s, go_id FROM %s INNER JOIN go_%s_all USING (_id) WHERE %s IN (%s) AND go_id IN (%s) AND evidence NOT IN (%s)"
inClauseGO <- toSQLStringSet(wantedGO)
inClause1 <- toSQLStringSet(universe)
SQL <- sprintf(SQL, GENEIDS, TABLENAME, ontology(p), GENEIDS, inClause1, inClauseGO, inClause2)
ans <- dbGetQuery(db, SQL)
#----------------------------------------------0--------------------------------------
detach("package:GOstats", unload = TRUE)
detach("package:Category", unload = TRUE)
install.packages("~/doctorado/programacion/test_de_sobrerepresentacion/Category.tar.gz", type="source", repos = NULL)
#install.packages("~/doctorado/programacion/test_de_sobrerepresentacion/GOstats.tar.gz", type="source", repos = NULL)
require(GOstats)
hyperGTest(params)

phyper(1-1L, 3, 1627, 192, lower.tail=FALSE)

limite <- nrow(a)
V <- unique(c(a[, "go_idh"], a[, "go_idp"]))
edL <- vector("list", length=length(terminos_go))
names(edL) <- terminos_go
l<-aggregate(go_idp~go_idh, data=a, c)
rownames(l)<-l[, "go_idh"]
l <- l[terminos_go, ]
gr <- graphNEL(nodes=V, edgemode="directed")
for(i in rownames(l)){
  #edL[[i]] <- list(edges=l[i, "go_idp"][[1]])
  for(j in l[i, "go_idp"][[1]]){
    addEdge(from = i, to = j, gr)
  }
}
gR <- graphNEL(nodes=V, edgeL=edL, edgemode = "directed")
edges(gR)
edgeWeights(gR)

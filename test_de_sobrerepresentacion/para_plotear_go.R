cat(nombres_de_genes[rownames(clusters_finales[[1]])], sep="\n")

genes<-clusters_finales[[1]]

p                            <- hypertest(universo, ontologia = ontologia, ontologia_a_filtrar = "BP", evidencia_a_filtrar)
terminos_sobrerepresentadoss <- names(p$pvalues[cbh(p$pvalues, 0.05)])  
ic                           <- -log2(p$anotaciones/p$anotaciones["GO:0008150"])


terminos_sobrerepresentadoss <- terminos_sobrerepresentados[[i]]
icc <- ic[[i]]
ancestros<-as.list(GOBPANCESTOR)

mica<-names(which.max(icc[unique(c(Reduce(intersect, ancestros[terminos_sobrerepresentadoss]), terminos_sobrerepresentadoss))]))


distRes[[i]] <- 1/(1+similaridad_semantica(ancestros[p[[i]]$terminos_go][terminos_sobrerepresentados[[i]]], ic[[i]]))

simRes2<-mgoSim(terminos_sobrerepresentadoss, terminos_sobrerepresentadoss, ont="BP", organism="arabidopsis", measure = "Resnik", combine=NULL)
distRes[[i]] <- 1/(1+simRes2)

simRes2<-mgoSim(terminos_sobrerepresentadoss, terminos_sobrerepresentadoss, ont="BP", organism="arabidopsis", measure = "Lin", combine=NULL)
distRes[[i]] <- 1-simRes2

simRes2<-mgoSim(terminos_sobrerepresentadoss, terminos_sobrerepresentadoss, ont="BP", organism="arabidopsis", measure = "Rel", combine=NULL)
distRes[[i]] <- 1-simRes2

simRes2<-mgoSim(terminos_sobrerepresentadoss, terminos_sobrerepresentadoss, ont="BP", organism="arabidopsis", measure = "Jiang", combine=NULL)
distRes[[i]] <- 1-simRes2

simRes2<-mgoSim(terminos_sobrerepresentadoss, terminos_sobrerepresentadoss, ont="BP", organism="arabidopsis", measure = "Wang", combine=NULL)
distRes[[i]] <- 1-simRes2

h<-hclust(as.dist(distRes[[i]][terminos_sobrerepresentadoss, terminos_sobrerepresentadoss]))
plot(h)
#abline(h=0.3, col="red")
#clus<-cutree(h, h = 0.3)

clus<-cutreeDynamic(h, distM = distRes[[i]], minClusterSize = 3)
table(clus)


nAttrs<-list()
nAttrs$label <- unique(c(unlist(ancestros[terminos_sobrerepresentados]), terminos_sobrerepresentados))#p$dag@nodes#terminos#
nAttrs$label <- nAttrs$label[which(nAttrs$label != "all")]
nAttrs$fillcolor <- rep("#FFFFFFFF", length(nAttrs$label))
names(nAttrs$fillcolor) <- nAttrs$label
#nAttrs$fillcolor[terminos_sobrerepresentados] <- c("#FF0000FF", "#00FF00FF", "#0000FFFF", "#FF00FFFF")[clus]
nAttrs$fillcolor[terminos_sobrerepresentados] <- rainbow(length(table(clus)))[clus]
nAttrs$fontsize<-rep(60, length(nAttrs$label))
names(nAttrs$label) <- names(nAttrs$fontsize) <- names(nAttrs$fillcolor)
plot(subGraph(names(nAttrs$label), p$dag), nodeAttrs = nAttrs)

title(main = "Wang")


terminos_sobrerepresentados

p$mi_ontologia

genes_a<- p[[i]]$genes_blancos[p[[i]]$genes_blancos %in% unique(p[[i]]$mi_ontologia[which(p[[i]]$mi_ontologia[, 2] %in% terminos_sobrerepresentados[[i]][clus==j]), 1])]

histograma <- c()
for(g in universo){
  histograma <- c(histograma, max(ic[ontologia[which(ontologia[, "gene_id"] == g), "go_id"]]))
}
plot(ecdf(histograma))
abline(h=0.1, col="red")

hist(ic)
plot(ecdf(ic))



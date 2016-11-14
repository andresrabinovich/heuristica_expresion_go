genes_cluster = nombres_de_genes[rownames(genes)][clusters == 1]

#Los genes blancos son los genes de mi cluster que tienen anotaciones, es decir, los que están en la sub_ontologia
genes_blancos         <- p$genes_blancos

#Los genes negros son los del universo que no son blancos
#genes_negros              <- p$genes_negros

#phyper((sum(genes_cluster %in% genes_termino) - 1L), length(genes_blancos), length(genes_negros), length(genes_termino), lower.tail=FALSE)

n <- max(unique(clus))
vee<-vector("list", 15)
g<-vector("list", 2)
ll = 1
aaa<- c()
for(k in 1:(n-1)){
  g[[1]] = genes_blancos[genes_blancos %in% unique(p$mi_ontologia[which(p$mi_ontologia[, 2] %in% terminos_sobrerepresentados[clus==k]), 1])]
  aaa<-unique(c(aaa, g[[1]]))
  for(l in (k+1):n){
    names(g)<-c(k, l)
    g[[2]] = genes_blancos[genes_blancos %in% unique(p$mi_ontologia[which(p$mi_ontologia[, 2] %in% terminos_sobrerepresentados[clus==l]), 1])]
    #grid.newpage()
    #grid.draw(venn.diagram(g, filename = NULL))  
    vee[[ll]]<-gTree(children=venn.diagram(g, filename = NULL, main=paste(k, "(", colores[k], ")", "vs", l, "(", colores[l], ")")))
    ll = ll + 1
    aaa<-unique(c(aaa, g[[2]]))
  }
}
grid.newpage()
grid.arrange(vee[[1]], vee[[2]], vee[[3]],
             vee[[4]], vee[[5]], vee[[6]],
             vee[[7]], vee[[8]], vee[[9]],
             vee[[10]], vee[[11]], vee[[12]],
             vee[[13]], vee[[14]], vee[[15]],
             ncol=4)

p<-hypertest(genes_cluster[genes_cluster %in% unique(p$mi_ontologia[which(p$mi_ontologia[, 2] %in% terminos_sobrerepresentados[clus==6]), 1])], universe, ontologia_a_filtrar = "BP", evidencia_a_filtrar)

k = 5
cat(genes_cluster[genes_cluster %in% unique(p$mi_ontologia[which(p$mi_ontologia[, 2] %in% terminos_sobrerepresentados[clus==k]), 1])], sep="\n")
k = 6
cat(genes_cluster[genes_cluster %in% unique(p$mi_ontologia[which(p$mi_ontologia[, 2] %in% terminos_sobrerepresentados[clus==k]), 1])], sep="\n")
cat(genes_cluster, sep="\n")


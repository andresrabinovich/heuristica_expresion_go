media_armonica_h0<-function(p, terminos, ontologia = "BP"){
  g<-unique(p$mi_ontologia[, "gene_id"])
  a<-p$mi_ontologia[which(p$mi_ontologia[, "go_id"] %in% terminos), ]
  b<-unique(a[, "gene_id"])
  
  
  
  if(ontologia == "CC"){
    ontologia_a_filtrar <- "CC"
    #Guardamos los filtros anteriores
    evidencia_a_filtrar_anterior <- org.At.tairGO2ALLTAIRS@L2Rchain[[2]]@filter
    ontologia_a_filtrar_anterior <- org.At.tairGO2ALLTAIRS@L2Rchain[[1]]@filter
    
    #Aplicamos filtros si los hay
    org.At.tairGO2ALLTAIRS@L2Rchain[[2]]@filter <- "1"
    org.At.tairGO2ALLTAIRS@L2Rchain[[1]]@filter <- "1"
    if(is.vector(evidencia_a_filtrar)){
      org.At.tairGO2ALLTAIRS@L2Rchain[[2]]@filter <- paste("evidence NOT IN (", toSQLStringSet(evidencia_a_filtrar), ")", sep="")
    }
    if(!is.null(ontologia_a_filtrar)){ 
      if(ontologia_a_filtrar %in% c("BP", "CC", "MF")) {
        org.At.tairGO2ALLTAIRS@L2Rchain[[1]]@filter <- paste("Ontology = '", ontologia_a_filtrar, "'", sep="")
      }
    }
    print("FILTROS APLICADOS:")
    print(getBimapFilters(org.At.tairGO2ALLTAIRS))
    #-------------------------------FILTROS--------------------------------------
    
    
    #Traigo todos los genes y los go donde están anotados
    g2at                  <- as.list(org.At.tairGO2ALLTAIRS)
    
    #Volvemos a los filtros anteriores
    org.At.tairGO2ALLTAIRS@L2Rchain[[2]]@filter <- evidencia_a_filtrar_anterior
    org.At.tairGO2ALLTAIRS@L2Rchain[[1]]@filter <- ontologia_a_filtrar_anterior
    
    #g2at                  <- unique(g2at[, c("gene_id", "go_id")])
    
    ggg <- lapply(g2at, function(x){x[x %in% g]})
    gggg<-ggg[lapply(ggg, length)>1]
  }
  #a<-sort(p$anotaciones[terminos])
  #for(i in as.integer(names(table(a)[table(a)>1]))){
  #  names(a)[a==i]<-names(sort(unlist(lapply(ancestros[names(a)[a==i]], function(x) {
  #    sum(x %in% names(a)[a==i])
  #  })), decreasing = TRUE))
  #  n<-names(a)[a==i]
  #  for(j in 1:length(n)){
  #    for(k in j:length(n)){
  #      intersect(p$mi_ontologia[which(p$mi_ontologia[, "go_id"] %in% n[j]), ], p$mi_ontologia[which(p$mi_ontologia[, "go_id"] %in% n[k]), ])
  #    }
  #  }
  #}
  
  #g<-unique(p$mi_ontologia[, "gene_id"])
  pv2<-0
  iterac<-1000
  for(r in 1:iterac){
    if(ontologia == "CC") g <- gggg[[sample(which(lapply(gggg, function(x){length(x) >= length(b)}) == TRUE), 1)]]
    e<-a
    d<-sample(g, length(b), replace = FALSE)
    names(d)<-b
    for(i in 1:length(d)){
      e[which(a[, "gene_id"] == names(d)[i]), "gene_id"]<-d[i]
    }
    
    #l<-vector(mode="list", length=length(names(a)))
    #names(l)<-names(a)
    
    #for(i in 1:2){#length(l)){
    #  b<-sample(g[!(g %in% l[[i]])], (a[i]-length(l[[i]])), replace = FALSE);
    #  if(length(b) > 0){
    #    for(j in c(names(a[i]), intersect(names(a), ancestros[names(a)[i]][[1]]))){
    #      l[[j]]<-unique(c(l[[j]], b))
    #    }
    #  }
    #}
    #genes_blancos_sacados_en_cada_go <- unlist(lapply(l, function(x){length(intersect(x, p$genes_blancos))}))
    genes_blancos_sacados_en_cada_go <- unlist(lapply(terminos, function(x){length(intersect(unique(e[which(e[, "go_id"] == x), "gene_id"]), p$genes_blancos))}))
    genes_blancos <- p$genes_blancos
    genes_negros <- p$genes_negros 
    genes_sacados_en_cada_go <-  p$anotaciones[terminos]
    
    #Usamos phyper para hacer el test hypergeométrico. Hay que restarle 1 a los genes blancos sacados si usamos sobrerepresentación porque el test es para p<=p0 o p > p0
    pv2<-pv2+phyper((genes_blancos_sacados_en_cada_go - 1L), length(genes_blancos), length(genes_negros), genes_sacados_en_cada_go, lower.tail=FALSE)
  }
  pv2<-pv2/iterac
  names(pv2)<-terminos
  return (list(pvalues=pv2, mh=exp(mean(log(-log(pv2))))))
}
#p$pvalues[names(a)]
#pv2

#proc.time() - ptm

#anotaciones_random
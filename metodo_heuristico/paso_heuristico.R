paso_heuristico <- function(genes, ontologia, arbol, padre){

  distancia             <- as.dist(0.5*(1 - (cor(t(genes)))))
  clusters_jerarquicos  <- hclust(distancia, method = metodo)  
  clusters              <- cutree(tree = clusters_jerarquicos, h = (clusters_jerarquicos$height[length(clusters_jerarquicos$height)]+clusters_jerarquicos$height[length(clusters_jerarquicos$height)-1])/2)
  
  print(table(clusters))

  for(i in unique(clusters)){
    p                           <- hypertest(rownames(genes)[clusters == i], ontologia=ontologia)
    terminos_sobrerepresentados <- names(p$pvalues[cbh(p$pvalues, 0.05)])
    ic                          <- -log2(p$anotaciones/p$anotaciones["GO:0008150"])
    #descartamos los terminos sobrerepresentados con poco contenido de información
    terminos_sobrerepresentados <- terminos_sobrerepresentados[ic[terminos_sobrerepresentados] > 4]
    
    if(length(terminos_sobrerepresentados) > 1){
      #distRes <- 1/(1+similaridad_semantica(ancestros[p$terminos_go][terminos_sobrerepresentados], ic))
      bla <- terminos_sobrerepresentados
      m <- matrix(0, nrow = length(bla), ncol = length(bla))
      colnames(m) <- rownames(m) <- bla
      for (x in 1:length(bla)) {
        #print(bla[x])
        for(y in x:length(bla)){
          #cat(bla[x], bla[y], "\n", sep=" ")
          
          m[x, y] <- dishin(bla[x], bla[y], ic)
        }
      }
      m[lower.tri(m)] <- t(m)[lower.tri(m)]
      m      
      distRes <- 1/(1+m[terminos_sobrerepresentados, terminos_sobrerepresentados])
    }else{
      distRes <- NA
    }
    arbol[[length(arbol) + 1]] <- list(genes = rownames(genes)[clusters == i], terminos_sobrerepresentados = terminos_sobrerepresentados, padre = padre, terminado = FALSE)
    terminado <- FALSE
    
    if(TRUE){
    if(nrow(genes[clusters == i, ]) < 50){
      terminado <- TRUE
    }else{
      if(length(terminos_sobrerepresentados) < 2){
        if(length(terminos_sobrerepresentados) == 1){
          if(ic[terminos_sobrerepresentados] > quantile(ic, c(0.5))){
            terminado <- TRUE
          }
        }
      }else{
        distancia<-distRes
        h<-hclust(as.dist(distancia[terminos_sobrerepresentados, terminos_sobrerepresentados]))
        clus<-cutreeDynamic(h, distM = distancia[terminos_sobrerepresentados, terminos_sobrerepresentados], minClusterSize = 3)
        #clus<-cutree(tree = h, h = 0.3) 
        a<-c()
        lim <- max(as.integer(names(table(clus))))
        
        for(j in 1:(lim-1)){
          genes_a <- as.vector(t(unique(subset(ontologia, go_id %in% terminos_sobrerepresentados[clus == j] & gene_id %in% rownames(genes)[clusters == i], "gene_id"))))
          #genes_a<- p$genes_blancos[p$genes_blancos %in% unique(p$mi_ontologia[which(p$mi_ontologia[, 2] %in% terminos_sobrerepresentados[clus==j]), 1])]
          for(k in (j+1):lim){
            genes_b <- as.vector(t(unique(subset(ontologia, go_id %in% terminos_sobrerepresentados[clus == k] & gene_id %in% rownames(genes)[clusters == i], "gene_id"))))
            #genes_b<- p$genes_blancos[p$genes_blancos %in% unique(p$mi_ontologia[which(p$mi_ontologia[, 2] %in% terminos_sobrerepresentados[clus==k]), 1])]
            
            #modificación para que la intersección la calcule con intersección de conjuntos          
            #a<-c(a, p_interseccion_conjuntos(genes_a, genes_b, p$genes_blancos))
            #modificación para que la intersección la calcule con un test de Fisher          
            #a<-c(a, phyper((sum(genes_b %in% genes_a)-1L), length(genes_a), (sum(clusters==i) - length(genes_a)), length(genes_b), lower.tail = F))
            a<-c(a, phyper((sum(genes_b %in% genes_a)-1L), length(genes_a), (nrow(genes) - length(genes_a)), length(genes_b), lower.tail = F))
          }
        }
        
        a
        sum(a<0.05)/length(a)
        if(sum(a<0.05)/length(a) >= 0.8) {
          terminado <- TRUE
        }
      }
    }
    }else{
      if(is.null(nrow(genes[clusters == i, ])) || nrow(genes[clusters == i, ]) < 50){
        terminado = TRUE
      }
    }
    arbol[[length(arbol)]]$terminado <- terminado       
    if(!terminado){
      arbol <- paso_heuristico(genes[clusters == i, ], ontologia, arbol, (padre + 1))
    }
  }
  return(arbol)
  #arbol[[length(arbol)]]$hermano   <- (length(arbol) - 1)
  #arbol[[length(arbol)-1]]$hermano <- length(arbol)
}
ph<-paso_heuristico(genes, ontologia, arbol = list(), 1)

for(i in 1:length(ph)){
  if(ph[[i]]$terminado == TRUE){
    cat(i, length(ph[[i]]$genes), ph[[i]]$padre, ph[[i]]$terminado, "\n")
  }
}

cat(ph[[18]]$genes, sep="\n")

AT5G61380
AT2G25930
AT5G60100

AT5G23730
AT5G11260
AT5G62430
AT5G13930

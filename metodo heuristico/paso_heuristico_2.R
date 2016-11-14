combinar_conjuntos <- function(ontologia, terminos_sobrerepresentados, clus, genes){
  terminado = FALSE
  if(max(clus) > 0){
    while(terminado == FALSE){
      #lim <- max(as.integer(names(table(clus))))
      lim <- 1:(max(clus)-1)
      salir_del_for = FALSE
      for(j in 1:(max(clus)-1)){
        genes_a <- as.vector(t(unique(subset(ontologia, go_id %in% terminos_sobrerepresentados[clus == j] & gene_id %in% rownames(genes), "gene_id"))))
        if(length(genes_a) > 0){
          #genes_a
          #cat(genes_a, sep="\n")#genes_a<- p$genes_blancos[p$genes_blancos %in% unique(p$mi_ontologia[which(p$mi_ontologia[, 2] %in% terminos_sobrerepresentados[clus==j]), 1])]
          for(k in (j+1):max(clus)){
            genes_b <- as.vector(t(unique(subset(ontologia, go_id %in% terminos_sobrerepresentados[clus == k] & gene_id %in% rownames(genes), "gene_id"))))
            if(length(genes_b) > 0){
              genes_b
              sum(genes_a %in% genes_b)/max(length(genes_a), length(genes_b))
              a<-phyper((sum(genes_b %in% genes_a)-1L), length(genes_a), (140 - length(genes_a)), length(genes_b), lower.tail = F)
              #print(a)
              k = k + 1
              if(a < 0.05) { #Los dos grupos son iguales, los junto
                clus[clus == k] <- j
                #combinar_conjuntos(ontologia, terminos_sobrerepresentados, clus, genes, clusters)
                salir_del_for = TRUE
                break
              }
            }
          }
        }
        if(salir_del_for == TRUE) break      
      }
      terminado = TRUE
    }
  }
  #REVISAR QUE HACEMOS CON LOS GENES QUE CAYERON EN UN TERMINO QUE TIENE CLUS == 0
  #POR AHORA LOS METO EN LOS QUE QUEDARON SUELTOS
  #TAMBIEN PENSAR QUE PASA CON TODOS LOS GENES QUE DESCARTAMOS PORQUE NO ESTAN EN UN TERMINO SOBREREPRESENTADO. LOS JUNTO A LOS DEL CLUSTER "0"
  #POR AHORA Y SE SIGUEN PROCESANDO
  lim <- names(table(clus))
  if(sum(lim == "0") > 0){
    lim <- lim[-which(lim == "0")]
  }
  nuevos_conjuntos <- vector("list", (length(lim) + 1))
  if(length(lim) > 1){
    for(j in 1:length(lim)){
      nuevos_conjuntos[[j]] <- as.vector(t(unique(subset(ontologia, go_id %in% terminos_sobrerepresentados[clus == lim[j]] & gene_id %in% rownames(genes), "gene_id"))))
    }
  }
  nuevos_conjuntos[[(j+1)]] <- setdiff(rownames(genes), unique(unlist(nuevos_conjuntos)))
  if(FALSE){
    for(j in 1:length(lim)){
      a<-cor(t(genes))
      b<-nuevos_conjuntos[[j]]
      a<-a[b, b]
      print(mean(a[upper.tri(a)]))  
    }
    for(j in 1:(max(clus))){
    genes_a <- as.vector(t(unique(subset(ontologia, go_id %in% terminos_sobrerepresentados[clus == j] & gene_id %in% rownames(genes), "gene_id"))))
    cat(genes_a, sep="\n")#genes_a<- p$genes_blancos[p$genes_blancos %in% unique(p$mi_ontologia[which(p$mi_ontologia[, 2] %in% terminos_sobrerepresentados[clus==j]), 1])]
      
    a<-cor(t(genes))
    b<-genes_a
    a<-a[b, b]
    print(mean(a[upper.tri(a)]))   
    }
  }
  return(nuevos_conjuntos)
}

paso_heuristico_2 <- function(genes, ontologia, arbol, padre){
  
  distancia             <- as.dist(0.5*(1 - (cor(t(genes)))))
  clusters_jerarquicos  <- hclust(distancia, method = metodo)  
  #plot(clusters_jerarquicos)
  h <- (clusters_jerarquicos$height[length(clusters_jerarquicos$height)]+clusters_jerarquicos$height[length(clusters_jerarquicos$height)-1])/2
  #abline(h = h, col="red")
  clusters <- cutree(tree = clusters_jerarquicos, h = h)
  #abline(h=0.9, col="red")
  #clusters <- cutree(tree = clusters_jerarquicos, h = 0.9)
  #clusters <-cutreeDynamic(clusters_jerarquicos, distM = as.matrix(distancia), cutHeight = 0.75)#minClusterSize = 250)
  
  print(table(clusters))
  for(i in 1:max(clusters)){
    #if(mean(as.matrix(distancia)[rownames(genes[clusters == i, ]), rownames(genes[clusters == i, ])][upper.tri(as.matrix(distancia)[rownames(genes[clusters == i, ]), rownames(genes[clusters == i, ])])]) > 0.25){
    #  print(mean(as.matrix(distancia)[rownames(genes[clusters == i, ]), rownames(genes[clusters == i, ])][upper.tri(as.matrix(distancia)[rownames(genes[clusters == i, ]), rownames(genes[clusters == i, ])])]))
    #}
    #if(mean( cor(t(genes))[rownames(genes[clusters == i, ]), rownames(genes[clusters == i, ])][upper.tri( cor(t(genes))[rownames(genes[clusters == i, ]), rownames(genes[clusters == i, ])])]) < 0.5){
      a<-cor(t(genes))
      b<-rownames(genes[clusters == i, ])
      a<-a[b, b]
      print(mean(a[upper.tri(a)]))    
    #}    
  }
  
  for(i in unique(clusters)){
    p                           <- hypertest(rownames(genes)[clusters == i], ontologia=ontologia)
    terminos_sobrerepresentados <- names(p$pvalues[cbh(p$pvalues, 0.05)])
    ic                          <- -log2(p$anotaciones/p$anotaciones["GO:0008150"])
    #descartamos los terminos sobrerepresentados con poco contenido de información
    terminos_sobrerepresentados <- terminos_sobrerepresentados[ic[terminos_sobrerepresentados] > 4]
    print(length(terminos_sobrerepresentados))
    if(length(terminos_sobrerepresentados) > 1){
      distRes <- 1/(1+similaridad_semantica(ancestros[p$terminos_go][terminos_sobrerepresentados], ic))
    }else{
      distRes <- NA
    }
    arbol[[length(arbol) + 1]] <- list(genes = rownames(genes)[clusters == i], terminos_sobrerepresentados = terminos_sobrerepresentados, padre = padre, terminado = FALSE)
    terminado <- FALSE
    
    if(TRUE){
      if(nrow(genes[clusters == i, ]) < 15){
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
          plot(h)
          clus<-cutreeDynamic(h, distM = distancia[terminos_sobrerepresentados, terminos_sobrerepresentados], minClusterSize = 3)
          #clus <- cutree(tree=h, h=0.1)
          print(table(clus))
          #clus<-cutree(tree = h, h = 0.3)
          #Combinamos los grupos que salgan sobrerepresentados
          nuevos_conjuntos<-combinar_conjuntos(ontologia, terminos_sobrerepresentados, clus, genes[clusters == i, ])
        }
      }
    }else{
      if(is.null(nrow(genes[clusters == i, ])) || nrow(genes[clusters == i, ]) < 50){
        terminado = TRUE
      }
    }
    arbol[[length(arbol)]]$terminado <- terminado       
    if(!terminado){
      for(i in 1:length(nuevos_conjuntos)){
        arbol <- paso_heuristico_2(genes[nuevos_conjuntos[[i]], ], ontologia, arbol, (padre + 1))
      }
    }
  }
  return(arbol)
  #arbol[[length(arbol)]]$hermano   <- (length(arbol) - 1)
  #arbol[[length(arbol)-1]]$hermano <- length(arbol)
}
ph<-paso_heuristico_2(genes, ontologia, arbol = list(), 1)

for(i in 1:length(ph)){
  if(ph[[i]]$terminado == TRUE){
    cat(i, length(ph[[i]]$genes), ph[[i]]$padre, ph[[i]]$terminado, "\n")
  }
}

cat(ph[[18]]$genes, sep="\n")

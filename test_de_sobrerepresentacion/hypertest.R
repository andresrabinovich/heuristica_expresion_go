#Función para calcular la sobrerepresentación de genes de una lista en los términos de una ontología
#Recibe la lista de genes (parametro genes) y todo el universo de genes de background (universe)
#Puede recibir la ontología y un vector de tipos de evidencia que queremos sacar


hypertest <- function(genes, universe=NULL, ontologia=NULL, ontologia_a_filtrar=NULL, evidencia_a_filtrar=NULL){

  if(is.null(ontologia)){
    ontologia <- g2at(ontologia_a_filtrar, evidencia_a_filtrar) 
    ontologia <- unique(ontologia[, c("gene_id", "go_id")])
    
    #Me quedo únicamente con los genes de mi universo que tengan alguna anotacion y al conjunto de genes y terminos go los llamo mi ontologia
    ontologia             <- subset(ontologia, gene_id %in% universe, select=c(gene_id, go_id))
    universo              <- unique(ontologia[, "gene_id"])
  }else{
    universo              <- unique(ontologia[, "gene_id"])
  }    
  #Hago lo mismo que con el universo, pero solamente con los genes del cluster que quiero ver (los genes que serían mis bolas blancas) y lo llamo mi sub_ontologia
  sub_ontologia         <- subset(ontologia, gene_id %in% genes, select=c(gene_id, go_id))
  #Los términos go que me interesan son solamente los que están en mi sub_ontologia, es decir, los que tienen al menos uno de los genes del cluster anotados. Si no inmediatamente da 1 el pvalue  
  terminos_go           <- unique(sub_ontologia[, "go_id"])
  
  #Los genes blancos son los genes de mi cluster que tienen anotaciones, es decir, los que están en la sub_ontologia
  genes_blancos         <- unique(sub_ontologia[, "gene_id"])
  
  #Los genes negros son los del universo que no son blancos
  genes_negros          <- setdiff(universo, genes_blancos)
  
  #Calculo la cantidad de anotados que tiene cada termino de mi ontologia.
  #Recordar que la sub_ontologia tiene menos genes anotados en cada termino (tiene solo los blancos), mientras que la ontologia tiene todos los del universo 
  genes_sacados_en_cada_go  <- table(ontologia[, "go_id"])[terminos_go]
  genes_blancos_sacados_en_cada_go <- table(sub_ontologia[, "go_id"])[terminos_go]
  
  #Usamos phyper para hacer el test hypergeométrico. Hay que restarle 1 a los genes blancos sacados si usamos sobrerepresentación porque el test es para p<=p0 o p > p0
  p<-phyper((genes_blancos_sacados_en_cada_go - 1L), length(genes_blancos), length(genes_negros), genes_sacados_en_cada_go, lower.tail=FALSE)
  
  #Armamos el dag
  gobpc<-toTable(GOBPCHILDREN)
  names(gobpc)<-c("go_idp", "go_idc", "RelationshipType")
  gobpc <- subset(gobpc, go_idp %in% terminos_go, select=c(go_idc, go_idp))
  gr<-ftM2graphNEL(ft = as.matrix(gobpc[order(gobpc[,1]), ]), edgemode = "directed")
  names(gr@nodes) <- gr@nodes
  
  return(list(pvalues=sort(p), anotaciones=genes_sacados_en_cada_go, genes_blancos = genes_blancos, genes_negros=genes_negros, genes_blancos_sacados_en_cada_go=genes_blancos_sacados_en_cada_go, mi_ontologia=ontologia, terminos_go=terminos_go, dag=gr, g2at=g2at))
}

p_interseccion_conjuntos<-function(conjunto_a, conjunto_b, universo){
  n <- length(universo)
  a <- length(conjunto_a)
  b <- length(conjunto_b)
  l <- length(intersect(conjunto_a, conjunto_b))
  s <- 0
  for(i in l:min(a, b)){
    s <- s + choose(a,i)*choose((n-a), (b-i))
  }
  return(s/choose(n, b))
}

#Trae de GO2ALLTAIRS todos los genes filtrados por evidencia y ontología
g2at <- function(ontologia_a_filtrar=NULL, evidencia_a_filtrar=NULL){
  
  #-------------------------------FILTROS--------------------------------------
  
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
  g2at                  <- toTable(org.At.tairGO2ALLTAIRS)
  
  #Volvemos a los filtros anteriores
  org.At.tairGO2ALLTAIRS@L2Rchain[[2]]@filter <- evidencia_a_filtrar_anterior
  org.At.tairGO2ALLTAIRS@L2Rchain[[1]]@filter <- ontologia_a_filtrar_anterior
  
  return(g2at)
}

#N<-LETTERS[seq( from = 1, to = 10 )]
#A<-c("A", "C", "E")
#B<-c("B", "D", "F")
#p_interseccion_conjuntos(A, B, N)

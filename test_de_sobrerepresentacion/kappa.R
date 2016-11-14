#Función para calcular la distancia kappa entre términos de la ontología
kappa<-function(ontologia, terminos){

  m<-matrix(0, ncol=length(unique(ontologia[, "go_id"])), nrow=length(unique(ontologia[, "gene_id"])))
  rownames(m)<-unique(ontologia[, "gene_id"])
  colnames(m)<-unique(ontologia[, "go_id"])
  for(i in 1:nrow(ontologia)){
    m[ontologia[i,"gene_id"], ontologia[i,"go_id"]]<-1
  }
  m<-t(m)
  
  k <- matrix(0, nrow=length(terminos), ncol=length(terminos))
  rownames(k)<-colnames(k)<-terminos
  
  for(i in terminos){
    for(j in terminos){
      sisi<-sum(m[i, ] & m[j, ])
      nono<-sum(!m[i, ] & !m[j, ])
      sino<-sum(m[i, ] & !m[j, ])
      nosi<-sum(!m[i, ] & m[j, ])
      
      total <- sisi+nono+sino+nosi
      pra<-(sisi+nono)/total
      pre<-(sum(m[i, ])*sum(m[j, ])+sum(!m[i, ])*sum(!m[j, ]))/(total*total)
      k[i, j]<-(pra-pre)/(1-pre)
    }
  }
  return(k)
}

#Ejemplo de uso
#k <- kappa(p$mi_ontologia, names(which(p$pvalues < 0.05)))
#dk<- 0.5*(1-k)
#d<-as.dist(dk)
#h<-hclust(d)
#clus<-cutreeDynamic(h, distM = dk, minClusterSize = 5)

#save(ancestros, p, ic, file="1.RData")
#setwd("~/doctorado/programacion/test_de_sobrerepresentacion/")
#(load("1.RData"))
#l<-ancestros[p$terminos_go]
#n<-names(l)
#go <- sapply(n, function(x) {which(names(ic) %in% c(x, l[[x]]))})
#ic[length(ic)+1]<-0
#names(ic)[length(ic)]<-"all"
similaridad_semantica <- function(dag, ic){
  ic[length(ic)+1]<-0
  names(ic)[length(ic)]<-"all"
  l<-dag
  n<-names(l)
  m<-matrix(0, ncol=length(n), nrow=length(n))
  rownames(m)<-colnames(m)<-n
  nn<-length(n)
  for(i in 1:nn){
    for(j in i:nn){
      m[n[i], n[j]] <- max(ic[intersect(c(n[i], l[[n[i]]]), c(n[j], l[[n[j]]]))])
    }
  }
  m[lower.tri(m)]<-t(m)[lower.tri(m)]
  return(m)
}

#dyn.load("c/semsimv1.01.so")
#simRes<-.Call("semsim", go, ic)





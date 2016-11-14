# Implementa funcionalidad para setear filtros utilizados en objetos 
# AnnotationDbi (org.XX.XX)
# se necesita especificar 'org' y ec.out
# por ejemplo: org  <-"org.Sc.sgd"  
#             ec.out<- "IEA"

#funcion (interna) para ejecutar codigo generando script y sourceandolo
execRcode<-function(statment){
  cat(statment,sep="\n",file=".aux.R")
  source(".aux.R")
}

#reset filtro para los bimaps orgGO orgGO2ORF orgGO2ALLORFS
resetECfilter<-function(){
 statment<-c(
        paste(orgGO        ,"@L2Rchain[[2]]@filter <- \"1\"",sep=""),
        paste(orgGO2ORF    ,"@L2Rchain[[2]]@filter <- \"1\"",sep=""),
        paste(orgGO2ALLORFS,"@L2Rchain[[2]]@filter <- \"1\"",sep=""))
 execRcode(statment)     
}

#seteo filtro para los bimaps orgGO orgGO2ORF orgGO2ALLORFS
# ecOut: vector de codigos
#si ecOut="1" no se filtra nada
#       sino contiene los codigos a filtrar: por ejemplo ecOUT=c("IEA","NAS")    
setECfilter<-function(ecOut,verbose=FALSE){
 if(length(ecOut)==1){
  saux <- paste("evidence NOT IN ('",ecOut,"')",sep="")
 }else{
  saux<-c()
  for(i in seq_along(ecOut)){
    saux<-c(saux,paste("'",ecOut[i],"'",sep=""))
  }
  saux<-paste(saux,collapse=",")
  saux<-paste("evidence NOT IN (",saux,")",sep="") 
 } 

 statment<-c(
        paste(orgGO        ,"@L2Rchain[[2]]@filter <- \"",saux,"\"",sep=""),
        paste(orgGO2ORF    ,"@L2Rchain[[2]]@filter <- \"",saux,"\"",sep=""),
        paste(orgGO2ALLORFS,"@L2Rchain[[2]]@filter <- \"",saux,"\"",sep=""))
 execRcode(statment)   
 if(verbose)print(statment)
}

# a partir de un string org (e.g. org<-"org.Sc.sgd") se generan los objetos character 
# orgdb, orgGO, orgGO2ALLORFS, orgGO2ORF
initECfilter<-function(org,verbose=FALSE){
 orgdb<<-paste(org,".db",sep="")
 if(org=="org.Sc.sgd"){
  orgGO<<-paste(org,"GO",sep="")
  orgGO2ORF<<-paste(org,"GO2ORF",sep="")
  orgGO2ALLORFS<<-paste(org,"GO2ALLORFS",sep="")
 }else{
  if(org=="org.At.tair"){
   orgGO<<-paste(org,"GO",sep="")
   orgGO2ORF<<-paste(org,"GO2TAIR",sep="")
   orgGO2ALLORFS<<-paste(org,"GO2ALLTAIRS",sep="")
  }else{
   if(org=="ath1121501"){
    orgGO<<-paste(org,"GO",sep="")
    orgGO2ORF<<-paste(org,"GO2PROBE",sep="")
    orgGO2ALLORFS<<-paste(org,"GO2ALLPROBES",sep="")   
   }else{
    orgGO<<-paste(org,"GO",sep="")
    orgGO2ORF<<-paste(org,"GO2EG",sep="")
    orgGO2ALLORFS<<-paste(org,"GO2ALLEGS",sep="")
   } 
  }
 } 

 if(verbose) cat(paste(" orgGO <-",orgGO,"\n"),
                 paste("orgGO2ORF <-",orgGO2ORF,"\n"),
                 paste("orgGO2ALLORFS <-",orgGO2ALLORFS,"\n"))
} 


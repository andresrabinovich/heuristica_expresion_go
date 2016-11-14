#seteo parametros HyperG
library(GO.db)
require(GOstats)
require(ath1121501.db)
library(igraph)
library(GOSemSim)
library(ggplot2)
library(grid)
library(gridExtra)
#Cargamos librerías
library(genefilter) #Para el filtrado de genes K over A
library(dynamicTreeCut) #Para clusterizar con dynamic tree cut
#Traemos la base de datos del chipset que nos dice a que gen corresponde cada probeset
require(ath1121501.db)
library(VennDiagram)
library(STRINGdb)

### R code from vignette source 'GOvis.Rnw'

###################################################
### code chunk number 1: Setup
###################################################
library("Biobase")
library("annotate")
library("GOstats")
library("xtable")
library("multtest")
library("Rgraphviz")
library("hgu95av2.db")
library("GO.db")
library("genefilter")
library("ALL")


###################################################
### code chunk number 2: Example
###################################################
## subset of interest: 37+42 samples
data(ALL)
eset <- ALL[, intersect(grep("^B", as.character(ALL$BT)),
          which(as.character(ALL$mol) %in% c("BCR/ABL","NEG", "ALL1/AF4")))]

MB = factor(eset$mol.bio)

## intensities above 100 in at least 7 of the samples
f1 <- kOverA(7, log2(100))
f2 = function(x) {gps = split(2^x, MB); mns=sapply(gps, mean);
          if(max(mns) - min(mns) < 100) FALSE else TRUE}
ff <- filterfun(f1, f2)
selected <- genefilter(eset, ff)
sum(selected)
esetSub <- eset[selected,]



###################################################
### code chunk number 3: mtFstats
###################################################

 pvCutOff = 0.05

 Fstats = mt.maxT(exprs(esetSub), as.numeric(MB)-1,
                   B=1000, test="f")

 eSet = esetSub[Fstats$index[Fstats$adjp<pvCutOff],]

 gN = featureNames(eSet)
 lls = unlist(mget(gN, hgu95av2ENTREZID, ifnotfound=NA))

 eS = eSet[!duplicated(lls),]
 gN = featureNames(eS)
 lls = unlist(mget(gN, hgu95av2ENTREZID, ifnotfound=NA))
 lls = as.character(lls[!is.na(lls)])
 syms = unlist(mget(gN, hgu95av2SYMBOL, ifnotfound=NA))



###################################################
### code chunk number 4: inducedGO
###################################################
 gMF = makeGOGraph(lls, "MF", chip="hgu95av2")
 gBP = makeGOGraph(lls, "BP", chip="hgu95av2")
 gCC = makeGOGraph(lls, "CC", chip="hgu95av2")


 nMF = list()
 lbs = rep("", length(nodes(gMF)))
 names(lbs) = nodes(gMF)
 nMF$label = lbs

 nBP = list()
 lbs = rep("", length(nodes(gBP)))
 names(lbs) = nodes(gBP)
 nBP$label = lbs

 nCC = list()
 lbs = rep("", length(nodes(gCC)))
 names(lbs) = nodes(gCC)
 nCC$label = lbs

if( require("Rgraphviz") ) {
  gMFlo = agopen(gMF, "gMF")
  gBPlo = agopen(gBP, "gBP")
  gCClo = agopen(gCC, "gCC")
}



###################################################
### code chunk number 5: GOMF
###################################################
  if (require("Rgraphviz"))
     plot(gMFlo, nodeAttrs=nMF, main="Molecular Function")


###################################################
### code chunk number 6: GOBP
###################################################
  if (require("Rgraphviz"))
    plot(gBPlo, nodeAttrs=nBP, main="Biological Process")


###################################################
### code chunk number 7: GOCC
###################################################
  if (require("Rgraphviz"))
    plot(gCClo, nodeAttrs=nCC, main="Cellular Component")


###################################################
### code chunk number 8: findhigestmeans
###################################################
mns = apply(exprs(eS), 1, function(x) sapply(split(x, MB), mean))
whismax = apply(mns, 2, function(x) match(max(x), x))
maxNames = names(mns[,1])[whismax]
names(maxNames) = names(whismax)
table(maxNames)


###################################################
### code chunk number 9: genesToGO
###################################################

 CCnodes = nodes(gCC)
 nodes2affy = mget(CCnodes, hgu95av2GO2ALLPROBES, ifnotfound=NA)

 cts = sapply(nodes2affy, function(x) {
       wh = names(whismax) %in% x
       table(maxNames[wh])
   })

 all1cts = sapply(cts, function(x) x["ALL1/AF4"])
 all1cts = ifelse(is.na(all1cts), 0, all1cts)
 ##put nice names on
 names(all1cts) = names(cts)

 bcrcts = sapply(cts, function(x) x["BCR/ABL"])
 bcrcts = ifelse(is.na(bcrcts), 0, bcrcts)
 names(bcrcts) = names(cts)

 negcts = sapply(cts, function(x) x["NEG"])
 negcts = ifelse(is.na(negcts), 0, negcts)
 names(negcts) = names(cts)

 ctmat = cbind(all1cts, bcrcts, negcts)


###################################################
### code chunk number 10: layoutandrender
###################################################

if( require(Rgraphviz) ) {
opar = par(xpd = NA)
plotPieChart <- function(curPlot, counts, main) {
    renderNode <- function(x) {
        force(x)
        y <- x*100+1
        function(node, ur, attrs=list(), radConv=1) {
            nodeCenter <- getNodeCenter(node)
            pieGlyph(y, xpos=getX(nodeCenter),
                     ypos=getY(nodeCenter),
                     radius=getNodeRW(node),
                     col=c("blue", "green", "red"))
        }
    }
    drawing <- vector(mode="list", length=nrow(counts))
    for (i in 1:length(drawing)) {
        drawing[[i]] <- renderNode(counts[i,])
    }
    if( missing(main) )
       main="Example Pie Chart Plot"

    plot(curPlot, drawNode=drawing, main=main)

    legend(x="bottomleft", legend=c("ALL1/AF4", "BCR/ABL", "NEG"),
          fill=c("blue", "green", "red"))
}
plotPieChart(gCClo, ctmat)
par(opar)
}


###################################################
### code chunk number 11: imageMap
###################################################

 getNodeBB = function(ingraph) {
     xypos = getNodeXY(ingraph)
     hts = getNodeHeight(ingraph)
     wdsR = getNodeRW(ingraph)
     wdsL = getNodeLW(ingraph)
     llx = xypos$x - wdsL
     lly = xypos$y - (hts/2)
     urx = xypos$x + wdsR
     ury = xypos$y + (hts/2)
     cbind(llx, lly, urx, ury)
 }

 CCbb = getNodeBB(gCClo)

 ##now we need to adjust the graphBB to the plot region
 bbG = boundBox(gCClo)
 llx = getX(botLeft(bbG))
 lly = getY(botLeft(bbG))
 urx = getX(upRight(bbG))
 ury = getY(upRight(bbG))


##FIXME: work in progress here - we need to do a bit of mapping
##to get the right user coords on the jpg file...
if (interactive()) {
  jpeg("mygraph.jpg", quality=100, width=urx-llx, height=ury-lly)
  opar = par(plt=c(0,1,0,1), xpd=NA)
  plotPieChart(gCClo, ctmat)
  par(opar)
 dev.off()
}

## FIXME (wh 20.1.2005): better to use imageMap function for Ragraph objects
## in the package Rgraphviz instead.
if(require("geneplotter")) {
  con = openHtmlPage("example", "An Image Map")

  imageMap(CCbb, con, tags=list(
       TITLE = getNodeNames(gCClo),
       HREF  = rep("", length(nodes(gCC)))),
     imgname="mygraph.jpg")

  closeHtmlPage(con)
}



###################################################
### code chunk number 12: multip
###################################################

 igenes = unlist(mget(gN[1:10], hgu95av2ENTREZID))

 igenes = igenes[!is.na(igenes)]
 igenes = as.character(igenes)

 psets = mget(igenes, revmap(hgu95av2ENTREZID))
 psets = sapply(psets, function(x) x[1])
 GOs = mget(psets, hgu95av2GO, ifnotfound=NA)
 gBP = sapply(GOs, getOntology, "BP")
 ##how many terms

 sum(sapply(gBP, length))

 ##drop those with no BP annotations
 hasBP = sapply(gBP, function(x) length(x) > 0 && !is.na(x[1]))

 gBP = gBP[hasBP]

 ggs = lapply(igenes[hasBP], makeGOGraph, "BP", chip="hgu95av2.db")

## you could also do ggx = lapply(gBP, GOGraph, GOBPPARENTS)
## and should get the same thing

 simatM1 = matrix(1, nr=sum(hasBP), nc=sum(hasBP))
 for(i in 1:sum(hasBP))
    for( j in 1:sum(hasBP) )
        if( i== j ) next else
          simatM1[i,j] = simUI(ggs[[i]], ggs[[j]])

 library("RBGL")

 simatM2 = matrix(1, nr=sum(hasBP), nc=sum(hasBP))
 for(i in 1:sum(hasBP))
    for( j in 1:sum(hasBP) )
        if( i== j ) next else {
          simatM2[i,j] = simLP(ggs[[i]], ggs[[j]])
      }



###################################################
### code chunk number 13: <leavesEx
###################################################
gCC.leaves = leaves(gCC, "in")
length(gCC.leaves)
gCC.leaves[1:5]



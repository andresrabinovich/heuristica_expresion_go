source("http://bioconductor.org/biocLite.R")
biocLite("topGO")
biocLite("ALL")
biocLite("hgu95av2.db")
biocLite("multtest")
biocLite("lattice")
library(topGO)
library(ALL)
library(genefilter)
library(igraph)
library(hgu95av2.db, character.only = TRUE)
library(multtest)
library(lattice)
data(ALL)
data(geneList)
affyLib <- paste(annotation(ALL), "db", sep = ".")
sum(topDiffGenes(geneList))
sampleGOdata <- new("topGOdata",
                      description = "Simple session", ontology = "BP",
                      allGenes = geneList, geneSel = topDiffGenes,
                      nodeSize = 10,
                      annot = annFUN.db, affyLib = affyLib)
resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")

resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")

allRes <- GenTable(sampleGOdata, classicFisher = resultFisher,
                     classicKS = resultKS, elimKS = resultKS.elim,
                     orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 10)

pValue.classic <- score(resultKS)
pValue.elim <- score(resultKS.elim)[names(pValue.classic)]
gstat <- termStat(sampleGOdata, names(pValue.classic))
gSize <- gstat$Annotated / max(gstat$Annotated) * 4
gCol <- colMap(gstat$Significant)
colMap <- function(x) {
  .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
  return(.col[match(1:length(x), order(x))])
}
plot(pValue.classic, pValue.elim, xlab = "p-value classic", ylab = "p-value elim", pch = 19, cex = gSize, col = gCol)
sel.go <- names(pValue.classic)[pValue.elim < pValue.classic]
cbind(termStat(sampleGOdata, sel.go),
        elim = pValue.elim[sel.go],
        classic = pValue.classic[sel.go])

showSigOfNodes(sampleGOdata, score(resultKS.elim), firstSigNodes = 5, useInfo = ' all ' )

selProbes <- genefilter(ALL, filterfun(pOverA(0.20, log2(100)), function(x) (IQR(x) > 0.25)))
eset <- ALL[selProbes, ]


geneID2GO <- readMappings(file = system.file("examples/geneid2go.map", package = "topGO"))
str(head(geneID2GO))
geneNames <- names(geneID2GO)
head(geneNames)
myInterestingGenes <- sample(geneNames, length(geneNames) / 10)
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames
str(geneList)
GOdata <- new("topGOdata", ontology = "MF", allGenes = geneList,
                annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata
y <- as.integer(sapply(eset$BT, function(x) return(substr(x, 1, 1) == 'T')))
table(y)
geneList <- getPvalues(exprs(eset), classlabel = y, alternative = "greater")

topDiffGenes <- function(allScore) {
  return(allScore < 0.01)
}
x <- topDiffGenes(geneList)
sum(x) ## the number of selected genes

GOdata <- new("topGOdata",
                description = "GO analysis of ALL data; B-cell vs T-cell",
                ontology = "BP",
                allGenes = geneList,
                geneSel = topDiffGenes,
                annot = annFUN.db,
                nodeSize = 5,
                affyLib = affyLib)

allProb <- featureNames(ALL)
groupProb <- integer(length(allProb)) + 1
groupProb[allProb %in% genes(GOdata)] <- 0
groupProb[!selProbes] <- 2
groupProb <- factor(groupProb, labels = c("Used", "Not annotated", "Filtered"))
tt <- table(groupProb)
tt

pValue <- getPvalues(exprs(ALL), classlabel = y, alternative = "greater")
geneVar <- apply(exprs(ALL), 1, var)
dd <- data.frame(x = geneVar[allProb], y = log10(pValue[allProb]), groups = groupProb)
xyplot(y ~ x | groups, data = dd, groups = groups)

plot(density(a))
     
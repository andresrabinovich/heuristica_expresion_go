setMethod("categoryToEntrezBuilder",
          signature(p="KEGGHyperGParams"),
          function(p) {
              getKeggToEntrezMap(p)
          })

setMethod("categoryToEntrezBuilder",
          signature(p="GOHyperGParams"),
          function(p) {
              if (isDBDatPkg(p@datPkg))
                getGoToEntrezMap_db(p)
              else
                getGoToEntrezMap(p)
          })

setMethod("categoryToEntrezBuilder",
          signature(p="PFAMHyperGParams"),
          function(p) {
              getPfamToEntrezMap(p)
          })

getGoToEntrezMap_db <- function(p) {
    keep.all <- switch(testDirection(p),
                       over=FALSE,
                       under=TRUE,
                       stop("Bad testDirection slot"))
    #annPkgNS <- getNamespace(annotation(p))
    #db <- get("db_conn", annPkgNS)
    db <- do.call(paste(p@annotation, "dbconn", sep="_"), list())
    univ <- unlist(universeGeneIds(p), use.names=FALSE)

    ## For over representation:
    ## Obtain all unique GO IDs from specified ontology that have at
    ## least one of the genes from geneIds(p) annotated at it.
    ##
    ## For under representation:
    ## Obtain all unique GO IDs from specified ontology that have at
    ## least one of the genes from univ annotated at it.
    ##
    ## These are the GO IDs that form the keys in our GO_to_Entrez map.
    ## First we need to handle the fact that different species have different
    ## mappings for their names.
    if( is(p@datPkg, "YeastDatPkg") || is(p@datPkg, "Org.Sc.sgdDatPkg") ) {
       TABLENAME = "sgd"; GENEIDS="systematic_name"
    } else {
       TABLENAME = "genes"; GENEIDS="gene_id"
    }
    if(exists("evidencia")){
	cat("USANDO evidencia para el filtrado: ", str(evidencia));
    SQL <- "SELECT DISTINCT go_id
FROM %s INNER JOIN go_%s_all USING (_id)
WHERE %s IN (%s) AND evidence NOT IN (%s)"
    inClause2 <- toSQLStringSet(evidencia) # may get reused below
}else{
    SQL <- "SELECT DISTINCT go_id
FROM %s INNER JOIN go_%s_all USING (_id)
WHERE %s IN (%s)"
}
    inClause1 <- if (!keep.all)
      geneIds(p)
    else
      univ
    inClause1 <- toSQLStringSet(inClause1) # may get reused below
if(exists("evidencia")){
    SQL <- sprintf(SQL, TABLENAME, ontology(p), GENEIDS, inClause1, inClause2)
}else{
    SQL <- sprintf(SQL, TABLENAME, ontology(p), GENEIDS, inClause1)
}
    ##Check to make sure that we have some geneIds(p) to actually get GO terms FOR...
    ##basically, the "gene universe" at this point, is only those genes that
    ##have representation for on this ontology, and none of your genes were in
    ##THAT list, so there are not going to be any GO terms to find in this query.
    if(inClause1=="")
        stop("genes being tested do not have corresponding GO terms")
    wantedGO <- dbGetQuery(db, SQL)[[1]]
    ## Now collect the Entrez IDs annotated at our wantedGO IDs making
    ## sure to only keep those that are in the gene ID universe
    ## specified in p.
if(exists("evidencia")){
    SQL <- "SELECT DISTINCT %s, go_id
FROM %s INNER JOIN go_%s_all USING (_id)
WHERE %s IN (%s) AND go_id IN (%s) AND evidence NOT IN (%s)"
}else{
    SQL <- "SELECT DISTINCT %s, go_id
FROM %s INNER JOIN go_%s_all USING (_id)
WHERE %s IN (%s) AND go_id IN (%s)"
}
    inClauseGO <- toSQLStringSet(wantedGO)
    if (!keep.all)                      # avoid recomputing
      inClause1 <- toSQLStringSet(univ)
if(exists("evidencia")){    
SQL <- sprintf(SQL, GENEIDS, TABLENAME, ontology(p), GENEIDS, inClause1, 
                   inClauseGO, inClause2)
}else{
SQL <- sprintf(SQL, GENEIDS, TABLENAME, ontology(p), GENEIDS, inClause1, 
                   inClauseGO)
}
    ans <- dbGetQuery(db, SQL)
    if (nrow(ans) == 0)
        list()
    else 
        split(ans[[GENEIDS]], ans[["go_id"]])
}

getGoToEntrezMap <- function(p) {
    keep.all <- switch(testDirection(p),
                       over=FALSE,
                       under=TRUE,
                       stop("Bad testDirection slot"))
    lib <- p@datPkg
    ontology <- ontology(p)
    ## Return a list mapping GO ids to the Entrez Gene ids annotated
    ## at the GO id.  Only those GO ids that are in the specified
    ## ontology and have at least one annotation in the set of 
    ## Entrez Gene ids specified by 'selected' are included.
    go2allprobes <- GO2AllProbes(lib, ontology)
    probeAnnot <- getGoToProbeMap(go2allprobes, ontology)
    ## Map to Entrez Gene and flag GO ids that don't have any
    ## annotations in our selected set.  No sense testing these.
    probeToEntrezMapHelper(probeAnnot, geneIds(p), lib, universeGeneIds(p),
                           keep.all=keep.all)
}


getKeggToEntrezMap <- function(p) {
    keep.all <- switch(testDirection(p),
                       over=FALSE,
                       under=TRUE,
                       stop("Bad testDirection slot"))
##     lib <- annotation(p) ##not unless we want to dispatch on 'character'

##  Need to instead get the value from datPkg like this:
  lib <- p@datPkg
    
    ##kegg2allprobes <- revmap(getDataEnv("PATH", lib))  ## this means I should make a generic method for this step that will dispatch on the type and return something that can be "as.list()'ed"
    kegg2allprobes <- KEGG2AllProbes(lib)
    
    probeAnnot <- getKeggToProbeMap(kegg2allprobes)    
    probeToEntrezMapHelper(probeAnnot, geneIds(p), lib, universeGeneIds(p),
                           keep.all=keep.all)
}

.setKeytype <- function(p){
    keytype <- switch(class(p@datPkg),
                      YeastDatPkg = "ORF",
                      Org.XX.egDatPkg = "ENTREZID",
                      "PROBEID")
     keytype
}

## basically this function returns a list of PFAM IDs (names) where each element contains a vector of the entrez gene IDs that go with thos PFAM IDs.
getPfamToEntrezMap <- function(p) {
    keep.all <- switch(testDirection(p),
                       over=FALSE,
                       under=TRUE,
                       stop("Bad testDirection slot"))
    keytype <- .setKeytype(p)
    ##alteration:
    obj <- get(paste0(annotation(p),'.db'))
    probes <- keys(obj, keytype=keytype, column='PFAM')
    suppressWarnings(
       tab <- select(obj, keys=probes, columns='PFAM',keytype=keytype))
    probeAnnot <- split(tab[[keytype]], f=as.factor(tab$PFAM))
    probeToEntrezMapHelper(probeAnnot, geneIds(p), p@datPkg, universeGeneIds(p),
                           keep.all=keep.all)
}


probeToEntrezMapHelper <- function(probeAnnot, selected, lib, universe,
                                   keep.all=FALSE) {
    ## Given a list 'probeAnnot' mapping category => probe, convert the probe
    ## ids to Entrez Gene ids (unless we are using YEAST, in which case we skip
    ## this step).  Then reduce the Entrez Gene ids to the specified universe
    ## and only keep those entries in the list that have a non-empty
    ## intersection with the selected genes.  If keep.all is TRUE, then we keep
    ## entries even if the list of gene IDs includes no gene ID from the
    ## selected list.
    id2entrezEnv <- ID2EntrezID(lib)
    egAnnot <- lapply(probeAnnot, function(x) {
        z <- unique(unlist(mget(unique(x), id2entrezEnv)))
        z  <- intersect(z, universe)
        ## would be nice to have a short-circuiting way to do this
        if (length(z) > 0 && (keep.all || any(selected %in% z))) {
              return(z)
        }
        NULL
    })
    notNull <- sapply(egAnnot, function(x) !is.null(x))
    egAnnot[notNull]
}


getGoToProbeMap <- function(go2allprobes, ontology, goids) {
    ## Return a list with one element for each GO id in the specified 
    ## ontology that has at least one Probe probe id annotated at it.  
    ## The elements are vectors of Probe ids.  Names are GO ids.
    ##
    probeAnnot = as.list(go2allprobes)
    if (!missing(goids))
      probeAnnot = probeAnnot[goids]

    ## Remove "all" artifact GO nodes
    whAll = match("all", names(probeAnnot), 0)
    ## This is a tad faster than if (any(whAll > 0)) 
    ## because we stop early
    for (el in whAll) {
        if (el > 0) {
            probeAnnot = probeAnnot[-whAll]
            break
        }
    }

    goids = names(probeAnnot)

    ## filter on desired GO ontology
    inOnt <- filterGOByOntology(goids, ontology)
    probeAnnot = probeAnnot[inOnt]

    ## remove any GO ids that don't map to an probe id (NA)
    removeLengthZeroAndMissing(probeAnnot)
}


getKeggToProbeMap <- function(kegg2allprobes, keggIds) {
    probeAnnot = as.list(kegg2allprobes)
    if (!missing(keggIds))
      probeAnnot = probeAnnot[keggIds]
    removeLengthZeroAndMissing(probeAnnot)
}

getPfamToProbeMap <- function(pfam2allprobes, pfamIds) {
    probeAnnot = as.list(pfam2allprobes)
    if (!missing(pfamIds))
      probeAnnot = probeAnnot[pfamIds]
    removeLengthZeroAndMissing(probeAnnot)
}



removeLengthZeroAndMissing <- function(map) {
    notNA = sapply(map, function(x) {
        len = length(x)
        !(len == 0 || (len == 1 && is.na(x)))
    })
    map <- map[notNA]
}    

splitOrfByPfam <- function(ypfEnv){
  ## There should be a PFAM to ORF data set in addition to YEASTPFAM...a
  ## YEASTPFAM2PROBE, but for now we invert
  probe2pfam <- as.list(ypfEnv)
  probe2pfam <- removeLengthZeroAndMissing(probe2pfam)
  pfamlen <- listLen(probe2pfam)
  orf <- rep(names(probe2pfam), pfamlen)
  pfam <- unlist(probe2pfam)
  orfByPfam <- split(orf, pfam)
  orfByPfam

}

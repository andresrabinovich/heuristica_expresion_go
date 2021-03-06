\name{GSEAGOHyperGParams}
\alias{GSEAGOHyperGParams}
\alias{GSEAKEGGHyperGParams}

\title{Helper function for constructing a GOHyperGParams objects or
  KEGGHyperGParams objects from a GeneSetCollection}

\description{
  Helps to create A parameter class for representing all parameters
  needed for running the \code{hyperGTest} method.  If it is a
  GOHyperGParams object, being made, then with one of the GO ontologies
  (BP, CC, MF) as the category. This function will construct the
  parameter object from a GeneSetCollection object and if necessary will
  also try to check to make sure that the object is based on a GO2ALL
  mapping.
}

\usage{
  GSEAGOHyperGParams(name, geneSetCollection, geneIds, universeGeneIds,
  ontology, pvalueCutoff, conditional, testDirection, ...)

  GSEAKEGGHyperGParams(name, geneSetCollection, geneIds, universeGeneIds,
  pvalueCutoff, testDirection, ...)
}

\arguments{
    \item{name}{String specifying name of the GeneSetCollection.}
    \item{geneSetCollection}{A GeneSetCollection Object.  If a
      GOHyperGParams object is sought, then this GeneSetCollection
      should be based on a GO2ALLFrame object and so the idType of that
      GeneSetCollection should be GOAllFrameIdentifier. If a
      KEGGHyperGParams object is sought then a GeneSetCollection based
      on a KEGGFrame object should be used and the idType will be a
      KEGGFrameIdentifier.}
    \item{geneIds}{Object of class \code{"ANY"}: A vector of
      gene identifiers.  Numeric and character vectors are probably the
      only things that make sense.  These are the gene ids for the
      selected gene set.}
    \item{universeGeneIds}{Object of class \code{"ANY"}: A
      vector of gene ids in the same format as \code{geneIds} defining a
      subset of the gene ids on the chip that will be used as the
      universe for the hypergeometric calculation.  If this is
      \code{NULL} or has length zero, then all gene ids on the chip will
      be used.}
    \item{ontology}{A string specifying the GO ontology to use.
      Must be one of "BP", "CC", or "MF". (used with GO only)}
    \item{pvalueCutoff}{A numeric values between zero and one
      used as a p-value cutoff for p-values generated by the
      Hypergeometric test. When the test being performed is
      non-conditional, this is only used as a default value for
      printing and summarizing the results. For a conditional
      analysis, the cutoff is used during the computation to determine
      perform the conditioning: child terms with a p-value less than
      pvalueCutoff are conditioned out of the test for their parent
      term.}
    \item{conditional}{A logical indicating whether the
      calculation should condition on the GO structure. (GO only)}
    \item{testDirection}{A string which can be either "over"
      or "under". This determines whether the test performed detects
      over or under represented GO terms.}
    \item{...}{
      optional arguments to configure the GOHyperGParams object.
    }
    
}
  

\author{M. Carlson}

\seealso{
  \code{\link{HyperGResult-class}}
  \code{\link{GOHyperGParams-class}}
  \code{\link{hyperGTest}}
}

\keyword{classes}

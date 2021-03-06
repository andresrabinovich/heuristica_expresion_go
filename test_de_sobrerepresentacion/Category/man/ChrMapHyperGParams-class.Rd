\name{ChrMapHyperGParams-class}
\docType{class}
\alias{ChrMapHyperGParams-class}

\alias{conditional,ChrMapHyperGParams-method}
\alias{conditional<-,ChrMapHyperGParams,logical-method}

\title{Class "ChrMapHyperGParams"}
\description{
  This class encapsulates parameters needed for Hypergeometric testing
  of over or under representation of chromosome bands among a selected
  gene list using \code{\link{hyperGTest}}.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{new("ChrMapHyperGParams", ...)}.
}
\section{Slots}{
  \describe{

    \item{\code{chrGraph}:}{Object of class \code{"graph"}.  The nodes
      are the chromosome bands and the edges describe the tree structure
      of the bands.  Each node has a "geneIds" node attributes (see
      \code{nodeData}) which contains a vector of gene IDs annotated at
      the given band.}

    \item{\code{conditional}:}{Object of class \code{"logical"},
      indicating whether the test performed should be a conditional
      test.  }

    \item{\code{geneIds}:}{Object of class \code{"ANY"}: A vector of
      gene identifiers.  Numeric and character vectors are probably the
      only things that make sense.  These are the gene ids for the
      selected gene set.}

    \item{\code{universeGeneIds}:}{Object of class \code{"ANY"}: A
      vector of gene ids in the same format as \code{geneIds} defining a
      subset of the gene ids on the chip that will be used as the
      universe for the hypergeometric calculation.  If this is
      \code{NULL} or has length zero, then all gene ids on the chip will
      be used.}

    \item{\code{annotation}:}{A string giving the name of the
      annotation data package for the chip used to generate the data.}

    \item{\code{categorySubsetIds}:}{Object of class \code{"ANY"}:
      If the test method supports it, can be used to specify a subset of
      category ids to include in the test instead of all possible
      category ids.}

    \item{\code{categoryName}:}{A string describing the category.
      Usually set automatically by subclasses.  For example
      \code{"GO"}.}

    \item{\code{pvalueCutoff}:}{The p-value to use as a cutoff for
      significance for testing methods that require it.  This value will
      also be passed on to the result instance and used for display and
      counting of significant results.  The default is 0.01.}

    \item{\code{testDirection}:}{A string indicating whether the test
      should be for overrepresentation (\code{"over"}) or
      underrepresentation (\code{"under"}).}

    \item{\code{datPkg}:}{Object of class \code{"DatPkg"} used to assist
      with dispatch based on type of annotation data available. }

  }
}
\section{Extends}{
Class \code{"\linkS4class{HyperGParams}"}, directly.
}
\section{Methods}{
No methods defined with class "ChrMapHyperGParams" in the signature.
}

\author{Seth Falcon}

\examples{
showClass("ChrMapHyperGParams")
}
\keyword{classes}

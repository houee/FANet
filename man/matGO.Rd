\name{matGO}
\alias{matGO}
\docType{data}
\title{GO terms binary matrix
%%   ~~ data name/kind ... ~~
}
\description{Matrix of GO terms
Binary matrix (0/1) where rows corresponds to the 328 genes related to adiposity and columns to the 
5977 GO terms that are associated with at least one of these genes.
}
\usage{data(matGO)}
\format{
  A matrix with 328 genes and 5977 GO terms.
  }
%\source{UMR Genetique Animale - INRA/AGROCAMPUS OUEST - Rennes, France.}
\references{Blum Y., Le Mignon G., Lagarrigue S. and Causeur S. (2010) - A factor model to analyze heterogeneity in gene expression, BMC Bioinformatics, 11:368.
}
\examples{
data(matGO)
dim(matGO) 
}
%\keyword{datasets}

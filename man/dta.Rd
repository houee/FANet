\name{dta}
\alias{dta}
\docType{data}
\title{Gene expressions dataframe
%%   ~~ data name/kind ... ~~
}
\description{
This dataset concerns hepatic transciptome profiles of 45 half sib male chickens selected for their adipose variability. Rows correspond to the annotated genes having their expression correlated to adiposity (328 genes) and columns to the individuals (45 animals).
}
\usage{data(dta)}
\format{

  A dataframe with 328 genes and 45 individuals.
}
\source{UMR Genetique Animale - INRA/AGROCAMPUS OUEST - Rennes, France.}
\references{Blum Y., Le Mignon G., Lagarrigue S. and Causeur S. (2010) - A factor model to analyze heterogeneity in gene expression, BMC Bioinformatics, 11:368.
}
\examples{
data(dta)
dim(dta)
summary(dta)
}
\keyword{dataset}

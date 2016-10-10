\name{nbfactors}
\alias{nbfactors}
\title{nbfactors}
\description{Number of factors of the FA model}
\usage{nbfactors(dta, nb.iter=50, sv.sig=0.10)}
\arguments{
\item{dta}{A dataframe of the gene expression dataset where rows correspond
  to genes and columns to individuals.}
\item{nb.iter}{The number of null iterations to perform.}
\item{sv.sig}{The significance cutoff for eigengenes.}
}
\value{
\item{n.sv}{The optimal number of factors.}
}
\details{The function estimates the number of factors, based on the permutation test of Buja and Eyuboglu's method, and adapted from the package sva.}
\references{
Buja A and Eyuboglu N. (1992) Remarks on parallel analysis. Multivariate Behavioral Research, 27(4), 509-540
Leek JT and Storey JD. (2008) A general framework for multiple testing dependence. Proceedings of the National Academy of Sciences, 105: 18718-18723. http://www.biostat.jhsph.edu/~jleek/publications.html 

}
\author{David Causeur, Yuna Blum, Magalie Houee-Bigot.}
\keyword{methods}
\examples{
data(dta)
nbfopt = nbfactors(dta)
nbfopt 
}

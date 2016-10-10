\name{emfas.ccdlasso}
\alias{emfas.ccdlasso}
\title{Sparse Factor analysis model}
\description{A function that fits a cyclic coordinate descent algorithm for the LASSO estimation of the Factor model.}
\usage{emfas.ccdlasso(dta, nbf, lambda=0, zeros=NULL, min.err=1e-06, max.iter=100) }
\arguments{
\item{dta}{A dataframe of the gene expression dataset where rows correspond
  to genes and columns to individuals.}
\item{nbf}{Number of factors in the FA model, see \code{\link{nbfactors}}.}
\item{lambda}{Tuning parameter controlling the shrinkage (0 by default). See \code{\link{bicfa}} to choose the optimal value of lambda.}
\item{zeros}{Loadings equal to zero (NULL by default).}
\item{min.err}{Stopping criterion value for iterations in cyclic coordinate descent algorithm (default value: 1e-06).}
\item{max.iter}{Maximun number of iterations in the cyclic coordinate descent algorithm (default value: 100).}
}
\value{
\item{S}{Estimation of the correlation matrix.}	
\item{B}{Estimation of the loadings.}
\item{Psi}{Estimation of Psi, which corresponds to the specific variance component.}
}
\author{David Causeur, Yuna Blum, Magalie Houee-Bigot}
\keyword{methods}
\details{In order to use this function, the number of factors is requiered. You can use \code{\link{nbfactors}}.}
\examples{
data(dta)
nbfopt = nbfactors(dta) 
optimlambda = 0.46 # see the function bicfa
falasso = emfas.ccdlasso(dta, nbf=nbfopt, lambda=optimlambda)
}

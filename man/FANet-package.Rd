\name{FANet-package}
\alias{FANet-package}
\alias{FANet}
\docType{package}
\title{
Factor Analysis model for Network modeling (FANet)
}
\description{
The method proposed in this package takes advantage of a low-dimensional latent linear structure
of dependence to improve the stability of correlation estimations. An EM algorithm is used to fit 
a sparse factor model for correlations. Two functions are available in this package in order to 
introduce sparsity in the network estimation. The function \code{emfas.ccdlasso} is based 
on a  LASSO estimation using a cyclic coordinate descent algorithm. As an alternative, the second 
function \code{SparseEmfaBio} is based on biological knowledge integration as Gene Ontology annotation.
}

\details{
\tabular{ll}{
Package: \tab FANet\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-06-25\cr
License: \tab GPL>2\cr
}

}

\author{
Yuna Blum, Magalie Houee-Bigot, David Causeur.
Maintainer: Yuna Blum <yuna.blum[at]gmail.com>
}
% \references{
% ~~ Literature or other references for background information ~~
% }
% ~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
% ~~ the R documentation directory ~~
% \keyword{ package }
% \seealso{
% ~~ Optional links to other man pages, e.g. ~~
% ~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
% }
\examples{
data(dta)

# estimation of the correlation matrix using a factor model where sparsity is introduced using 
# LASSO regularization
data(dta)
nbf = 7 #see the function 'nbfactors'
lambda = 0.46 #see the function 'bicfa'
falasso = emfas.ccdlasso(dta,nbf=nbf,lambda=lambda,min.err=1e-01)
Slasso = falasso$S

# estimation of the correlation matrix using a factor model where sparsity is introduced using 
# biological knowledge
data(matGO) #see the function 'Construct.BioMatrix' to automatically retrieve GO annotations matrix from HGNC geneID.
zeros=MatBioPrep(matGO, optimal.nbclasses=FALSE,nbclasses=4)$zeros
Sbio=SparseEmfaBio(dta, zeros=zeros)
# see the function 'BioContribution' to evaluate the biological information contribution.
}


\keyword{Factor analysis, Co-expression networks, High dimension, Sparsity.}

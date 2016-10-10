\name{bicfa}
\alias{bicfa}
\title{bicfa}
\description{A function for choosing the optimal Tuning parameter.}
\usage{bicfa(dta,B)}
\arguments{
\item{dta}{A dataframe of gene expression where rows correspond
  to genes and columns to individuals.}
\item{B}{Matrix of loadings.}
}
\value{
\item{bic}{BIC criterium.}
}
\author{David Causeur, Yuna Blum, Magalie Houee-Bigot.}
\examples{
\dontrun{
data(dta)

nbfoptim = nbfactors(dta)
vlambda = seq(0.01,1,0.01)
vbic = rep(0,length(vlambda))
for (k in 1:length(vlambda)){
   falasso = emfas.ccdlasso(dta, nbf=nbfoptim, lambda=vlambda[k], min.err=1e-01)
   vbic[k] = bicfa(dta,falasso$B)$bic ; print(vbic)
}
optimlambda = vlambda[which.min(vbic)]
optimlambda
}
}
\keyword{methods}

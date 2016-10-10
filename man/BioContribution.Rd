\name{BioContribution}
\alias{BioContribution}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Contribution of \code{SparseEmfaBio} algorithm with biological constrains for gene module detection.
}
\description{
WARNING: this function should be considered experimental.
It compares empirically the \code{SparseEmfaBio} algorithm with biological contrains to approaches 
using biological information only or expression profiles only. These different approaches are compared in 
terms of gene module detection using 'WGCNA' method (see Horvath et al 2005 for more details).
}
\usage{
BioContribution(dta, Sbio, zeros, beta = 6, cutHeight = 0.998)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dta}{A dataframe corresponding to the gene expression dataset where rows correspond
  to genes and columns to individuals.
%%     ~~Describe \code{dta} here~~
}
  \item{Sbio }{
Correlation matrix calculated using the \code{SparseEmfaBio} function) with biological 
constrains. 
}
  \item{zeros }{
  Biological constrain matrix where row names correspond to gene identifiers and are exactly in the same order than in the expression dataset.
}
  \item{beta }{
Soft-thresholding power for network construction (see Horvath et al 2005 for more details).
 By default beta = 6. 
}
  \item{cutHeight}{Dendrogram cut height for module detection. See
          'cutreeDynamic' function in WGCNA package for more details.
%%     ~~Describe \code{cutHeight} here~~
}
}
%\details{
%%  
%}
\value{
Returns a list including:
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{BioContrib.Score}{Value which evaluates the contribution of biological knowledge 
and expression profiles in our algorithm. It corresponds to the ratio between the Rand 
Index comparing clustering using biological information only and \code{SparseEmfaBio} with 
biological constrains and Rand Index comparing clustering using \code{SparseEmfaBio} with 
biological constrains and using expression profiles only. If >1, the biological information has a greater impact on the network estimation than the gene expression itself and inversely.}
\item{Clust.OnlyBioinfo.PriorBioInfo }{Results comparing clustering using biological 
information only and \code{SparseEmfaBio} with biological constrains. It includes the contingency
 table couting the number of overlapping genes between modules, the corresponding -log10(P-value) using a Fisher's exact Test and the Rand Index.}
\item{Clust.WithoutBioinfo.PriorBioInfo }{Results comparing clustering using \code{SparseEmfaBio}
 with biological constrains and expression profiles only. It includes the contingency
 table couting the number of overlapping genes between modules, the corresponding -log10(P-value) using a Fisher's exact Test and the Rand Index.}
%% ...BioContrib.Score=BioContrib.Criteria, Clust.OnlyBioinfo.PriorBioInfo=res.GO.bio, Clust.WithoutBioinfo.PriorBioInfo=res.bio.without
The function returns a plot for both contingency tables where cells are colored according to the 
corresponding -log10(P-value) (Fisher's exact test).

}
\references{
Bin Zhang and Steve Horvath (2005) "A General Framework for
     Weighted Gene Co-Expression Network Analysis", Statistical
     Applications in Genetics and Molecular Biology: Vol. 4: No. 1,
     Article 17.
}
\author{
Yuna Blum, Magalie Houee-Bigot, David Causeur.
}
%\note{
%%  maybe a few words about beta....
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
 data(dta)
 data(matGO)
 
# Preprocess of the biological information matrix:
 Res.preparation=MatBioPrep(matGO, optimal.nbclasses=FALSE, nbclasses=4)
 
# Retrieve the constrain matrix:
 zeros=Res.preparation$zeros
 
# SparseEmfaBio algorithm with biological information to estimate the 
# correlations between gene expression:
 Sbio=SparseEmfaBio(dta, zeros=zeros)
 
# Contribution of SparseEmfaBio with biological constrains:
 Contri=BioContribution(dta, Sbio, zeros=zeros, beta=3)
 Contri$BioContrib.Score 
# The score here is close to 1: biological information and expression
# profiles seem to equally contribute to the correlation estimations using \code{SparseEmfaBio} algorithm.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

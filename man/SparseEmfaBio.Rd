\name{SparseEmfaBio}
\alias{SparseEmfaBio}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of the correlation matrix of gene expression using a sparse EM algorithm and
 biological constrains.
}
\description{
This function estimates the correlation profiles using a sparse factor analysis model. 
An EM alogithm is used to estimate the parameters of the model and sparsity is introduced using the 
biological constrain matrix in the M step. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SparseEmfaBio(dta, zeros, min.err = 0.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dta}{A dataframe of the gene expression dataset where rows correspond
  to genes and columns to individuals.
}
  \item{zeros}{
  Biological constrain matrix where row names correspond to gene identifiers and are exactly in the same order than in the expression dataset.
}
  \item{min.err}{Stopping criterion value for iterations (default value: 0.001).
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Returns the estimated correlation matrix.
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Yuna Blum, Magalie Houee-Bigot, David Causeur.
}
\note{
The expression dataset and the constrain matrix should have exactly the same row names which correspond to the gene identifiers.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data(dta)
data(matGO)
res = MatBioPrep(matGO, optimal.nbclasses=FALSE, nbclasses=4)
zeros = res$zeros
Sbio = SparseEmfaBio(dta, zeros=zeros)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

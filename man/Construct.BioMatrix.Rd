\name{Construct.BioMatrix}
\alias{Construct.BioMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates a matrix of biological information using GO annotation.
}
\description{
This function creates a matrix of biological information using GO annotation. For a given list of
 gene symbols (HGNC symbols) provided by the user, this function automatically retrieves the associated GO terms
 using Ensembl database and creates a binary matrix where rows correspond to genes and columns to 
 GO terms. Presence of 1 means that the gene is characterized by the corresponding GO term.
}
\usage{
Construct.BioMatrix(IDhgnc, GO.classe = "biological_process")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{IDhgnc}{A dataframe with two columns: the first one corresponds to gene identifiers that 
  are used in the expression dataset (in the same order) and the second column to the corresponding gene 
  symbol (HGNC symbol).
%%     ~~Describe \code{IDhgnc} here~~
}
  \item{GO.classe}{Desired class of GO annotation. Default is "biological_process". Other options
  are "molecular_function" or "cellular_component".
%%     ~~Describe \code{GO.classe} here~~
}
}
%\details{
%}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{matrixGO }{Binary matrix (0/1) where rows correspond to gene symbols (HGNC symbols) and
 columns to GO terms of the specified GO class. Presence of 1 means that the gene is characterized
 by the corresponding GO term.}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%BioMart?
%}
\author{
Yuna Blum, Magalie Houee-Bigot, David Causeur
}
\note{
This function uses BioMart package which requires an internet connection. Moreover, all gene
 symbols may not have available GO annotation: these genes should then be removed from the 
 expression dataset.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
## This function requieres an internet connection and may take a few minutes depending on the
# size of the genelist in the table 'IDhgnc'
## Don't run
# data(IDhgncFake)
# matGOfake = Construct.BioMatrix(IDhgncFake)
# dim(matGOfake)
# Among the 50 genes of this example, 45 have a GO annotation available. There are 321 GO terms 
# that are at least associated to one of the 50 genes.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

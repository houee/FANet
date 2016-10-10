\name{table2binmatrix}
\alias{table2binmatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transform a 2 columns table into a binary matrix. 
}
\description{
Transform a table with 2 columns corresponding respectively to gene symbols and G0 terms into 
a binary matrix (0/1)  where rows correspond to genes and columns to GO terms. Presence 
of 1 means that a gene is characterized by the corresponding GO term.
}
\usage{
table2binmatrix(tableGeneGO)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tableGeneGO}{A dataframe with 2 columns corresponding respectively to gene symbols 
  and G0 terms.
}
}
% \details{
%  ~~ If necessary, more details than the description above ~~
% }
\value{
Returns a binary matrix (0/1)  where rows correspond to gene symbols and columns 
to GO terms. Presence of 1 means that the gene is characterized by the corresponding GO term.
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
% \references{
% ~put references to the literature/web site here ~
% }
 \author{
Yuna Blum, Magalie Houee-Bigot, David Causeur.
 }
 \note{
This function can also be used for other kind of biological information in the same 2 columns 
format where the first column correponds to the gene identifiers.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

% \seealso{
% ~~objects to See Also as \code{\link{help}}, ~~~
% }
\examples{
# Toy example
gene.symbol = c(rep("gene1",3),rep("gene2",2),rep("gene3",1),rep("gene4",2),rep("gene5",2),
  rep("gene6",1),rep("gene7",1))
GOterm = c("GO1","GO2","GO3","GO1","GO3","GO2","GO1","GO2","GO1","GO2","GO2","GO2")
tableGeneGO = cbind.data.frame(gene.symbol,GOterm)
matrixGO = table2binmatrix(tableGeneGO)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

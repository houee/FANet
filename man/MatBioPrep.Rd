\name{MatBioPrep}
\alias{MatBioPrep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Preparation of the Biological Information Matrix
%%  ~~function to do ... ~~
}
\description{This function prepares the Biological Information Matrix before
 using it as the constrain matrix (called 'zeros' matrix) in the \code{SparseEmfaBio} algorithm. 
The process consists in filtering and clustering the columns corresponding to biological
information available for the genes. In the case where columns correspond to Gene Ontology (GO)
 terms, this function is useful to keep representative terms and deal with redundancy of
 the GO terms by clustering them.

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MatBioPrep(matrixBio, filtersize = TRUE, nbgenes.min = 50, nbgenes.max = 200, 
  optimal.nbclasses = TRUE, nbclasses.min = 5, nbclasses.max = 20, nbclasses = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matrixBio}{Binary matrix (0/1) where rows correspond to genes and columns to the 
  biological information as GO terms. Presence of 1 means that the gene is characterized by the
  biological information of the corresponding column. Gene identifiers (row names of the matrix)
  should be in the same order than in the expression dataset.
}
  \item{filtersize}{Logical. If TRUE, columns will be filtered according to 'nbgenes.min' and 
  'nbgenes.max' arguments.
}
  \item{nbgenes.min}{Only used if 'filtersize==TRUE', minimal number of genes associated to
  a biological term (column). 
}
  \item{nbgenes.max}{Only used if 'filtersize==TRUE', maximal number of genes associated to
  a biological term (column).
}
  \item{optimal.nbclasses}{Logical. If TRUE, optimal number of classes for the clustering will be
  determined using 'best.cutree' function.
}
  \item{nbclasses.min}{Only used if 'optimal.nbclasses==TRUE', minimum number of terms in one 
  cluster.
}
  \item{nbclasses.max}{Only used if 'optimal.nbclasses==TRUE', maximal number of terms in one 
  cluster.
}
  \item{nbclasses}{Only used if 'optimal.nbclasses==FALSE', an integer scalar with the desired 
  number of classes.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
 \item{zeros}{Constrain matrix that can be directly used in the \code{SparseEmfaBio} function.}
 \item{tree}{A tree of the biological terms clustering produced by \code{hclust}.}
 \item{bestcut }{Scalar corresponding to the optimal number of classes when 'best.cutree' function
 is used.}
 \item{MaxClassMembership }{Names of the biological terms the most correlated to each
 of the merged terms (columns of the 'zeros' matrix).}
}
% \references{
% ~put references to the literature/web site here ~
% }
\author{
Yuna Blum, Magalie Houee-Bigot, David Causeur.
}
\note{
The function \code{MatBioPrep} is particularly useful when the number of biological terms is high
which is typically the case with GO annotation.
}


% \seealso{
% ~~objects to See Also as \code{\link{help}}, ~~~
% }
\examples{
data(matGO)
Res.preparation=MatBioPrep(matGO, optimal.nbclasses=FALSE,  nbclasses=4)

# To vizualize tree of the biological terms. If you use 'optimal.nbclasses=FALSE', vizualize the 
# tree may help to define the desired number of classes:
plot(Res.preparation$tree)

# Preprocess biological information matrix where each column corresponds to a merged term:
zeros=Res.preparation$zeros

# Cluster biological interpretation: gives the biological term that is the most correlated to 
# each of the merged terms
Res.preparation$MaxClassMembership
# For example the GO term  "GO.0065007" is the GOterm the most correlated to the second merged term and 
# corresponds to 'biological regulation'. This information helps to interpret the biological meaning of 
# a merged term.
}
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

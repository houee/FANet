\name{IDhgncFake}
\alias{IDhgncFake}
\docType{data}
\title{Table of gene identifiers}
\description{Example of a dataframe with probe IDs and their corresponding gene symbols (HGNC).
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(IDhgncFake)}
\format{
  A dataframe of identifiers for 50 genes.
  \describe{
    \item{\code{gene_ID}}{probe identifiers.}
    \item{\code{HGNC_homosapiens}}{HGNC gene symbol.}
  }
}
% \details{
%  ~~ If necessary, more details than the __description__ above ~~
% }
% \source{
%  ~~ reference to a publication or URL from which the data were obtained ~~
% }
% \references{
%  ~~ possibly secondary sources and usages ~~
% }
\examples{
data(IDhgncFake)
head(IDhgncFake)
}
%\keyword{datasets}

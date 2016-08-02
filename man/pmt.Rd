\name{pmt}
\alias{pmt}
\title{Poisson Margin Test}
\description{
  Computes p-values for two sets of counts using the Poisson Margin Test.
}
\usage{
pmt(a, b)
}
\arguments{
  \item{a}{
    Vector of counts
}
  \item{b}{
    Vector of counts
}
}
\details{
  This function computes the Poisson Margin Test (PMT) for two sets of counts.
}
\value{
  A vector of log10 p-values is returned, corresponding to the input vectors.
}
\references{
Kowalczyk, A., Bedő, J., Conway, T., Beresford-Smith, B., 2011. The Poisson margin test for normalization-free significance analysis of NGS data. Journal of Computational Biology 18, 391–400. doi:10.1089/cmb.2010.0272
}
\author{
  Justin Bedő
}
%\examples{
%}
\keyword{NGS}
\keyword{Poisson}

\name{bnp4c-term}
\alias{bnp4c}
\alias{b1np4c}
\alias{b2np4c}
\docType{methods}
\title{Bipartite nodepower four-cycles terms (b1np4c and b2np4c)}
\description{
  Nodewise four-cycle counts raised to a power (inside-\eqn{\alpha} weighting) for first and second modes of bipartite networks.
}
\usage{
# binary: b1np4c(alpha=0.5, fixed=TRUE)
# binary: b2np4c(alpha=0.5, fixed=TRUE)
}
\arguments{
\item{alpha}{value of \eqn{\alpha}, the power to which the four-cycle count at each node is raised.}
\item{fixed}{optional argument indicating whether the \eqn{\alpha} parameter is fixed at the given value, or is to be fit as a curved epxonential-family model (see Hunter and Handcock, 2006). The default value is TRUE, which means that the \eqn{\alpha} parameter is fixed. Currently this parameter must be set to TRUE.}
} 

\details{

TODO: PUT DESCRIPTION HERE.
The argument fixed indicates whether the scale parameter lambda is to be fit as a curved exponential-family model (see Hunter and Handcock, 2006). Currently, the fixed argument must be set to TRUE, which means the weight parameter alpha is not separately estimated and thus the model is not a CEF model. This term can only be used with undirected bipartite networks.

\code{b2np4c(alpha=0.5 fixed=TRUE) (binary) (bipartite) (undirected)}


TODO: PUT DESCRIPTION HERE.

Note that allthough the default is \code{fixed=FALSE} the terms are not yet able to handle a non-fixed decay term so it must be set to TRUE.
}

\seealso{
  See also \code{\link[ergm]{'gwb1dsp-ergmTerm'}} term. 
}


%\value{\code{\link{ergm}} returns an object of class 'ergm' that is a list.
%}
\references{
Hunter, D. R. and M. S. Handcock (2006). Inference in curved exponential family models for networks. Journal of Computational and Graphical Statistics, 15: 565-583.

Stivala, A., Wang., P., and Lomi, A. (2024). Improving exponential-family random graph models for bipartite networks. Unpublished manusript.

Wang, P., Sharpe, K., Robins, G. L., and Pattison, P. E. (2009). Exponential random graph (p*) models for affiliation networks. Social Networks. 31(1): 12-25.

}

\author{
Alex Stivala \email{alex.d.stivala@gmail.com}
}

\examples{
require(ergm.terms.contrib)
# TODO: PUT EXAMPLE HERE

}
\keyword{models}

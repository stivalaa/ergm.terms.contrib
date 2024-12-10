\name{bnp4c-term}
\alias{bnp4c}
\alias{b1np4c}
\alias{b2np4c}
\docType{methods}
\title{Bipartite nodepower four-cycles terms (b1np4c and b2np4c)}
\description{
  Nodewise four-cycle counts raised to a power (\eqn{\alpha}-inside weighting) for first and second modes of bipartite networks.
}
\usage{
# binary: b1np4c(alpha=0.5, fixed=TRUE)
# binary: b2np4c(alpha=0.5, fixed=TRUE)
}
\arguments{
\item{alpha}{value of \eqn{\alpha}, the power to which the four-cycle count at each node is raised, where \eqn{0 < \alpha \leq 1}.}
\item{fixed}{optional argument indicating whether the \eqn{\alpha} parameter is fixed at the given value, or is to be fit as a curved exponential-family model (see Hunter and Handcock, 2006). The default value is TRUE, which means that the \eqn{\alpha} parameter is fixed. Currently this parameter must be set to TRUE.}
} 

\details{

\code{b1np4c(alpha=0.5 fixed=TRUE) (binary) (bipartite) (undirected)}

\code{b2np4c(alpha=0.5 fixed=TRUE) (binary) (bipartite) (undirected)}


These statistics are described in Stivala et al. (2024).
The \code{b1np4c} and \code{b2np4c} statistics were first implemented in EstimNetDirected (\url{https://github.com/stivalaa/EstimNetDirected}) as BipartiteFourCyclesNodePowerA and BipartiteFourCyclesNodePowerB.

These statistics sum the number of four-cycles at each node in one
bipartition (the first bipartition for \code{b1np4c} and the second
bipartition for \code{b2np4c}) raised to the power \eqn{\alpha} (\eqn{0 <
  \alpha \leq 1}). That is, the count at each node is raised to the power
\eqn{\alpha}, and then they are summed. In the terminology of
Wilson et al. (2017), this is an "\eqn{\alpha}-inside" weighting.

The argument \code{fixed} indicates whether the parameter \code{alpha}
is to be fit as a curved exponential-family model (see Hunter and
Handcock, 2006). Currently, the fixed argument must be set to TRUE,
which means the weight parameter alpha is not separately estimated and
thus the model is not a CEF model.

This term can only be used with undirected bipartite networks.

Note that allthough there is a \code{fixed} paramter, the terms are not yet able to handle a non-fixed \code{alpha} term so it must be set to the default value \code{fixed=TRUE}.
}

\seealso{
  See also \code{\link[ergm]{gwb1dsp-ergmTerm}} term. 
}


\references{
Hunter, D. R. and M. S. Handcock (2006). Inference in curved exponential family models for networks. Journal of Computational and Graphical Statistics, 15: 565-583.

Stivala, A., Wang., P., and Lomi, A. (2024). Improving exponential-family random graph models for bipartite networks. Unpublished manusript.

Wang, P., Sharpe, K., Robins, G. L., and Pattison, P. E. (2009). Exponential random graph (p*) models for affiliation networks. Social Networks. 31(1): 12-25.

Wilson, J. D., Denny, M. J., Bhamidi, S., Cranmer, S. J., & Desmarais,
B. A. (2017). Stochastic weighted graphs: Flexible model specification and simulation. Social Networks, 49, 37-47.
}

\author{
Alex Stivala \email{alex.d.stivala@gmail.com}
}

\examples{
require(ergm.terms.contrib)
# TODO: PUT EXAMPLE HERE

}
\keyword{models}

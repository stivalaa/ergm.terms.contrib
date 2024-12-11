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
which means the weight parameter \code{alpha} is not separately estimated and
thus the model is not a CEF model.

This term can only be used with undirected bipartite networks.

Note that although there is a \code{fixed} parameter, the terms are not yet able to handle a non-fixed \code{alpha} term so it must be set to the default value \code{fixed=TRUE}.
}

\seealso{
  See also \code{\link[ergm]{gwb1dsp-ergmTerm}} term. 
}


\references{
Hunter, D. R. and M. S. Handcock (2006). Inference in curved exponential family models for networks. Journal of Computational and Graphical Statistics, 15: 565-583.

Stivala, A., Wang., P., and Lomi, A. (2024). Improving exponential-family random graph models for bipartite networks. Unpublished manuscript.

Wang, P., Sharpe, K., Robins, G. L., and Pattison, P. E. (2009). Exponential random graph (p*) models for affiliation networks. Social Networks. 31(1): 12-25.

Wilson, J. D., Denny, M. J., Bhamidi, S., Cranmer, S. J., & Desmarais,
B. A. (2017). Stochastic weighted graphs: Flexible model specification and simulation. Social Networks, 49, 37-47.
}

\author{
Alex Stivala \email{alex.d.stivala@gmail.com}
}

\examples{
library(ergm.terms.contrib)

## Construct the Four-fan-3 graph
fourfan.3.df <- data.frame(A = c(1, 1, 1, 1, 1, 1, 3, 3, 6, 6, 9, 9),
                           B = c(2, 4, 5, 7, 8, 10,2, 4, 5, 7, 8, 10))
fourfan.3.net <- network(fourfan.3.df, bipartite=TRUE, directed=FALSE)

## From Stivala et al. (2024):
##
##   ... in this graph, the nodes in mode B contribute more to the
##   total as each one (of the six) is involved in exactly one
##   four-cycle (and hence raising to the power of α still
##   contributes one to the sum), while of the four nodes in mode A,
##   three are involved in only one four-cycle, while the fourth is
##   involved in three four-cycles and hence contributes only 3^α ≈
##   1.73205 (when α = 0.5).
##
summary(fourfan.3.net ~ b1np4c(0.5, TRUE) + b2np4c(0.5, TRUE))

## generate some graphs similar to 'zero.pos' in Fig. 5 and Fig. 6
## of Stivala et al. (2024), with the parameter for b1np4c zero
## and the parameter for b2np4c positive
## (warning: takes about 30 seconds to run)
system.time( g.zero.pos <- simulate(network(150, bipartite=100,
                                           directed=FALSE) ~
  edges + gwb1degree(.5,TRUE) + gwb2degree(.5, TRUE) + b2np4c(1/5),
  coef=c(-6.0, 0.6, -0.5, 3.0),
  nsim=100,
  control = control.simulate(MCMC.burnin = 100000,
                             MCMC.interval=10000)) )

## plot graph statistics to check sufficient burnin and interval
par(mfrow=c(2,3))
plot(summary(g.zero.pos ~ edges))
plot(summary(g.zero.pos ~ gwb1degree(.5,TRUE)))
plot(summary(g.zero.pos ~ gwb2degree(.5,TRUE)))
plot(summary(g.zero.pos ~ b2np4c(1/5)))
plot(summary(g.zero.pos ~ b1np4c(1/5)))
plot(summary(g.zero.pos ~ cycle(4)))

## plot visualization of one graph similar to
## Fig. 6 of Stivala et al. (2024)
plot(g.zero.pos[[88]],
     vertex.col= c(rep('red',100), rep('blue',50)),
     vertex.shape = c(rep(50, 100), rep(4,50)))
}
\keyword{models}

#  This software is distributed under the GPL-3 license.  

#' @templateVar name b1np4c
#'
#' @title Bipartite nodepower four-cycles terms (b1np4c and b2np4c)
#'
#' @description Nodewise four-cycle counts raised to a power
#'   (\eqn{\alpha}-inside weighting) for first and second modes of
#'   bipartite networks.
#'
#' @usage
#' #binary: b1np4c(alpha=0.5, fixed=TRUE)
#'
#' @template ergmTerm-general
#'
#' @param alpha value of \eqn{\alpha}, the power to which the
#'   four-cycle count at each node is raised, where \eqn{0 < \alpha
#'   \leq 1}.
#'
#' @param fixed optional argument indicating whether the \eqn{\alpha}
#' parameter is fixed at the given value, or is to be fit as a curved
#' exponential-family model (see Hunter and Handcock, 2006). The
#' default value is TRUE, which means that the \eqn{\alpha} parameter
#' is fixed. Currently this parameter must be set to TRUE.
#'
#' @details
#'
#' These statistics are described in Stivala et al. (2025).  The
#' \code{b1np4c} and \code{b2np4c} statistics were first implemented
#' in EstimNetDirected
#' (\url{https://github.com/stivalaa/EstimNetDirected}) as
#' BipartiteFourCyclesNodePowerA and BipartiteFourCyclesNodePowerB.
#'
#' These statistics sum the number of four-cycles at each node in one
#' bipartition (the first bipartition for \code{b1np4c} and the second
#' bipartition for \code{b2np4c}) raised to the power \eqn{\alpha} (\eqn{0 <
#'   \alpha \leq 1}). That is, the count at each node is raised to the power
#' \eqn{\alpha}, and then they are summed. In the terminology of
#' Wilson et al. (2017), this is an "\eqn{\alpha}-inside" weighting.
#'
#' The argument \code{fixed} indicates whether the parameter
#' \code{alpha} is to be fit as a curved exponential-family model (see
#' Hunter and Handcock, 2006). Currently, the fixed argument must be
#' set to TRUE, which means the weight parameter \code{alpha} is not
#' separately estimated and thus the model is not a CEF model.
#'
#' This term can only be used with undirected bipartite networks.
#'
#' Note that although there is a \code{fixed} parameter, the terms are
#' not yet able to handle a non-fixed \code{alpha} term so it must be
#' set to the default value \code{fixed=TRUE}.
#'
#' @references
#'
#' Hunter, D. R. and M. S. Handcock (2006). Inference in curved
#' exponential family models for networks. Journal of Computational
#' and Graphical Statistics, 15: 565-583.
#'
#' Stivala, A., Wang., P., and Lomi, A. (2025). Improving
#' exponential-family random graph models for bipartite
#' networks. arXiv preprint
#' 2502.01892. \url{https://arxiv.org/abs/2502.01892}
#'
#' Wang, P., Sharpe, K., Robins, G. L., and Pattison,
#' P. E. (2009). Exponential random graph (p*) models for affiliation
#' networks. Social Networks. 31(1): 12-25.
#'
#' Wilson, J. D., Denny, M. J., Bhamidi, S., Cranmer, S. J., &
#' Desmarais, B. A. (2017). Stochastic weighted graphs: Flexible model
#' specification and simulation. Social Networks, 49, 37-47.
#'
#' @author
#' Alex Stivala \email{alex.d.stivala@gmail.com}
#'
#' @examples
#'
#' library(ergm.terms.contrib)
#'
#' ## Construct the Four-fan-3 graph
#' fourfan.3.df <- data.frame(A = c(1, 1, 1, 1, 1, 1, 3, 3, 6, 6, 9, 9),
#'                            B = c(2, 4, 5, 7, 8, 10,2, 4, 5, 7, 8, 10))
#' fourfan.3.net <- network(fourfan.3.df, bipartite=TRUE, directed=FALSE)
#'
#' ## From Stivala et al. (2025):
#' ##
#' ##   ... in this graph, the nodes in mode B contribute more to the
#' ##   total as each one (of the six) is involved in exactly one
#' ##   four-cycle (and hence raising to the power of $\alpha$ still
#' ##   contributes one to the sum), while of the four nodes in mode A,
#' ##   three are involved in only one four-cycle, while the fourth is
#' ##   involved in three four-cycles and hence contributes only
#' ##   $3^\alpha \approx 1.73205$ (when $\alpha=0.5$).
#' ##
#' summary(fourfan.3.net ~ b1np4c(0.5, TRUE) + b2np4c(0.5, TRUE))
#'
#' \dontrun{
#' ## generate some graphs similar to 'zero.pos' in Fig. 5 and Fig. 6
#' ## of Stivala et al. (2025), with the parameter for b1np4c zero
#' ## and the parameter for b2np4c positive
#' ## (warning: takes about 2 minutes to run)
#' system.time( g.zero.pos <- simulate(network(75, bipartite=50,
#'                                             directed=FALSE) ~
#'               edges + gwb1degree(.5,TRUE) + gwb2degree(.5, TRUE) + b2np4c(1/5),
#'               coef=c(-5.0, 1.0, -0.5, 6.5),
#'               nsim=100,
#'               control = control.simulate(MCMC.burnin   = 100000,
#'                                          MCMC.interval = 100000)) )
#'
#' ## plot graph statistics to check sufficient burnin and interval
#' par(mfrow=c(2,3))
#' plot(summary(g.zero.pos ~ edges))
#' plot(summary(g.zero.pos ~ gwb1degree(.5,TRUE)))
#' plot(summary(g.zero.pos ~ gwb2degree(.5,TRUE)))
#' plot(summary(g.zero.pos ~ b2np4c(1/5)))
#' plot(summary(g.zero.pos ~ b1np4c(1/5)))
#' plot(summary(g.zero.pos ~ cycle(4)))
#'
#' ## plot visualization of one graph similar to
#' ## Fig. 6 of Stivala et al. (2025)
#' i <- 88
#' N.network <- network.size(g.zero.pos[[i]])
#' N.b1 <- get.network.attribute(g.zero.pos[[i]], 'bipartite')
#' N.b2 <- N.network - N.b1
#' par(mfrow=c(1,1))
#' plot(g.zero.pos[[i]],
#'      vertex.col= c(rep('red',N.b1), rep('blue',N.b2)),
#'      vertex.sides = c(rep(50, N.b1), rep(4, N.b2)))
#' }
#'
#' @concept undirected
#' @concept bipartite

InitErgmTerm.b1np4c <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed = FALSE, bipartite = TRUE, varnames = c("alpha", "fixed"), vartypes = c("numeric", "logical"), defaultvalues = list(0.5, TRUE), required = c(FALSE, FALSE))
  if(length(a$alpha) > 1)
    stop("The argument alpha to b1np4c expected a vector of length ",
         "1, but received a vector of length ",length(a$alpha))
  alpha = a$alpha[1]
  if (alpha <= 0.0 || alpha > 1)
    stop("The argument alpha to b1np4c must be 0 < alpha <= 1 ",
         "but received the value ", alpha)
  fixed <- a$fixed
  if (!fixed) {
    stop("The b1np4c term is not yet able to handle a ", 
         "non-fixed decay term.", call. = FALSE)
  }
  list(name = "b1np4c",
       coef.names = paste("b1np4c.fixed.", alpha, sep=""),
       inputs = c(alpha),
       dependence = TRUE,
       pkgname = "ergm.terms.contrib",
       auxiliaries = ~.spcache.net("UTP"))
}


#' @templateVar name b1np4c
#' @template ergmTerm-rdname
#' @aliases b2np4c-ergmTerm
#' @usage
#' #binary: b2np4c(alpha=0.5, fixed=TRUE)
#'

InitErgmTerm.b2np4c <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed = FALSE, bipartite = TRUE, varnames = c("alpha", "fixed"), vartypes = c("numeric", "logical"), defaultvalues = list(0.5, TRUE), required = c(FALSE, FALSE))
  if(length(a$alpha) > 1)
    stop("The argument alpha to b2np4c expected a vector of length ",
         "1, but received a vector of length ",length(a$alpha))
  alpha = a$alpha[1]
  if (alpha <= 0.0 || alpha > 1)
    stop("The argument alpha to b1np4c must be 0 < alpha <= 1 ",
         "but received the value ", alpha)
  fixed <- a$fixed
  if (!fixed) {
    stop("The b2np4c term is not yet able to handle a ", 
         "non-fixed decay term.", call. = FALSE)
  }
  list(name = "b2np4c",
       coef.names = paste("b2np4c.fixed.", alpha, sep=""),
       inputs = c(alpha),
       dependence = TRUE,
       pkgname = "ergm.terms.contrib",
       auxiliaries = ~.spcache.net("UTP"))
}


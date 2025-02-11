#  This software is distributed under the GPL-3 license.  


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
       pkgname = "ergm.terms.contrib")
}


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
       pkgname = "ergm.terms.contrib")
}


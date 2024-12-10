##
## File:    bnp4c.test.R
## Author:  Alex Stivala
## Created: December 2024
##
## Tests for the contributed ergm terms b1np4c and b2np4c, change
## statistics for number of 4-cycles at each node raised to a power
## alpha (0 < alpha <= 1) before summing (i.e. the "alpha-inside"
## weighting according to the description of Wilson et al. (2017)).
##
## b1np4c is for nodes in the first bipartition, and b2np4c is for nodes
## in the second bipartition.
##
## These statistics are described in Stivala et al. (2024). The b1np4c
## and b2np4c statistics were first implemented in EstimNetDirected
## (https://github.com/stivalaa/EstimNetDirected} as
## BipartiteFourCyclesNodePowerA and BipartiteFourCyclesNodePowerB.
##
## These statistics are only for undirected bipartite networks.
##
## For more details on the statistics and test data see
## Stivla et al. (2024) and
## https://github.com/stivalaa/EstimNetDirected/tree/master/Test/TestChangeStatsBipartite
##
##
## References:
##
## Stivala, A., Wang., P., and Lomi, A. (2024). Improving
## exponential-family random graph models for bipartite
## networks. Unpublished manusript.
##
## Wilson, J. D., Denny, M. J., Bhamidi, S., Cranmer, S. J., &
## Desmarais, B. A. (2017). Stochastic weighted graphs: Flexible model
## specification and simulation. Social Networks, 49, 37-47.
##

library(testthat)
library(network)
library(ergm)
library(ergm.terms.contrib)

##
## Read example networks. Note use of multiline strings as here document.
## These networks were originally constructed using igraph in
## by makeExampleNetworks.R in the Tests/TestChangeStatsBipartite directory
## of EstimNetDirected (see GitHub location in file header comment),
## so rather than adding a dependency on igraph here, easier to just
## read them in from the Pajek .net files created by makeExamplesNetworks.R
## (and unsure if I should add data files to the tests directory, easier
## to embed them as here documents in this script).
##

fourcycle.pajek.text <- 
'*Vertices 4 2
1 "1"
2 "3"
3 "2"
4 "4"
*Edges
1 3
3 2
2 4
1 4'
fourcycle.net <- read.paj(textConnection(fourcycle.pajek.text))

fourcycles3.pajek.text <-
'*Vertices 5 2
1 "1"
2 "3"
3 "2"
4 "4"
5 "5"
*Edges
1 3
3 2
2 4
1 4
1 5
2 5'
fourcycles3.net <- read.paj(textConnection(fourcycles3.pajek.text))

fourcycles3.revmode.pajek.text <-
'*Vertices 5 3
1 "2"
2 "4"
3 "5"
4 "1"
5 "3"
*Edges
4 1
1 5
5 2
4 2
4 3
5 3'
fourcycles3.revmode.net <- read.paj(textConnection(fourcycles3.revmode.pajek.text))

fourcycles6.pajek.text <-
'*Vertices 8 2
1 "1"
2 "3"
3 "2"
4 "4"
5 "5"
6 "6"
7 "7"
8 "8"
*Edges
1 3
3 2
2 4
1 4
1 5
2 5
1 6
2 6
1 7
2 7
1 8
2 8'
fourcycles6.net <- read.paj(textConnection(fourcycles6.pajek.text))


##
## Tests
##

test_that('bnp4c terms input validation', {
  ## alpha must be in (0, 1]
  expect_error(summary(fourcycle.net ~ b1np4c(-0.1)))
  expect_error(summary(fourcycle.net ~ b2np4c(-0.1)))
  expect_error(summary(fourcycle.net ~ b1np4c(0.0)))
  expect_error(summary(fourcycle.net ~ b2np4c(0.0)))
  expect_error(summary(fourcycle.net ~ b1np4c(1.01)))
  expect_error(summary(fourcycle.net ~ b2np4c(1.01)))
  expect_no_error(summary(fourcycle.net ~ b1np4c(1.0)))
  expect_no_error(summary(fourcycle.net ~ b2np4c(1.0)))

  ## can use a default value of alpha (0.5)
  expect_no_error(summary(fourcycle.net ~ b1np4c))
  expect_no_error(summary(fourcycle.net ~ b2np4c))
  

  ## fixed=FALSE is not (yet) supported
  expect_error(summary(fourcycle.net ~ b1np4c(0.4, FALSE)))
  expect_error(summary(fourcycle.net ~ b1np4c(fixed=FALSE)))
  expect_error(summary(fourcycle.net ~ b1np4c(alpha = 0.1,fixed=FALSE)))
  expect_no_error(summary(fourcycle.net ~ b1np4c(0.4, TRUE)))
  
})



test_that('bnp4c terms', {
  ##
  ## Tests on the example/test networks.
  ## See Table 4 of Stivala et al. (2024) and Tests/TestChangeStatsBipartite/
  ## in EstimNetDirected GitHub repository.
  ##

  ## Four-cycle
  expect_equal(as.vector(summary(fourcycle.net ~  b1np4c(0.5, TRUE))), 2)
  expect_equal(as.vector(summary(fourcycle.net ~  b2np4c(0.5, TRUE))), 2)

  ## Four-cycles-3
  expect_equal(as.vector(summary(fourcycles3.net ~ b1np4c(0.5, TRUE))), 3.4641, tolerance=1e-04)
  expect_equal(as.vector(summary(fourcycles3.net ~ b2np4c(0.5, TRUE))), 4.24264, tolerance=1e-05)

  ## Four-cycles-3-B
  expect_equal(as.vector(summary(fourcycles3.revmode.net ~ b1np4c(0.5, TRUE))), 4.24264, tolerance=1e-05)
  expect_equal(as.vector(summary(fourcycles3.revmode.net ~ b2np4c(0.5, TRUE))), 3.4641, tolerance=1e-04)


  ## Four-cycles-6
  expect_equal(as.vector(summary(fourcycles6.net ~ b1np4c(0.5, TRUE))), 7.74597, tolerance=1e-05)
  expect_equal(as.vector(summary(fourcycles6.net ~ b2np4c(0.5, TRUE))), 13.4164, tolerance=1e-04)

  ## Ten-cycle

  ## Eight-path

  ## Nine-star

  ## Nine-star-B

  ## Four-fan-3
  
})



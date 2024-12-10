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



##
## Tests on the example/test networks
##

test_that('bnp4c terms', {

  expect_equal(as.vector(summary(fourcycle.net ~  b1np4c(0.5, TRUE))), 2)
  expect_equal(as.vector(summary(fourcycle.net ~  b2np4c(0.5, TRUE))), 2)
  
})



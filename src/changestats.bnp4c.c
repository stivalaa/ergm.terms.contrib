/*  
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 */
/*****************************************************************************
 *
 * File:    changestats.bnp4c.c
 * Author:  Alex Stivala
 * Created: December 2024
 *
 * Change statistics for number of 4-cycles at each node raised to a
 * power alpha (0 < alpha <= 1) before summing (i.e. the "alpha-inside"
 * weighting according to the description of Wilson et al. (2017)).
 *
 * b1np4c is for nodes in the first bipartition, and b2np4c is for nodes
 * in the second bipartition.
 *
 * These statistics are described in Stivala et al. (2025). The b1np4c
 * and b2np4c statistics were first implemented in EstimNetDirected
 * (https://github.com/stivalaa/EstimNetDirected) as
 * BipartiteFourCyclesNodePowerA and BipartiteFourCyclesNodePowerB.
 *
 * These statistics are only for undirected bipartite networks.
 *
 *
 * References:
 *
 * Stivala, A., Wang., P., and Lomi, A. (2025). Improving
 * exponential-family random graph models for bipartite
 * networks. arXiv preprint 2502.01892. https://arxiv.org/abs/2502.01892
 *
 * Wilson, J. D., Denny, M. J., Bhamidi, S., Cranmer, S. J., &
 * Desmarais, B. A. (2017). Stochastic weighted graphs: Flexible model
 * specification and simulation. Social Networks, 49, 37-47.
 *
 *****************************************************************************/

#include <assert.h>

#include "ergm_changestat.h"
#include "ergm_storage.h"

/*****************************************************************************
 *
 * local functions
 *
 *****************************************************************************/

/*
 * Binomial coefficient n choose 2
 * simple and efficient special case instead of using CHOOSE(n, 2)
 */
static unsigned long n_choose_2(int n)
{
  if (n < 2) {
    return 0;
  }
  return (unsigned long)n * (n - 1) / 2;
}

/*
 * number of undirected two-paths for (i, j): paths i -- v -- j for some v
 *
 */
static unsigned int twopaths(Network *nwp, Vertex i, Vertex j)  {
  /* Note Network *nwp parameter has to be called nwp for use of macros */

  Vertex vnode, wnode;
  Edge edge1, edge2;
  unsigned int count = 0;

  /* It would be better if we could precompute or cache two-path counts like
     EstimNetDirected or statnet gwesp etc. cache */

  /* In an undirected network, each edge is only stored as (tail, head) where
     tail < head, so to step through all edges of a node it is necessary
     to step through all outedges and also through all inedges */
  STEP_THROUGH_OUTEDGES(i, edge1, vnode) {     /* i -- v */
    if (vnode == i || vnode == j)
      continue;
    STEP_THROUGH_OUTEDGES(j, edge2, wnode) {
      if (wnode == vnode)     /* v -- j */
        count++;
    }
    STEP_THROUGH_INEDGES(j, edge2, wnode) {
      if (wnode == vnode)      /* v -- j */
        count++;
    }
  }
  STEP_THROUGH_INEDGES(i, edge1, vnode) {     /* i -- v */
    if (vnode == i || vnode == j)
      continue;
    STEP_THROUGH_OUTEDGES(j, edge2, wnode) {
      if (wnode == vnode)     /* v -- j */
        count++;
    }
    STEP_THROUGH_INEDGES(j, edge2, wnode) {
      if (wnode == vnode)      /* v -- j */
        count++;
    }
  }
  return count;
}


/*
 * number of four-cycles that node unode is involved in
 *
 * the visited parameter is an array of N_NODES int allocated by caller
 */

static unsigned long num_fourcycles_node(Network *nwp, Vertex unode,
                                         int *visited)  {
  /* Note Network *nwp parameter has to be called nwp for use of macros */

  Vertex vnode, wnode;
  Edge edge1, edge2;
  unsigned long fourcycle_count = 0;

  memset(visited, 0, N_NODES);

  /*
    Note it seems that the Vertex is an int from 1 .. N_NODES
     (R indexing from 1 rather than C indexing from 0)
     so we subtract one when indexing the visited array.
     This code depends on this being the case.
  */
  
  /* this involves iterating over all nodes that are distance 2 from unode */
  /* In an undirected network, each edge is only stored as (tail, head) where
     tail < head, so to step through all edges of a node it is necessary
     to step through all outedges and also through all inedges */
  STEP_THROUGH_OUTEDGES(unode, edge1, vnode) {
    STEP_THROUGH_OUTEDGES(vnode, edge2, wnode) {
      if (wnode != unode && !visited[wnode-1]) {
        visited[wnode-1] = 1;
        fourcycle_count += n_choose_2(twopaths(nwp, unode, wnode));
      }
    }
    STEP_THROUGH_INEDGES(vnode, edge2, wnode) {
      if (wnode != unode && !visited[wnode-1]) {
        visited[wnode-1] = 1;
        fourcycle_count += n_choose_2(twopaths(nwp, unode, wnode));
      }
    }
  }
  STEP_THROUGH_INEDGES(unode, edge1, vnode) {
    STEP_THROUGH_OUTEDGES(vnode, edge2, wnode) {
      if (wnode != unode && !visited[wnode-1]) {
        visited[wnode-1] = 1;
        fourcycle_count += n_choose_2(twopaths(nwp, unode, wnode));
      }
    }
    STEP_THROUGH_INEDGES(vnode, edge2, wnode) {
      if (wnode != unode && !visited[wnode-1]) {
        visited[wnode-1] = 1;
        fourcycle_count += n_choose_2(twopaths(nwp, unode, wnode));
      }
    }
  }
  return fourcycle_count;
}


/*
 * Change statistic for number of four-cycles in undirected network
 * when edge i -- j is added
 */
static unsigned long change_fourcycles(Network *nwp, Vertex i, Vertex j) {
  /* Note Network *nwp parameter has to be called nwp for use of macros */
  Vertex vnode;
  Edge edge;
  unsigned long delta = 0;

  /* In an undirected network, each edge is only stored as (tail, head) where
     tail < head, so to step through all edges of a node it is necessary
     to step through all outedges and also through all inedges */
  STEP_THROUGH_OUTEDGES(i, edge, vnode) {
    delta += twopaths(nwp, vnode, j);
  }
  STEP_THROUGH_INEDGES(i, edge, vnode) {
    delta += twopaths(nwp, vnode, j);
  }
  return delta;
}




/*****************************************************************************
 *
 * change statistics functions
 *
 *****************************************************************************/


/*
 * From changestats.c:
 *
 *   For all bipartite networks:
 * It is assumed that in this bipartite network, the only edges are
 * of the form (b1, b2), where b1 is always strictly less
 * than b2.  In other words, the degree of a b1 is equivalent
 *  to its outdegree and the degree of a b2 is equivalent to its
 * indegree.
 */



/* Initializer: allocate private storage. */
I_CHANGESTAT_FN(i_b1np4c) {
  ALLOC_STORAGE(N_NODES, int, visited1); /* array of visited node flags */
}
I_CHANGESTAT_FN(i_b2np4c) {
  ALLOC_STORAGE(N_NODES, int, visited2); /* array of visited node flags */
}


/*
 * Change statistic function for b1np4c
 */
C_CHANGESTAT_FN(c_b1np4c) {
  double change, alpha, delta;
  unsigned long count, vcount;
  Vertex b1, b2;
  int is_delete;

  GET_STORAGE(int, visited1); /* Obtain a pointer to private storage
                                       and cast it to the correct type. */

  alpha = INPUT_PARAM[0];

  b1 = tail;
  b2 = head;
  is_delete = IS_UNDIRECTED_EDGE(b1, b2);
  /* NOTE: For a delete move, we actually toggle the edge ourselves here so
   * that the proposed edge is always NOT present for all the calculations
   * as they involve counting two-paths and four-cycles, on the assumption
   * that the proposed edge does not (yet) exist. We must therefore
   * also add it back afterwards to fit in with the standard logic.
   */
  if (is_delete) {
    TOGGLE(b1, b2);
  }
  if (IS_UNDIRECTED_EDGE(b1, b2)) error("Edge must not exist\n");

  /* Number of four-cycles the node is already involved in */
  count = num_fourcycles_node(nwp, b1, visited1);

  /* change statistic for four-cycles */
  delta = change_fourcycles(nwp, b1, b2);
  change = pow(count + delta, alpha) - pow(count, alpha);

  /* add contribution from sum over neighbours of b2 */
  EXEC_THROUGH_EDGES(b2, edge, vnode,  { /* step through edges of b2 */
    vcount = num_fourcycles_node(nwp, vnode, visited1);
    delta = twopaths(nwp, vnode, b1);
    change += pow(vcount + delta, alpha) - pow(vcount, alpha);
  });
  CHANGE_STAT[0] += IS_UNDIRECTED_EDGE(b1, b2) ? -change : change;
  /* For a delete move, we deleted the edge at the start, now add it again */
  if (is_delete) {
    TOGGLE(b1, b2);
    if (!IS_UNDIRECTED_EDGE(b1, b2)) error("Edge must exist\n");
  }
}



/*
 * Change statistic function for b2np4c
 */

C_CHANGESTAT_FN(c_b2np4c) {
  double change, alpha, delta;
  unsigned long count, vcount;
  Vertex b1, b2;
  int is_delete;

  GET_STORAGE(int, visited2); /* Obtain a pointer to private storage
                                       and cast it to the correct type. */

  alpha = INPUT_PARAM[0];

  b1 = tail;
  b2 = head;
  is_delete = IS_UNDIRECTED_EDGE(b1, b2);
  /* NOTE: For a delete move, we actually toggle the edge ourselves here so
   * that the proposed edge is always NOT present for all the calculations
   * as they involve counting two-paths and four-cycles, on the assumption
   * that the proposed edge does not (yet) exist. We must therefore
   * also add it back afterwards to fit in with the standard logic.
   */
  if (is_delete) {
    TOGGLE(b1, b2);
  }
  if (IS_UNDIRECTED_EDGE(b1, b2)) error("Edge must not exist\n");

  /* Number of four-cycles the node is already involved in */
  count = num_fourcycles_node(nwp, b2, visited2);

  /* change statistic for four-cycles */
  delta = change_fourcycles(nwp, b1, b2);
  change = pow(count + delta, alpha) - pow(count, alpha);

  /* add contribution from sum over neighbours of b1 */
  EXEC_THROUGH_EDGES(b1, edge, vnode, { /* step through edges of b1 */
    vcount = num_fourcycles_node(nwp, vnode, visited2);
    delta = twopaths(nwp, vnode, b2);
    change += pow(vcount + delta, alpha) - pow(vcount, alpha);
  })
  CHANGE_STAT[0] += is_delete ? -change : change;
  /* For a delete move, we deleted the edge at the start, now add it again */
  if (is_delete) {
    TOGGLE(b1, b2);
    if (!IS_UNDIRECTED_EDGE(b1, b2)) error("Edge must exist\n");
  }
}


/* It seems there is no need to free private storage with a Finalizer
   (f_) function, it is done by statnet elsewhwere (see point 9
   ModelDestroy() is called in the API definition:
   https://cran.r-project.org/web/packages/ergm/vignettes/Terms-API.html
*/

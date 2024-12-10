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
 * weighting according to the description of Wilson et al. (2017).
 *
 * b1np4c is for nodes in the first bipartition, and b2np4c is for nodes
 * in the second bipartition.
 *
 * These statistics are described in Stivala et al. (2024).  The
 * b1np4c and b2np4c statistics were first implemented in
 * EstimNetDirected (https://github.com/stivalaa/EstimNetDirected}
 * as BipartiteFourCyclesNodePowerA and BipartiteFourCyclesNodePowerB.
 *
 * These statistics are only for undirected bipartite networks.
 *
 *
 * References:
 *
 * Stivala, A., Wang., P., and Lomi, A. (2024). Improving
 * exponential-family random graph models for bipartite
 * networks. Unpublished manusript.
 *
 * Wilson, J. D., Denny, M. J., Bhamidi, S., Cranmer, S. J., &
 * Desmarais, B. A. (2017). Stochastic weighted graphs: Flexible model
 * specification and simulation. Social Networks, 49, 37-47.
 *
 *****************************************************************************/

#include "changestats.bnp4c.h"


/*****************************************************************************
 *
 * local functions
 *
 *****************************************************************************/


/*
 * number of four-cycles that node unode is involved in
 *
 */
static unsigned long num_fourcycles_node(Network *nwp, Vertex unode)  {
  /* Note Network *nwp parameter has to be called nwp for use of macros */

  Vertex vnode, wnode;
  Edge edge1, edge2;
  int *visited = R_Calloc(N_NODES, int); /* list of visited nodes */
  unsigned long fourcycle_count = 0;

  /* this involves iterating over all nodes that are distance 2 from unode */
  /* In an undirected network, each edge is only stored as (tail, head) where
     tail < head, so to step through all edges of a node it is necessary
     to step through all outedges and also through all inedges */
  STEP_THROUGH_OUTEDGES(unode, edge1, vnode) {
    STEP_THROUGH_OUTEDGES(vnode, edge2, wnode) {
      if (wnode != unode && !visited[wnode]) {
        visited[wnode] = 1;
        fourcycle_count += 1; /*FIXME: CHOOSE(L2(unode,wnode), 2) */
      }
    }
    STEP_THROUGH_INEDGES(vnode, edge2, wnode) {
      if (wnode != unode && !visited[wnode]) {
        visited[wnode] = 1;
        fourcycle_count += 1; /*FIXME: CHOOSE(L2(unode,wnode), 2) */
      }
    }
  }
  STEP_THROUGH_INEDGES(unode, edge1, vnode) {
    STEP_THROUGH_OUTEDGES(vnode, edge2, wnode) {
      if (wnode != unode && !visited[wnode]) {
        visited[wnode] = 1;
        fourcycle_count += 1; /*FIXME: CHOOSE(L2(unode,wnode), 2) */
      }
    }
    STEP_THROUGH_INEDGES(vnode, edge2, wnode) {
      if (wnode != unode && !visited[wnode]) {
        visited[wnode] = 1;
        fourcycle_count += 1; /*FIXME: CHOOSE(L2(unode,wnode), 2) */
      }
    }
  }
  R_Free(visited);
  return fourcycle_count;
}



/*****************************************************************************
 *
 * change statistics functions
 *
 *****************************************************************************/


D_CHANGESTAT_FN(d_b1np4c) {
  double change, alpha;
  Vertex tail, head, unode, vnode, wnode;
  Edge edge1, edge2;
  int i;
  unsigned long fourcycle_count;
  

  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    tail = TAIL(i);
    head = HEAD(i);
    alpha = INPUT_PARAM[0];

    fourcycle_count = num_fourcycles_node(nwp, tail);
    change = fourcycle_count;

    CHANGE_STAT[0] += IS_OUTEDGE(tail, head) ? -change : change;
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}



D_CHANGESTAT_FN(d_b2np4c) { 
  int i;
  
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

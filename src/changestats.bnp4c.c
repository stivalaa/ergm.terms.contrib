/*  
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 */
#include "changestats.gwbnsp.h"

/*
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
 * Rereferences:
 *
 * Stivala, A., Wang., P., and Lomi, A. (2024). Improving
 * exponential-family random graph models for bipartite
 * networks. Unpublished manusript.
 *
 * Wilson, J. D., Denny, M. J., Bhamidi, S., Cranmer, S. J., &
 * Desmarais, B. A. (2017). Stochastic weighted graphs: Flexible model
 * specification and simulation. Social Networks, 49, 37-47.
 */

D_CHANGESTAT_FN(d_b1np4c) {
  int i;
  
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
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

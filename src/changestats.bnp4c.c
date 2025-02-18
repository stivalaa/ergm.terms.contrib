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

#include "ergm_changestat.h"
#include "ergm_storage.h"



#ifdef DEBUG
#define DEBUG_PRINT(x) printf("DEBUG BNP4C: "); printf x
#else
#define DEBUG_PRINT(x) /* nothing */
#endif

/*****************************************************************************
 *
 * type definitions
 *
 *****************************************************************************/


 /* private storage struct for b1np4c and b2np4c */
/* These are just allocated as arrays of length N_NODES for simplicity;
 * they could perhaps more appropriately be hash tables, but, unlike the case
 * for dyadic properties and so would require N^2 storage, it is not
 * too inefficient to require N storage here.
 * Note than indexiing in these two arrays is 0..N_NODES-1 so
 * must subrtact one when indexigin from a Vertex which
 * is 1..N_NODES (R style) */
typedef struct bnp4c_storage_s {
  int            *visited;         /* visited flag for each node
                                      (working storage) */
  unsigned long  *fourcycle_count; /* number of four-cycles at each node */
} bnp4c_storage_t;


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
 * If ignore1 and ignore2 are nonzero then the edge ignore1 -- ignore2
 * is not not included in traversals (treated as if it does not exist)
 *
 * (valid Vertices are numbered 1..N_NODES)
 */
static unsigned int twopaths(Network *nwp, Vertex i, Vertex j,
                             Vertex ignore1, Vertex ignore2)  {
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
    if (vnode == i || vnode == j ||
        (i == ignore1 && vnode == ignore2) || i == (ignore2 && vnode == ignore1))
      continue;
    STEP_THROUGH_OUTEDGES(j, edge2, wnode) {
      if ((j == ignore1 && wnode == ignore2) || (j == ignore2 && wnode == ignore1))
        continue;
      if (wnode == vnode)     /* v -- j */
        count++;
    }
    STEP_THROUGH_INEDGES(j, edge2, wnode) {
      if ((j == ignore1 && wnode == ignore2) || (j == ignore2 && wnode == ignore1))
        continue;
      if (wnode == vnode)      /* v -- j */
        count++;
    }
  }
  STEP_THROUGH_INEDGES(i, edge1, vnode) {     /* i -- v */
    if (vnode == i || vnode == j ||
        (i == ignore1 && vnode == ignore2) || (i == ignore2 && vnode == ignore1))
      continue;
    STEP_THROUGH_OUTEDGES(j, edge2, wnode) {
      if ((j == ignore1 && wnode == ignore2) || (j == ignore2 && wnode == ignore1))
        continue;
      if (wnode == vnode)     /* v -- j */
        count++;
    }
    STEP_THROUGH_INEDGES(j, edge2, wnode) {
      if ((j == ignore1 && wnode == ignore2) || (j == ignore2 && wnode == ignore1))
        continue;
      if (wnode == vnode)      /* v -- j */
        count++;
    }
  }
  return count;
}


/*
 * Number of four-cycles that node unode is involved in/
 *
 * The bnp4c_storage_t sto parameter containing the visited working
 * storage is allocated by the caller.
 *
 */
static unsigned long num_fourcycles_node(Network *nwp, Vertex unode,
                                         bnp4c_storage_t *sto)  {
  /* Note Network *nwp parameter has to be called nwp for use of macros */

  Vertex vnode, wnode;
  Edge edge1, edge2;
  unsigned long fourcycle_count = 0;
  int *visited = sto->visited;
  
  memset(visited, 0, N_NODES*sizeof(visited[0]));
  
  /* this involves iterating over all nodes that are distance 2 from unode */
  /* In an undirected network, each edge is only stored as (tail, head) where
     tail < head, so to step through all edges of a node it is necessary
     to step through all outedges and also through all inedges */
  STEP_THROUGH_OUTEDGES(unode, edge1, vnode) {
    STEP_THROUGH_OUTEDGES(vnode, edge2, wnode) {
      if (wnode != unode && !visited[wnode-1]) {
        visited[wnode-1] = 1;
        fourcycle_count += n_choose_2(twopaths(nwp, unode, wnode, 0, 0));
      }
    }
    STEP_THROUGH_INEDGES(vnode, edge2, wnode) {
      if (wnode != unode && !visited[wnode-1]) {
        visited[wnode-1] = 1;
        fourcycle_count += n_choose_2(twopaths(nwp, unode, wnode, 0, 0));
      }
    }
  }
  STEP_THROUGH_INEDGES(unode, edge1, vnode) {
    STEP_THROUGH_OUTEDGES(vnode, edge2, wnode) {
      if (wnode != unode && !visited[wnode-1]) {
        visited[wnode-1] = 1;
        fourcycle_count += n_choose_2(twopaths(nwp, unode, wnode, 0, 0));
      }
    }
    STEP_THROUGH_INEDGES(vnode, edge2, wnode) {
      if (wnode != unode && !visited[wnode-1]) {
        visited[wnode-1] = 1;
        fourcycle_count += n_choose_2(twopaths(nwp, unode, wnode, 0, 0));
      }
    }
  }
  return fourcycle_count;
}


/*
 * Change statistic for number of four-cycles in undirected network
 * when edge i -- j is added
 *
 * If ignore1 and ignore2 are nonzero then the edge ignore1 -- ignore2
 * is not not included in traversals (treated as if it does not exist)
 *
 * (valid Vertices are numbered 1..N_NODES)
 */
static unsigned long change_fourcycles(Network *nwp, Vertex i, Vertex j,
                                       Vertex ignore1, Vertex ignore2) {
  /* Note Network *nwp parameter has to be called nwp for use of macros */
  Vertex vnode;
  Edge edge;
  unsigned long delta = 0;

  /* In an undirected network, each edge is only stored as (tail, head) where
     tail < head, so to step through all edges of a node it is necessary
     to step through all outedges and also through all inedges */
  STEP_THROUGH_OUTEDGES(i, edge, vnode) {
    if ((i == ignore1 && vnode == ignore2) || (i == ignore2 && vnode == ignore1))
      continue;
    delta += twopaths(nwp, vnode, j, ignore1, ignore2);
  }
  STEP_THROUGH_INEDGES(i, edge, vnode) {
    if ((i == ignore1 && vnode == ignore2) || (i == ignore2 && vnode == ignore1))
      continue;
    delta += twopaths(nwp, vnode, j, ignore1, ignore2);
  }
  return delta;
}




/*****************************************************************************
 *
 * Exported: Initializer, updater, change statistics, and finalizer functions
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


/*
 * The change statitsics functions (c_) use private storage to communicate
 * with the update functions (u_). The c_ functions use the four-cycles
 * counts in fourcycle_count, which are updated by the u_ functions.
 */


/* Initializer: allocate private storage and store number of four-cycles
   at each node . */
I_CHANGESTAT_FN(i_b1np4c) {
  ALLOC_STORAGE(1, bnp4c_storage_t, sto1);
  sto1->visited = R_Calloc(N_NODES, int);
  sto1->fourcycle_count = R_Calloc(N_NODES, unsigned long);
  for (int i = 1; i <= N_NODES; i++) {
    sto1->fourcycle_count[i-1] = num_fourcycles_node(nwp, i, sto1);
     DEBUG_PRINT(("i_b1np4c %d set to %lu\n", i, sto1->fourcycle_count[i-1]));
  }
}

I_CHANGESTAT_FN(i_b2np4c) {
  ALLOC_STORAGE(1, bnp4c_storage_t, sto2);
  sto2->visited = R_Calloc(N_NODES, int);
  sto2->fourcycle_count = R_Calloc(N_NODES, unsigned long);
  for (int i = 1; i <= N_NODES; i++) {
    sto2->fourcycle_count[i-1] = num_fourcycles_node(nwp, i, sto2);
     DEBUG_PRINT(("i_b2np4c %d set to %lu\n", i, sto2->fourcycle_count[i-1]));
  }
}

/* Updater: will be called when toggling (tail, head) with state
 * edgestate is imminent.
 *
 * The updater function keeps the private storage up-to-date with
 * number of four-cycles at each node, which is used in the change
 * statistics function.
 *
 */
/* Note the u_ function is called AFTER the c_ function (and only
   if the move is accepted), see API definition
   https://cran.r-project.org/web/packages/ergm/vignettes/Terms-API.html
*/
U_CHANGESTAT_FN(u_b1np4c) {
  unsigned long delta;
  Vertex b1 = tail, b2 = head;
  int is_delete = edgestate;

  DEBUG_PRINT(("XXX u_b1np4c entered b1 = %d is_delete = %d\n", b1, is_delete));
  GET_STORAGE(bnp4c_storage_t, sto1); /* Obtain a pointer to private storage
                                         and cast it to the correct type. */
  /* change statistic for four-cycles */
  delta = change_fourcycles(nwp, b1, b2, b1, b2);
  sto1->fourcycle_count[b1-1] += is_delete ? -delta : delta;
  DEBUG_PRINT(("u_b1np4c for %d added %ld to get %lu\n", b1, delta,sto1->fourcycle_count[b1-1]));
  /* add also have to update neighbours of b2 */
  EXEC_THROUGH_EDGES(b2, edge, vnode,  { /* step through edges of b2 */
    if (vnode == b1) continue; /* except for b1 -- b2 edge (if is_delete) */   
    delta = twopaths(nwp, vnode, b1, b1, b2);
    sto1->fourcycle_count[vnode-1] += is_delete ? -delta : delta;
    DEBUG_PRINT(("u_b1np4c for %d added %ld to get %lu\n",vnode-1, delta,sto1->fourcycle_count[vnode-1]));
  });
  DEBUG_PRINT(("XXX u_b1np4c exit\n"));
}

U_CHANGESTAT_FN(u_b2np4c) {
  unsigned long delta;
  int is_delete = edgestate;
  Vertex b1 = tail, b2 = head;

  DEBUG_PRINT(("XXX u_b2np4c entered b2 = %d is_delete = %d\n", b2, is_delete));
  GET_STORAGE(bnp4c_storage_t, sto2); /* Obtain a pointer to private storage
                                         and cast it to the correct type. */
  /* change statistic for four-cycles */
  delta = change_fourcycles(nwp, b1, b2, b1, b2);
  sto2->fourcycle_count[b2-1] += is_delete ? -delta : delta;
  DEBUG_PRINT(("u_b2np4c [1] %d added %ld to get %lu\n", b2-1, is_delete ? -delta : delta, sto2->fourcycle_count[b2-1]));
  /* also changes four-cycle counts for neighbours of b1 */
  EXEC_THROUGH_EDGES(b1, edge, vnode, { /* step through edges of b1 */
    if (vnode == b2) continue; /* except for b1 -- b2 edge (if is_delete) */
    delta = twopaths(nwp, vnode, b2, b1, b2);
    sto2->fourcycle_count[vnode-1] += is_delete ? - delta : delta;
    DEBUG_PRINT(("u_b2np4c [2] %d added %ld to get %lu\n", vnode-1, is_delete ? -delta : delta, sto2->fourcycle_count[vnode-1]));    
  })
  DEBUG_PRINT(("XXX u_b2np4c exit\n"));  
}


/*
 * Change statistic function for b1np4c
 */
C_CHANGESTAT_FN(c_b1np4c) {
  double change, alpha, delta;
  unsigned long count, vcount;
  Vertex b1, b2;
  int is_delete;

  DEBUG_PRINT(("XXX c_b1np4c entered b1 = %d\n", tail));
  
  GET_STORAGE(bnp4c_storage_t, sto1); /* Obtain a pointer to private storage
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
  count = sto1->fourcycle_count[b1-1];
#ifdef DEBUG
  if (num_fourcycles_node(nwp, b1, sto1) != count) error("b1np4c incorrect fourcycle count [1] for %d correct %lu got %lu\n", b1, num_fourcycles_node(nwp, b1, sto1), count);
#endif /* DEBUG */
  /* change statistic for four-cycles */
  delta = change_fourcycles(nwp, b1, b2, 0, 0);
  change = pow(count + delta, alpha) - pow(count, alpha);

  /* add contribution from sum over neighbours of b2 */
  EXEC_THROUGH_EDGES(b2, edge, vnode,  { /* step through edges of b2 */
    vcount = sto1->fourcycle_count[vnode-1];
/* #ifdef DEBUG */
/*     /\* using #ifdef inside macro (EXEC_THROUGH_EDGES) gives compiler warning *\/ */
/*     if (num_fourcycles_node(nwp, vnode, sto1) != vcount) error("b1np4c incorrect fourcycle count [2] for %d correct %lu got %lu\n", vnode, num_fourcycles_node(nwp, vnode, sto1), vcount); */
/* #endif /\* DEBUG *\/ */
    delta = twopaths(nwp, vnode, b1, 0, 0);
    change += pow(vcount + delta, alpha) - pow(vcount, alpha);
  });
  CHANGE_STAT[0] += IS_UNDIRECTED_EDGE(b1, b2) ? -change : change;
  /* For a delete move, we deleted the edge at the start, now add it again */
  if (is_delete) {
    TOGGLE(b1, b2);
    if (!IS_UNDIRECTED_EDGE(b1, b2)) error("Edge must exist\n");
  }
  DEBUG_PRINT(("XXX c_b1np4c exit\n"));
}



/*
 * Change statistic function for b2np4c
 *
 */

C_CHANGESTAT_FN(c_b2np4c) {
  double change, alpha, delta;
  unsigned long count, vcount;
  Vertex b1, b2;
  int is_delete;

  DEBUG_PRINT(("XXX c_b2np4c entered b2 = %d edgestate = %d\n", head, edgestate));

  GET_STORAGE(bnp4c_storage_t, sto2); /* Obtain a pointer to private storage
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
    DEBUG_PRINT((" c_b2np4c is_delete TOGGLE\n"));
    TOGGLE(b1, b2);
  }
  if (IS_UNDIRECTED_EDGE(b1, b2)) error("Edge must not exist\n");

  /* Number of four-cycles the node is already involved in */
  count = sto2->fourcycle_count[b2-1];
#ifdef DEBUG
  if (num_fourcycles_node(nwp, b2, sto2) != count) error("b2np4c incorrect fourcycle count [1] for %d correct %lu got %lu\n", b2, num_fourcycles_node(nwp, b2, sto2), count);
#endif /* DEBUG */
  /* change statistic for four-cycles */
  delta = change_fourcycles(nwp, b1, b2, 0, 0);
  change = pow(count + delta, alpha) - pow(count, alpha);

  /* add contribution from sum over neighbours of b1 */
  EXEC_THROUGH_EDGES(b1, edge, vnode, { /* step through edges of b1 */
    vcount = sto2->fourcycle_count[vnode-1];
/* #ifdef DEBUG */
/*     /\* using #ifdef inside macro (EXEC_THROUGH_EDGES) gives compiler warning *\/ */
/*     if (num_fourcycles_node(nwp, vnode, sto2) != vcount) error("b2np4c incorrect fourcycle count [2] for %d correct %lu got %lu\n", vnode, num_fourcycles_node(nwp, vnode, sto2), vcount); */
/* #endif /\* DEBUG *\/ */
    delta = twopaths(nwp, vnode, b2, 0, 0);
    change += pow(vcount + delta, alpha) - pow(vcount, alpha);
  })
  CHANGE_STAT[0] += is_delete ? -change : change;
  /* For a delete move, we deleted the edge at the start, now add it again */
  if (is_delete) {
    TOGGLE(b1, b2);
    if (!IS_UNDIRECTED_EDGE(b1, b2)) error("Edge must exist\n");
  }
  DEBUG_PRINT(("XXX c_b2np4c exit\n"));
}


/* Finalizer: free private storage. */

/* There is no need to free private storage allocated by ALLOC_STORAGE
   with a Finalizer (f_) function, it is done by statnet elsewhwere,
   see point 9 ModelDestroy() is called in the API definition:
   https://cran.r-project.org/web/packages/ergm/vignettes/Terms-API.html

   However we do need to free the storage inside this that we allocated
   in the Initializer functions.
*/

F_CHANGESTAT_FN(f_b1np4c) {
  DEBUG_PRINT(("XXX Finalize b1np4c\n"));
  GET_STORAGE(bnp4c_storage_t, sto1);
  R_Free(sto1->visited);
  R_Free(sto1->fourcycle_count);
}

F_CHANGESTAT_FN(f_b2np4c) {
  DEBUG_PRINT(("XXX Finalize b2np4c\n"));
  GET_STORAGE(bnp4c_storage_t, sto2);
  R_Free(sto2->visited);
  R_Free(sto2->fourcycle_count);
}

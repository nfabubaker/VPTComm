/**
 * @author      : Nabil Abubaker (nabil.abubaker@bilkent.edu.tr)
 * @file        : initial
 * @created     : Perşembe Tem 16, 2020 17:23:02 +03
 */

#ifndef INITIAL_H

#define INITIAL_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <assert.h>
#include "vptcomm.h"
#include "stfw_def.h"
#include "vpt.h"
/*! \enum init_stratigies
 *
 *  Detailed description
 */

enum SparseOp {SPOP_EXPAND, SPOP_REDUCE};
enum init_stratigies {DRYRUN_HEAVY_P2P, DRYRUN_SPARSE_P2P  };

void init_cm_p(stfw_comm *cm);
void init_stfw(Comm *comm, vpt *gm,  int STFW_INIT_STRATEGY);
void init_stfw_ExpRed(Comm *comm, vpt *gm);
void init_stfw_nh(Comm *comm, vpt *gm);
void stfw(Comm *comm, vpt *gm);
void stfw_nh(Comm *comm, vpt *gm);
void stfw_ExpRed(Comm *comm, vpt *gm ,  void (*reduce_op)(vptcomm_real_t*, vptcomm_real_t*, int) , void (*reduce_op_init)(vptcomm_real_t *, int), enum SparseOp SpOp);
#endif /* end of include guard INITIAL_H */


/**
 * @author      : Nabil Abubaker (nabil.abubaker@bilkent.edu.tr)
 * @file        : vpt.h
 * @created     : Monday Sep 05, 2022 12:18:29 +03
 */

#ifndef GMESH_H

#define GMESH_H

#include "stfw_def.h"

/* communication order in the mesh */
/* should be parameter */
typedef enum _dim_order {
	D_LINEAR,
	D_LINEARREV,
	D_RANDOM
} dim_order;
/* n-dimensional generic mesh (process specific) */
typedef struct _vpt
{
	/* general info */
	int	 nd;					/* number of dimensions */
	int *order;					/* dimension ordering */
	int *size;					/* dimension size */
	int *csize;					/* cumulative size for convenience */
	int *nbit;					/* number of bits for each dim */
	int *cnbit;					/* cumulative number of bits for convenience */
	int *mask;
    int *map;                   /* vpt-aware mapping: map[i]=p maps MPI rank i to processor p in the vpt*/
    int *inv_map;              /* vpt-aware mapping: inv_map[p]=i maps processor p in the vpt to MPI rank i*/

	/* process specific */
	int			 *crd;			/* coordinates of the proc */
	stfw_comm	 *sc;			/* communicators for each dim */
	int 		 *comm_dim;		/* which proc will be communicated in which dim, size: nprocs */
} vpt;

void assign_crd(vpt *gm);
int get_nghbr(vpt *gm, int dim, /* on dimension d */ int x               /* xth order in d */);
void print_vpt(vpt *gm);
void get_all_crd(vpt *gm, int proc_rank, int *crd);
int get_crd(vpt *gm, int proc_rank, int d);
void init_comm_dim(vpt *gm);
void init_vpt(vpt *gm, dim_order ord, int *map);
void free_vpt(vpt *gm);
void gm_cstats_na(vpt *gm, int *maxSendVol, int *maxRecvVol, int *totalVol , int *maxSendMsgs, int *maxRecvMsgs, int *totalMsgs, int factor);
#endif /* end of include guard GMESH_H */

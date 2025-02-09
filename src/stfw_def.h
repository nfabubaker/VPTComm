/**
 * @author      : Nabil Abubaker (nabil.abubaker@bilkent.edu.tr)
 * @file        : initial
 * @created     : Perşembe Tem 16, 2020 17:23:02 +03
 */

#ifndef STFWDEF_H

#define STFWDEF_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <assert.h>
#include "vptcomm.h"
#include "util.h"
/* #define D_SPEC_STFW */

enum CommType {STFW_NAIVE, STFW_EXPRED};

#define DEFAULT_B_VALUE 1.0
#define DEFAULT_X_VALUE 1.0

#define	ROWDIM(proc,dimc) (proc/dimc)
#define COLDIM(proc,dimc) (proc%dimc)

#define TAG1 1 

#define MAX(x,y) ((x)>(y))?x:y
//#define sdbgfp stderr

extern int _DE; /* short for dual enabled, will be set by the user to enable dual communication during initialization */

typedef struct _comm{ /* for a processor P:  */
   enum CommType commType;
  int sno   ; /* number of processors that P communicates with */
  int rno   ;
  int* sprocs ; /* the processors that P communicates with */  
  int* rprocs ;
  int* ssizes ; /* the no of words to be sent */
  int* rsizes ; /* the no of words to be received */
  int **sendinds;
  int **recvinds;
  vptcomm_real_t** recv ; /* each (vptcomm_real_t*) in the array points to the corresponding
		     element of p vector */
  vptcomm_real_t** send ; /* each (vptcomm_real_t*) pointed in the array must be allocated in 
		     the beginning */
  int *indsmap;
  vptcomm_real_t **send_unit_ptr; /* pointer per unit */
  int* LGM; /* local to global map of indices in sendinds and recvinds */
  int maxunitID;
  int unitSize;

   /* for buffers, indexing starts from 0. For example, for send[i][j], i
    starts from 1, but j starts from 0, since send[i] is a buffer */
} Comm ;




/* per proc */
typedef struct _stfw_comm_stats_t
{
	int vsend, vrecv;
	int msend, mrecv;
	int hdr_send, hdr_recv;
	int vsend_with_hdr, vrecv_with_hdr;
} stfw_comm_stats_t;


#define NE -1




typedef struct _stfw_comm
{
	/* forward buffer */
	vptcomm_real_t	*fw_buf;
    Comm *p2pcomm; // pointer to the p2p comm 
	int			 fw_np;
	int			*fw_procs;
	int			*fw_locs;		/* size fw_np+1, cumulative */
    int         *fw_offset;  
    int         *fw_offset_d;   /* for dual comm */
	int 		*fw_imap;		/* of size d, points to indices of fw_procs */
    int         *fw_imap_d;     /* for dual communication, initialized according to st_buf */    


	/* store buffer */
	vptcomm_real_t	*st_buf;
	int			 st_np;
	int			*st_procs;
	int			*st_locs;		/* size st_np+1, cumulative */
	MPI_Request *st_reqs;		/* for nonblocking recvs */
	MPI_Status	*st_stts;		/* for nonblocking recvs */
	MPI_Request *st_reqs_d;		/* for dual communication */
	MPI_Status	*st_stts_d;		/* for dual comm */


	
	/* packet processing */
	int			 *pck_cnt;		/* cumulative, size st_np+1 */
	int			 *pck_sizes;	/* sizes of all packets */
	vptcomm_real_t	**pck_cp;		/* points to forwarding locations or input vector */

    /* unit processing for ExpRed */
    int unitSize;
    int     *unit_RepCnt; /* unit replication factor, used in expand */
    vptcomm_real_t **unit_ptr; /* pointer per unit */

	#ifdef D_SPEC_STFW
	/* debugging purposes */
	int *pck_fw_idx;
	int *pck_src;
	int	*pck_dest;
	#endif

	/* stfw_comm_stats_t *stfw_stats; */
} stfw_comm;





/* 
void init_gmesh(gmesh* gm, dim_order order); 

 */








#endif /* end of include guard INITIAL_H */


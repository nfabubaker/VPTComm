#ifndef DEF_H 

#define DEF_H
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <assert.h>
#include "vptcomm.h"
/* #define D_SPEC_STFW */


#define DEFAULT_B_VALUE 1.0
#define DEFAULT_X_VALUE 1.0

#define	ROWDIM(proc,dimc) (proc/dimc)
#define COLDIM(proc,dimc) (proc%dimc)

#define TAG1 1 

#define MAX(x,y) ((x)>(y))?x:y
//#define sdbgfp stderr

extern int _DE; /* short for dual enabled, will be set by the user to enable dual communication during initialization */

typedef struct{ /* for a processor P:  */
  int sno   ; /* number of processors that P communicates with */
  int rno   ;
  int* sprocs ; /* the processors that P communicates with */  
  int* rprocs ;
  int* ssizes ; /* the no of words to be sent */
  int* rsizes ; /* the no of words to be received */
  vptcomm_real_t** recv ; /* each (vptcomm_real_t*) in the array points to the corresponding
		     element of p vector */
  vptcomm_real_t** send ; /* each (vptcomm_real_t*) pointed in the array must be allocated in 
		     the beginning */
  int** send_map ; /* each send_map[i][j] holds the row number to be copied to send[i][j] */

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


typedef struct _stfw_comm_t
{
	vptcomm_real_t	*st_buf;
	int			 st_np;
	int			*st_procs;
	int			*st_locs;		/* size st_np+1, cumulative */
	
	vptcomm_real_t	*fw_buf;
	int			 fw_np;
	int			*fw_procs;
	int			*fw_locs;		/* size fw_np+1, cumulative */

	stfw_comm_stats_t *stfw_stats;
} stfw_comm_t;


#define NE -1
/* Create an n-dimensional mesh here in future */
/* now only for 2D - processor-specific information */
typedef struct _mesh2d_t
{
	int nrd, ncd;				/* global */
	int	rd, cd;					/* mine */

	/* communicators */
	stfw_comm_t *row_cm;
	stfw_comm_t *col_cm;

	/* inverse comm map - index is of comm->no */
	int *inv_row_map;
	int *inv_col_map;

	/* inverse comm map - dimension to procs in row/col dimension */
	int *row_fw_imap;
	int *col_fw_imap;

	/* inverse recv_map, index: proc, value idx in comm */
	int *inv_recv_map;
} mesh2d_t;



/* ========================================================================== */
/* ========================================================================== */
/* ========================================================================== */
/* ========================================================================== */
/* ========================================================================== */

/* communication order in the mesh */
/* should be parameter */
typedef enum _dim_order {
	D_LINEAR,
	D_RANDOM
} dim_order;



typedef struct _stfw_comm
{
	/* forward buffer */
	vptcomm_real_t	*fw_buf;
	int			 fw_np;
	int			*fw_procs;
	int			*fw_locs;		/* size fw_np+1, cumulative */
	int 		*fw_offset;		/* the most current locs, size: fw_np. equal to
								 * fw_locs at the beginning but updated when
								 * running the alg */
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
	vptcomm_real_t	**pck_cp;		/* points to forwarding locations or input
								 * vector */

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


vptcomm_real_t Innerp(vptcomm_real_t* v1, vptcomm_real_t* v2, int size) ;


int no_procs, no_right, no_left, no_up, no_down, halfX, halfY ;
int my_id, my_x, my_y, my_left, my_right, my_up, my_down ;



char inp_hg_name[100]  ;
char inp_b_name[100] ;
char inp_parts_name[100] ;

int  dimX ;
int  dimY ;



/* OGUZ-EDIT */
/* NABIL commented: */
 FILE *sdbgfp;
char dbg_fname[100];

 /**
 * @brief timer
 */
typedef struct _tmr_t
{
	struct timespec ts_beg;
	struct timespec ts_end;
	double			elapsed;	/* in nanosec */
} tmr_t;

void start_timer (tmr_t *t);
void stop_timer (tmr_t *t);
void init_timer (tmr_t *t);
#endif /* end of include guard DEF_H*/

/**
 * @author      : Nabil Abubaker (nabil.abubaker@bilkent.edu.tr)
 * @file        : libstfw
 * @created     : Perşembe Tem 16, 2020 14:07:00 +03
 * @TODO:- add vpt, comm, stfw_cm as global variables.
 *       - invistigate dual communication, do we need a dual stfw_cm
 * (re-initialize  ?
 */
#include "mpi.h"
#include "stfw.h"
#include "vptcomm.h"
#include "stfw_def.h"
#include "vpt.h"
#include <stdlib.h>
#include <string.h>

#ifdef D_SPEC_STFW
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

/* Thanks to USERNAME=(adamk) from stackoverflow */
#define swap(x, y)                                                             \
    do {                                                                         \
        unsigned char swap_temp[sizeof(x) == sizeof(y) ? (signed)sizeof(x) : -1];  \
        memcpy(swap_temp, &y, sizeof(x));                                          \
        memcpy(&y, &x, sizeof(x));                                                 \
        memcpy(&x, swap_temp, sizeof(x));                                          \
    } while (0)

int _ninst = 0;
vpt *gmG;
Comm *commG;
int _DE = 0;
/* #####   FUNCTION DEFINITIONS  -  EXPORTED FUNCTIONS
 * ############################ */

void _setdoublezero(double *arr, int size) {
    int i;
    for (i = 0; i < size; ++i) {
        arr[i] = 0.0;
    }
}

void STFW_init(int _nInstances, unsigned char dual_enabled) {
    if (_nInstances <= 0) {
        printf("Error: number of instances need to be > 0\n");
        exit(-1);
    }


    int rank, nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    _DE = (dual_enabled) ? 1 : 0;

    _ninst = _nInstances;
    assert(_ninst > 0);

#ifdef D_SPEC_STFW
    struct stat st = {0};                                                                                                                                                                                   
    if (stat("./stfw_dbg_logs", &st) == -1) {
        mkdir("./stfw_dbg_logs", 0700);
    }


    sprintf(dbg_fname, "stfw_dbg_logs/outfile-%d-%d", nprocs, rank);

    sdbgfp = fopen(dbg_fname, "w");
    fprintf(sdbgfp, "> in STFW first init, instances = %d\n", _nInstances);
#endif
    gmG = malloc(_ninst * sizeof(*gmG));
    commG = malloc(_ninst * sizeof(*commG));
}

/******************************************************************************
 * Function:
 * Description:
 * @args:
 *   _instance_id   :   id of VPT instance
 *   vpt_ndims      :   number of VPT dims
 *   vpt_dsizes     :   sizes of each dim
 *   npsend         :   #processors in my send list
 *   sendlist       :   processor IDs in my send list
 *   nprecv         :   #prcs in my recv list
 *   recvlist       : procs IDs in my recv list
 *   ssend          : size of message to be sent to sendproc[i]
 *   srecv          : size of message to be recvd from proc[i]
 *   sendp          : of size npsend, contains pointers to send buffers for each
 *p in sendlist recvp          : of size nprecv, contains pointer to recv
 *buffers for p in recvlist
 *
 * Where:
 * Return:
 * Error:
 *****************************************************************************/
void STFW_init_instance(int _instance_id, int vpt_ndims, int *vpt_dsizes,
        int npsend, int *sendlist, int nprecv, int *recvlist,
        int *ssend, vptcomm_real_t **sendp, int *srecv,
        vptcomm_real_t **recvp, int *topomap) {

    if (_ninst == 0) {
        fprintf(stderr, "Error: you need to call STFW_init firsrt\n");
        exit(1);
    }
    vpt *gm;
    Comm *comm;
    gm = &gmG[_instance_id];
    comm = &commG[_instance_id];
    comm->commType = STFW_NAIVE;
    comm->send_unit_ptr = NULL;
    gm->nd = vpt_ndims;
    gm->size = malloc(sizeof(*gm->size) * gm->nd);
    memcpy(gm->size, vpt_dsizes, sizeof(*gm->size) * gm->nd);

    assert(sizeof(vptcomm_real_t) == sizeof(**sendp));

#ifdef D_SPEC_STFW
    fprintf(sdbgfp, "> hello from STFW init instance %d , total instances = %d\n",
            _instance_id, _ninst);
    fflush(sdbgfp);
#endif
    init_vpt(gm, D_LINEAR, topomap);

#ifdef NA_DBG
    na_log(sdbgfp, "dbg p2.0.2 after vpt init\n");
#endif
    /* construct comm structure from args
     *
     * */
    // Comm comm;
    comm->sno = npsend;
    comm->rno = nprecv;
    comm->rprocs = recvlist;
    comm->sprocs = sendlist;
    comm->ssizes = ssend;
    comm->rsizes = srecv;
    comm->send = sendp;
    comm->recv = recvp;

    init_stfw(comm, gm, DRYRUN_SPARSE_P2P);

#ifdef NA_DBG
    na_log(sdbgfp, "dbg p2.0.3 after stfw instance %d init\n", _instance_id);
#endif
}

/******************************************************************************
 * Function:
 * Description:
 * @args:
 *   _instance_id   :   id of VPT instance
 *   vpt_ndims      :   number of VPT dims
 *   vpt_dsizes     :   sizes of each dim
 *   npsend         :   #processors in my send list
 *   sendlist       :   processor IDs in my send list
 *   nprecv         :   #prcs in my recv list
 *   recvlist       : procs IDs in my recv list
 *   ssend          : size of message to be sent to sendproc[i]
 *   srecv          : size of message to be recvd from proc[i]
 *   sendp          : of size npsend, contains pointers to send buffers for each
 *p in sendlist recvp          : of size nprecv, contains pointer to recv
 *buffers for p in recvlist
 *
 * Where:
 * Return:
 * Error:
 *****************************************************************************/
void STFW_init_ExpRed_instance(int _instance_id, int vpt_ndims, int *vpt_dsizes,
        int npsend, int *sendlist, int nprecv, int *recvlist,
        int *ssend, vptcomm_real_t **sendp, int *srecv,
        vptcomm_real_t **recvp, int **sendIndsPerMsg, int **recvIndsPerMsg, int *indsmap, int unitSize, int maxUnitID, int *topomap) {

    if (_ninst == 0) {
        fprintf(stderr, "Error: you need to call STFW_init firsrt\n");
        exit(1);
    }
    if (unitSize <= 0) {
        fprintf(stderr, "Error: Unit size must be > 0\n");
        exit(1);
    }
    vpt *gm;
    Comm *comm;
    gm = &gmG[_instance_id];
    comm = &commG[_instance_id];
    comm ->commType = STFW_EXPRED;
    comm->send_unit_ptr = NULL;
    gm->nd = vpt_ndims;
    gm->size = malloc(sizeof(*gm->size) * gm->nd);
    memcpy(gm->size, vpt_dsizes, sizeof(*gm->size) * gm->nd);

    assert(sizeof(vptcomm_real_t) == sizeof(**sendp));

#ifdef D_SPEC_STFW
    fprintf(sdbgfp, "> hello from STFW init instance %d , total instances = %d\n",
            _instance_id, _ninst);
    fflush(sdbgfp);
#endif
    init_vpt(gm, D_LINEARREV, topomap);

#ifdef NA_DBG
    na_log(sdbgfp, "dbg p2.0.2 after vpt init\n");
#endif
    /* construct comm structure from args
     *
     * */
    // Comm comm;
    comm->sno = npsend;
    comm->rno = nprecv;
    comm->rprocs = recvlist;
    comm->sprocs = sendlist;
    comm->ssizes = ssend;
    comm->rsizes = srecv;
    comm->send = sendp;
    comm->recv = recvp;
    comm->recvinds = recvIndsPerMsg;
    comm->sendinds = sendIndsPerMsg;
    comm->unitSize = unitSize;
    comm->maxunitID = maxUnitID;
    comm->LGM = indsmap;
    

    init_stfw_ExpRed(comm, gm);

#ifdef NA_DBG
    na_log(sdbgfp, "dbg p2.0.3 after stfw instance %d init\n", _instance_id);
#endif
}

/*
 * ===  FUNCTION
 * ====================================================================== Name:
 * _invert_comm Description:
 * =====================================================================================
 */
void _invert_comm(Comm *comm) {

    swap(comm->rno, comm->sno);
    swap(comm->rprocs, comm->sprocs);
    swap(comm->ssizes, comm->rsizes);
    swap(comm->send, comm->recv);
} /* -----  end of function _invert_comm  ----- */

void _invert_vpt(vpt *gm) {
    int i, tmp, n;
    n = gm->nd - 1;

    /* reverse the order array */
    for (i = 0; i < gm->nd / 2; ++i) {
        /* tmp = gm->order[n];
           gm->order[n--] = gm->order[i];
           gm->order[i] = tmp; */
        swap(gm->order[n], gm->order[i]);
        --n;
    }

    /* fix comm_dim according to the new order */
    init_comm_dim(gm);

    /* NABIL: Assuming stfw algorithm with headers
     * TODO: implement inversion for no-header version, the pointers in pck_cp
     * should be updated */
    for (i = 0; i < gm->nd; ++i) {
        stfw_comm *sc = &gm->sc[i];

        swap(sc->fw_buf, sc->st_buf);
        swap(sc->fw_np, sc->st_np);
        swap(sc->fw_procs, sc->st_procs);
        swap(sc->fw_imap, sc->fw_imap_d);
        swap(sc->fw_locs, sc->st_locs);
        swap(sc->fw_offset, sc->fw_offset_d);
        swap(sc->st_reqs, sc->st_reqs_d);
        swap(sc->st_stts, sc->st_stts_d);
    }
}

void STFW_Comm(int _instance_id) {

    vpt *gm = &gmG[_instance_id];
    Comm *comm = &commG[_instance_id];
    if(comm->commType != STFW_NAIVE){ fprintf(stderr,"Error: Use STFW_Comm_Expand or STFW_Comm_Reduce for this instance\n"); goto ERR;}

#ifdef D_SPEC_STFW
    fprintf(sdbgfp, "starting stfw_nh for instance %d mode %d\n", _instance_id,
            _instance_id / 2);
    fflush(sdbgfp);
#endif

    stfw_nh(comm, gm);

#ifdef D_SPEC_STFW
    fprintf(sdbgfp, "end of stfw_nh for instance %d mode %d\n", _instance_id,
            _instance_id / 2);
    fflush(sdbgfp);
#endif
    return;
ERR:
    MPI_Finalize();
    exit(1);
}

void STFW_Comm_Expand(int _instance_id) {

    vpt *gm = &gmG[_instance_id];
    Comm *comm = &commG[_instance_id];
    if(comm->commType != STFW_EXPRED){ printf("Error: Use STFW_Comm for this instance\n"); goto ERR;}
#ifdef D_SPEC_STFW
    fprintf(sdbgfp, "starting stfw Expand for instance %d\n", _instance_id);
    fflush(sdbgfp);
#endif

    stfw_ExpRed(comm, gm, NULL, NULL, SPOP_EXPAND);

#ifdef D_SPEC_STFW
    MPI_Barrier(MPI_COMM_WORLD);
    fprintf(sdbgfp, "end of stfw Expand for instance %d\n", _instance_id);
    fflush(sdbgfp);
#endif
    return;
ERR:
    MPI_Finalize();
    exit(1);
}

void STFW_Comm_Reduce(int _instance_id, void (*reduce_fun)(vptcomm_real_t *, vptcomm_real_t *, int), void (*reduce_init_fun)(vptcomm_real_t *, int)) {

    vpt *gm = &gmG[_instance_id];
    Comm *comm = &commG[_instance_id];
    if(comm->commType != STFW_EXPRED){ printf("Error: Use STFW_Comm for this instance\n"); goto ERR;}
#ifdef D_SPEC_STFW
    fprintf(sdbgfp, "starting stfw Expand for instance %d mode %d\n", _instance_id,
            _instance_id / 2);
    fflush(sdbgfp);
#endif

    stfw_ExpRed(comm, gm, reduce_fun, reduce_init_fun, SPOP_REDUCE);

#ifdef D_SPEC_STFW
    fprintf(sdbgfp, "end of stfw Expand for instance %d mode %d\n", _instance_id,
            _instance_id / 2);
    fflush(sdbgfp);
#endif
    return;
ERR:
    MPI_Finalize();
    exit(1);
}

void STFW_dual_Comm(int _instance_id) {

    if (!_DE) {
        fprintf(stderr, "Error: dual is not enabled, you should enable it with "
                "STFW_init\n Exitting\n");
        exit(-1);
    }
#ifdef NA_DBG
    fprintf(sdbgfp, "stfw_run_dbg: > start of dual comm\n");
    fflush(sdbgfp);
#endif
    vpt *gm = &gmG[_instance_id];
    Comm *comm = &commG[_instance_id];
    int i, j;

    _invert_vpt(gm);
    _invert_comm(comm);

#ifdef NA_DBG
    fprintf(sdbgfp, "stfw_run_dbg: > after first swaps\n");
    fflush(sdbgfp);
#endif

    /* set the "new" fw_buffs to zero */
    for (i = 0; i < gm->nd; ++i) {
        for (j = 0; j < gm->sc[i].fw_np; ++j) {
            gm->sc[i].fw_buf[gm->sc[i].fw_locs[j]] = 0.0;
        }
    }
    // stfw_na(comm, gm);
    //init_stfw_nh(comm, gm);

#ifdef NA_DBG
    fprintf(sdbgfp, "stfw_run_dbg: > after dual comm \n");
    fflush(sdbgfp);
#endif
    _invert_vpt(gm);
    _invert_comm(comm);

    /* set the "old" fw_buffs to zero again*/
    for (i = 0; i < gm->nd; ++i) {
        for (j = 0; j < gm->sc[i].fw_np; ++j) {
            gm->sc[i].fw_buf[gm->sc[i].fw_locs[j]] = 0.0;
        }
    }
#ifdef NA_DBG
    fprintf(sdbgfp, "stfw_run_dbg: > after second swap\n");
    fflush(sdbgfp);
#endif
}

void STFW_finalize() {
    int i;
#ifdef D_SPEC_STFW
    MPI_Barrier(MPI_COMM_WORLD);
    fprintf(sdbgfp, "starting STFW_finalize\n");
    fflush(sdbgfp);
#endif
    for (i = 0; i < _ninst; ++i) {
        free_vpt(&(gmG[i]));
        if(commG[i].send_unit_ptr) free(commG[i].send_unit_ptr);
#ifdef D_SPEC_STFW
    MPI_Barrier(MPI_COMM_WORLD);
    fprintf(sdbgfp, "STFW_finalize of instance %d done\n", i);
    fflush(sdbgfp);
#endif
    }
    if(_ninst){ free(gmG); free(commG);}
    _ninst = 0;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  STFW_stats
 *  Description:  gets actual communication stats using STFW VPT
 * =====================================================================================
 */
void STFW_stats(int _instance_id, int *maxSendVol, int *maxRecvVol, int *totalVol , int *maxSendMsgs, int *maxRecvMsgs, int *totalMsgs, int factor){
    gm_cstats_na(&gmG[_instance_id], maxSendVol, maxRecvVol, totalVol, maxSendMsgs, maxRecvMsgs, totalMsgs, factor);
#ifdef D_SPEC_STFW
    MPI_Barrier(MPI_COMM_WORLD);
    fprintf(sdbgfp, "end of gm_cstats_na for instance %d\n", _instance_id);
    fflush(sdbgfp);
#endif
}


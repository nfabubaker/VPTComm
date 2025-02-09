/**
 * @author      : Nabil Abubaker (nabil.abubaker@bilkent.edu.tr)
 * @file        : docomm
 * @created     : Pazartesi Tem 20, 2020 14:30:11 +03
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mpi.h"
#include "vptcomm.h"
#include "stfw.h"
#include "util.h"

void stfw_ExpRed(Comm *comm, vpt *gm ,  void (*reduce_op)(vptcomm_real_t*, vptcomm_real_t*, int) , void (*reduce_op_init)(vptcomm_real_t *, int), enum SparseOp SpOp) {
#define TAGA 77
#define TAGB 777

    int i, j, k;
    stfw_comm *cm;
    int firstd = gm->order[0];
    int unitSize = comm->unitSize;

    /* first, init fw_buffs */
    if(SpOp == SPOP_REDUCE){ for (i = 0; i < gm->nd; ++i) { cm = &gm->sc[i]; if(cm->fw_np) reduce_op_init(cm->fw_buf, cm->fw_locs[cm->fw_np]); } }

    cm = &gm->sc[firstd];
    for (i = 0; i < cm->st_np; ++i)
        MPI_Irecv(&cm->st_buf[cm->st_locs[i]], cm->st_locs[i + 1] - cm->st_locs[i],
                vptcomm_mpi_real_t, gm->inv_map[cm->st_procs[i]], TAGA, MPI_COMM_WORLD,
                &cm->st_reqs[i]);

    /* copy elements of p to corresponding locations */
    k = 0;
    for (i = 0; i < comm->sno; ++i) {
        int proc = gm->map[comm->sprocs[i]];
        vptcomm_real_t *tptr = comm->send[i];
        /* copy each unit to a sepcific location in fw_buf, useful in Expand */
#ifdef D_SPEC_STFW_L3
            na_log(sdbgfp, "copying");
#endif
        for(j = 0; j < comm->ssizes[i] / unitSize; ++j){
#ifdef D_SPEC_STFW_L3
            na_log(sdbgfp, " %f %f %f %f", tptr[0], tptr[1], tptr[2], tptr[3]);
#endif
            memcpy(comm->send_unit_ptr[k++], tptr, sizeof(vptcomm_real_t) * unitSize);
            tptr += unitSize;
        }
#ifdef D_SPEC_STFW_L3
         na_log(sdbgfp, " to sent ptr\n");
#endif
    }
#ifdef D_SPEC_STFW
            na_log(sdbgfp, "Irecvs of 1st dim issued\n");
#endif

    /* first dimension SENDs */
    cm = &gm->sc[firstd];
    for (i = 0; i < cm->fw_np; ++i) {
#ifdef D_SPEC_STFW_L3
            na_log(sdbgfp, "sending");
#endif
            vptcomm_real_t *tptr = &cm->fw_buf[cm->fw_locs[i]]; int ff=cm->fw_locs[i];
#ifdef D_SPEC_STFW_L3
            while(ff<cm->fw_locs[i+1]){
                na_log(sdbgfp, " %f", *tptr ); tptr+=unitSize; ff+=unitSize;}
#endif
        MPI_Send(&cm->fw_buf[cm->fw_locs[i]],
                cm->fw_locs[i + 1] - cm->fw_locs[i], vptcomm_mpi_real_t,
                gm->inv_map[cm->fw_procs[i]], TAGA, MPI_COMM_WORLD);
#ifdef D_SPEC_STFW_L3
        na_log(sdbgfp, " to %d at dim %d\n", gm->inv_map[cm->fw_procs[i]], firstd);
#endif

    }
#ifdef D_SPEC_STFW
            na_log(sdbgfp, "Sends of 1st dim Done\n");
#endif

    /* loop starting from second dim */
    for (i = 1; i < gm->nd; ++i) {
        int cur_dim = gm->order[i];
        cm = &gm->sc[cur_dim];

        /* RECVs of this dim */
        for (j = 0; j < cm->st_np; ++j)
            MPI_Irecv(&cm->st_buf[cm->st_locs[j]],
                    cm->st_locs[j + 1] - cm->st_locs[j], vptcomm_mpi_real_t,
                    gm->inv_map[cm->st_procs[j]], TAGB, MPI_COMM_WORLD,
                    &cm->st_reqs[j]);

        /* wait the prev dim's RECVs and process them */
        int prev_dim = gm->order[i - 1];
        cm = &gm->sc[prev_dim];
        int nrecv = 0;
        int recvid;
        MPI_Status stts;
        while (nrecv < cm->st_np) {
            MPI_Waitany(cm->st_np, cm->st_reqs, &recvid, &stts);
            if (recvid != MPI_UNDEFINED) {
                cm->st_reqs[recvid] = MPI_REQUEST_NULL;
                vptcomm_real_t *ptr = &cm->st_buf[cm->st_locs[recvid]];
                int starti = cm->st_locs[recvid] / unitSize; int endi=cm->st_locs[recvid+1]/unitSize;
#ifdef D_SPEC_STFW_L3
                na_log(sdbgfp, "recvd");
#endif
                for (k = starti; k< endi; ++k) {
#ifdef D_SPEC_STFW_L3
                    na_log(sdbgfp, " %f %f %f %f", ptr[0], ptr[1], ptr[2], ptr[3]);
#endif
                    if(SpOp == SPOP_EXPAND){
                        int l; for (l = cm->unit_RepCnt[k]; l < cm->unit_RepCnt[k+1]; ++l) {
                            memcpy(cm->unit_ptr[l], ptr, unitSize * sizeof(vptcomm_real_t));  } 
                    }
                    else{ /* REDUCE */
                        //if(!cm->unit_tmp_counter[k]){ reduce_op(cm->unit_ptr[k], ptr, unitSize); cm->unit_tmp_counter[k]++; }
                        //else{ memcpy(cm->unit_ptr[k], ptr, unitSize * sizeof(vptcomm_real_t)); cm->unit_tmp_counter[k] = 1; }
                        reduce_op(cm->unit_ptr[k], ptr, unitSize);
                    }
                    ptr += unitSize;
                }
#ifdef D_SPEC_STFW_L3
                na_log(sdbgfp, " from %d at dim %d\n", gm->inv_map[cm->st_procs[recvid]], prev_dim);
#endif
                ++nrecv;
            }
        }
#ifdef D_SPEC_STFW
            na_log(sdbgfp, "recvs of dim %d Done\n", prev_dim);
#endif

        /* SENDs of this dim */
        cm = &gm->sc[cur_dim];
        for (j = 0; j < cm->fw_np; ++j) {
#ifdef D_SPEC_STFW_L3
            na_log(sdbgfp, "sending");
#endif
            vptcomm_real_t *tptr = &cm->fw_buf[cm->fw_locs[j]]; int ff=cm->fw_locs[j];
#ifdef D_SPEC_STFW_L3
            while(ff<cm->fw_locs[j+1]){
                na_log(sdbgfp, " %f", *tptr ); tptr+=unitSize; ff+=unitSize;}
#endif
            MPI_Send(&cm->fw_buf[cm->fw_locs[j]],
                    cm->fw_locs[j + 1] - cm->fw_locs[j], vptcomm_mpi_real_t,
                    gm->inv_map[cm->fw_procs[j]], TAGB, MPI_COMM_WORLD);
#ifdef D_SPEC_STFW_L3
            na_log(sdbgfp, " to %d at dim %d\n", gm->inv_map[cm->fw_procs[j]], cur_dim);
#endif
        }
#ifdef D_SPEC_STFW
            na_log(sdbgfp, "sends of dim %d Done\n", cur_dim);
#endif
    }

    /* wait the last dim's RECVs and process them */
    int prev_dim = gm->order[i - 1];
    cm = &gm->sc[prev_dim];
    int nrecv = 0;
    int recvid;
    MPI_Status stts;
    while (nrecv < cm->st_np) {
        MPI_Waitany(cm->st_np, cm->st_reqs, &recvid, &stts);
        if (recvid != MPI_UNDEFINED) {
            cm->st_reqs[recvid] = MPI_REQUEST_NULL;
            vptcomm_real_t *ptr = &cm->st_buf[cm->st_locs[recvid]];
#ifdef D_SPEC_STFW_L3
            na_log(sdbgfp,"recvd "); 
#endif
            int starti = cm->st_locs[recvid] / unitSize; int endi=cm->st_locs[recvid+1]/unitSize;
            for (k = starti; k <endi; ++k) { 
#ifdef D_SPEC_STFW_L3
                na_log(sdbgfp, " %f %f %f %f", ptr[0], ptr[1], ptr[2], ptr[3]);
#endif
                memcpy(cm->unit_ptr[k], ptr, unitSize * sizeof(vptcomm_real_t)); ptr += unitSize;}
            ++nrecv;
#ifdef D_SPEC_STFW_L3
            na_log(sdbgfp, " from %d at dim %d\n", gm->inv_map[cm->st_procs[recvid]], prev_dim);
#endif
        }
    }
#ifdef D_SPEC_STFW
            na_log(sdbgfp, "recvs of dim %d Done\n", prev_dim);
#endif
}


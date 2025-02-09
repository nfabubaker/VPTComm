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

void stfw_nh(Comm *comm, vpt *gm
        /* vptcomm_real_t	*p */
        ) {
#define TAGA 77
#define TAGB 777

    int i, j, k;
    stfw_comm *cm;
    int firstd = gm->order[0];

    cm = &gm->sc[firstd];
    for (i = 0; i < cm->st_np; ++i)
        MPI_Irecv(&cm->st_buf[cm->st_locs[i]], cm->st_locs[i + 1] - cm->st_locs[i],
                vptcomm_mpi_real_t, gm->inv_map[cm->st_procs[i]], TAGA, MPI_COMM_WORLD,
                &cm->st_reqs[i]);

    /* copy elements of p to corresponding locations */
    for (i = 0; i < comm->sno; ++i) {
        int proc = gm->map[comm->sprocs[i]];
        int d = gm->comm_dim[proc];
        cm = &gm->sc[d];
        int proc_crd = get_crd(gm, proc, d);
        int iproc = cm->fw_imap[proc_crd];

#ifdef D_SPEC_STFW_L3
            na_log(sdbgfp, "copying");
            for(j=0; j<comm->ssizes[i]; ++j) 
                na_log(sdbgfp, " %0.2f", *((comm->send[i])+j));
#endif
        memcpy(&(cm->fw_buf[cm->fw_offset[iproc]]), comm->send[i], comm->ssizes[i] * sizeof(vptcomm_real_t));
        cm->fw_offset[iproc]+=comm->ssizes[i];
#ifdef D_SPEC_STFW_L3
         na_log(sdbgfp, " to sent ptr\n");
#endif
    }

    /* first dimension SENDs */
    cm = &gm->sc[firstd];
    for (i = 0; i < cm->fw_np; ++i) {
#ifdef D_SPEC_STFW_L3
            na_log(sdbgfp, "sending");
            for(k=0; k<cm->fw_locs[i + 1] - cm->fw_locs[i]; ++k) 
                na_log(sdbgfp, " %0.2f", cm->fw_buf[cm->fw_locs[i]+k]);
#endif
        MPI_Send(&cm->fw_buf[cm->fw_locs[i]],
                cm->fw_locs[i + 1] - cm->fw_locs[i], vptcomm_mpi_real_t,
                gm->inv_map[cm->fw_procs[i]], TAGA, MPI_COMM_WORLD);

#ifdef D_SPEC_STFW_L3
        na_log(sdbgfp, " to %d at dim %d\n", gm->inv_map[cm->fw_procs[i]], firstd);
#endif
        /* need to reset fw_offset */
        cm->fw_offset[i] = cm->fw_locs[i];
    }

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
#ifdef D_SPEC_STFW_L3
                na_log(sdbgfp, "recvd");
#endif
                for (k = cm->pck_cnt[recvid]; k < cm->pck_cnt[recvid + 1];
                        ++k) {
#ifdef D_SPEC_STFW_L3
                    int z; for(z=0; z < cm->pck_sizes[k]; ++z)
                        na_log(sdbgfp, " %0.2f", *(ptr+z));
#endif
                    memcpy(cm->pck_cp[k], ptr, cm->pck_sizes[k] * sizeof(vptcomm_real_t));
                    ptr += cm->pck_sizes[k];
                }
#ifdef D_SPEC_STFW_L3
                na_log(sdbgfp, " from %d at dim %d\n", gm->inv_map[cm->st_procs[recvid]], prev_dim);
#endif
                ++nrecv;
            }
        }

        /* SENDs of this dim */
        cm = &gm->sc[cur_dim];
        for (j = 0; j < cm->fw_np; ++j) {
#ifdef D_SPEC_STFW_L3
            na_log(sdbgfp, "sending");
            for(k=0; k<cm->fw_locs[j + 1] - cm->fw_locs[j]; ++k) 
                na_log(sdbgfp, " %0.2f", cm->fw_buf[cm->fw_locs[j]+k]);
#endif
            MPI_Send(&cm->fw_buf[cm->fw_locs[j]],
                    cm->fw_locs[j + 1] - cm->fw_locs[j], vptcomm_mpi_real_t,
                    gm->inv_map[cm->fw_procs[j]], TAGB, MPI_COMM_WORLD);
            /* need to reset fw_offset and fw_buf npackets fields for later
            */
            cm->fw_offset[j] = cm->fw_locs[j];
#ifdef D_SPEC_STFW_L3
        na_log(sdbgfp, " to %d at dim %d\n", gm->inv_map[cm->fw_procs[j]], cur_dim);
#endif
        }
    }

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
            for (k = cm->pck_cnt[recvid]; k < cm->pck_cnt[recvid + 1];
                    ++k) {
                memcpy(cm->pck_cp[k], ptr, cm->pck_sizes[k] * sizeof(vptcomm_real_t));
                ptr += cm->pck_sizes[k];
            }
            ++nrecv;
        }
    }
}

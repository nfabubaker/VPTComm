//include "dbg_print.h"
#include "vptcomm.h"
#include "stfw_def.h"
#include "util.h"
#include "vpt.h"
#include "stfw.h"
#include <assert.h>
#include <mpi.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
/* #define DEBUG */


typedef struct _unit{
    int src;
    int dst;
    int id;
    struct _unit *next;
} unit;


typedef struct _setupHelper{
    int **sbuff;
    int *sbuffSizes;
    int **rbuff;
    int *rbuffSizes;
    unit* *fw_list;
    /* gids of units in current fw_buff */
    int **fwbuff_gids;
    int *fw_GtoL;
    MPI_Request *reqs;
} setupH;

/* setup helper for comm */
typedef struct _commSH{
    int *inv_send_map;
    int *inv_recv_map;
    int **sendindsSorted;
    int **recvindsSorted;
    int **sendindsLocs;
    int **recvindsLocs;
} commSH;

void init_SH(setupH *sh){
    sh->sbuff = NULL; sh->sbuffSizes = NULL; sh->rbuff = NULL; sh->rbuffSizes = NULL;
    sh->fw_list = NULL; sh->fwbuff_gids = NULL; sh->fw_GtoL = NULL; sh->reqs = NULL;
}

void init_CSH(commSH *csh){
    csh->inv_recv_map = csh->inv_send_map = NULL;
    csh->recvindsLocs = csh->sendindsLocs = NULL;
    csh->sendindsSorted = csh->recvindsSorted = NULL;
}
void free_SH(setupH *sh, int dimsize){
    int i;
    for (i = 0; i < dimsize; ++i) {
        if(sh->fwbuff_gids[i]) free(sh->fwbuff_gids[i]);
        if(sh->rbuff[i])free(sh->rbuff[i]);
    }
    free(sh->rbuff); free(sh->rbuffSizes); free(sh->fw_list); free(sh->fw_GtoL); free(sh->reqs); free(sh->fwbuff_gids);
}

void free_CSH(commSH *csh, Comm *comm){
    int i;
    for (i = 0; i < comm->sno; ++i) {
        if(csh->sendindsLocs) if(csh->sendindsLocs[i]) free(csh->sendindsLocs[i]); 
        if(csh->sendindsSorted) if(csh->sendindsSorted[i]) free(csh->sendindsSorted[i]); 
    }
    for (i = 0; i < comm->rno; ++i) {
        if(csh->recvindsLocs) if(csh->recvindsLocs[i]) free(csh->recvindsLocs[i]); 
        if(csh->recvindsSorted) if(csh->recvindsSorted[i]) free(csh->recvindsSorted[i]); 
    }
    if(csh->sendindsLocs) free(csh->sendindsLocs); if(csh->sendindsSorted) free(csh->sendindsSorted);
    if(csh->recvindsLocs) free(csh->recvindsLocs); if(csh->recvindsSorted) free(csh->recvindsSorted); 
    if(csh->inv_send_map) free(csh->inv_send_map); if(csh->inv_recv_map) free(csh->inv_recv_map);
}

int compareints (const void * a, const void * b){ return ( *(int*)a - *(int*)b );}

void add_to_fw_list(unit **head, int src, int dst, int id){
    unit *u = malloc(sizeof(*u)); u->src=src; u->dst=dst; u->id=id; u->next= NULL; u->next = *head; *head = u;
}

void remove_from_fw_list(unit **head){ unit *tu; tu = *head; *head = (*head)->next; free(tu); }

/* should produce send buffers of the form
 * (id, repFactor, src, dst1, dst2, dst3, id,repFactor,srs,dst1, dst2, dst3, dst4 ...)
 *
 * Params:
 * [IN] int curd: current dimension
 * [IN] vpt* gm: pointer to the VPT
 * [IN] unit **fw_lists: pointer to the fw_list of this dim
 * [OUT] int ***cur_sbuff: pointer to the send buffers per neighbor of this dim
 * [OUT] int **cur_sbuffSizes: pointer to the send buffer sizes per neighbor of this dim
 * [OUT] int **cur_sbuffSizes: pointer to the send buffer sizes per neighbor of this dim
 *
 * */
void process_fw_list(int curd, vpt *gm, setupH *curr_sh, int unitSize,  int maxunitID){
    int i,j,k, dimsize = gm->size[curd], *perUnitCnt, *globalToLocal, *perCrdCnt, **cur_sbuff;
    int *uniqUnitCnt;
    perUnitCnt = malloc(sizeof(*perUnitCnt) * maxunitID);
    globalToLocal = malloc(sizeof(*globalToLocal) * maxunitID);
    curr_sh->sbuffSizes = calloc( dimsize, sizeof(*curr_sh->sbuffSizes));
    curr_sh->fw_GtoL = calloc( dimsize, sizeof(*curr_sh->fw_GtoL));
    perCrdCnt = curr_sh->sbuffSizes;
    uniqUnitCnt = calloc( dimsize, sizeof(*uniqUnitCnt));
    curr_sh->sbuff = malloc(sizeof(*curr_sh->sbuff) * dimsize);
    cur_sbuff = curr_sh->sbuff;
    curr_sh->fwbuff_gids = malloc(sizeof(*curr_sh->fwbuff_gids) * dimsize);
    gm->sc[curd].fw_np = 0;
    for(i = 0; i < dimsize; ++i){
        /* 1- initialize perUnitCnt and globalToLocal map for this coordinate */
        for(j =0; j < maxunitID; ++j){ perUnitCnt[j] = 0; globalToLocal[j]=-1;}
        /* 2 - Count number of unique occurances per unit  */
        unit *tu = curr_sh->fw_list[i];
        while(tu != NULL){ perUnitCnt[tu->id]++; tu = tu->next; }
        /* 3- Calculate the total send buff size of this coord, and # of uniqUnits to be sent */
        uniqUnitCnt[i] = 0;
        for(j =0; j < maxunitID; ++j) if( perUnitCnt[j] > 0){
            perCrdCnt[i]+= (3 + perUnitCnt[j]); 
            globalToLocal[j] = uniqUnitCnt[i]++;
        }
        if(uniqUnitCnt[i] > 0) gm->sc[curd].fw_np +=1;
        /* 4- allocate send buffer for this coord, and associated helpers */
        int *sbuff_offsets = calloc(uniqUnitCnt[i]+1, sizeof(*sbuff_offsets));
        char *sbuff_flags = calloc(uniqUnitCnt[i], sizeof(*sbuff_flags));
        cur_sbuff[i]=NULL; if(perCrdCnt[i] > 0) (cur_sbuff)[i] = malloc(sizeof(**cur_sbuff) * perCrdCnt[i]);
        curr_sh->fwbuff_gids[i]=NULL; if(uniqUnitCnt[i]) curr_sh->fwbuff_gids[i] = malloc(sizeof(*curr_sh->fwbuff_gids[i]) * uniqUnitCnt[i]);
        uniqUnitCnt[i] = 0;
        for(j =0; j < maxunitID; ++j) if( perUnitCnt[j] > 0)
        { sbuff_offsets[(uniqUnitCnt[i])+1] = (3 + perUnitCnt[j]); curr_sh->fwbuff_gids[i][uniqUnitCnt[i]++] = j;}
        for(j=1; j <= uniqUnitCnt[i]; ++j) sbuff_offsets[j] += sbuff_offsets[j-1];

        /*  5- Fill the send buffer for this coord */
        tu = curr_sh->fw_list[i];
        while (tu != NULL) {
            int idx = globalToLocal[tu->id];
            if(sbuff_flags[idx]){
                (cur_sbuff)[i][sbuff_offsets[idx]++] = tu->dst; }
            else{
                (cur_sbuff)[i][sbuff_offsets[idx]++] = tu->id;
                (cur_sbuff)[i][sbuff_offsets[idx]++] = perUnitCnt[tu->id];
                (cur_sbuff)[i][sbuff_offsets[idx]++] = tu->src;
                (cur_sbuff)[i][sbuff_offsets[idx]++] = tu->dst;
                sbuff_flags[idx] = 1;
            }
            unit* ttu = tu; tu = tu->next; free(ttu);
        }
        /* sanity check */
        if(uniqUnitCnt[i])
            assert(sbuff_offsets[uniqUnitCnt[i]-1] == sbuff_offsets[uniqUnitCnt[i]] && sbuff_offsets[uniqUnitCnt[i]] == perCrdCnt[i]);

        /* cleanup for this coord */
        free(sbuff_flags); free(sbuff_offsets);
    }

    /* setup fw info for this dim */
    stfw_comm *cm = &gm->sc[curd];
    cm->fw_buf = NULL;
    cm->fw_procs = NULL;
    cm->fw_locs = NULL;
    int total_fw_unit_cnt;
    total_fw_unit_cnt = 0;
    for(i = 0; i < dimsize; ++i){total_fw_unit_cnt += (uniqUnitCnt[i]); curr_sh->fw_GtoL[i] = -1;}
    if(total_fw_unit_cnt){ cm->fw_buf = malloc(total_fw_unit_cnt * unitSize * sizeof(*cm->fw_buf)); }
    if(cm->fw_np){ cm->fw_procs = malloc(sizeof(*cm->fw_procs) * cm->fw_np);}
    cm->fw_locs = malloc(sizeof(*cm->fw_locs) * (cm->fw_np+1));cm->fw_locs[0]=0; 
    int tcnt = 0; for(i = 0; i < dimsize; ++i) if(gm->crd[curd] != i && uniqUnitCnt[i] > 0)
    {
        curr_sh->fw_GtoL[i] = tcnt;
        cm->fw_procs[tcnt] = get_nghbr(gm, curd, i);
        cm->fw_locs[tcnt+1] = uniqUnitCnt[i]*unitSize; ++tcnt;
    } 
    for (i = 1; i <= cm->fw_np; ++i) cm->fw_locs[i] += cm->fw_locs[i-1];

    /* TODO FIXME : remove the following loop, no need unit_RepCnt should be for st_buff */
    /*     tcnt=0;
     *     for(i = 0; i < dimsize; ++i){
     *         for(j =0; j < maxunitID; ++j) perUnitCnt[j] = 0; 
     *         unit *tu = curr_sh->fw_list[i];
     *         while(tu != NULL){ perUnitCnt[tu->id]++; unit* ttu = tu; tu = tu->next; free(ttu); }
     *         for(j =0; j < maxunitID; ++j) if( perUnitCnt[j] > 0) cm->unit_RepCnt[tcnt++] = perUnitCnt[j];
     *     }
     */

    /* cleanup*/
    free(perUnitCnt); free(globalToLocal); free(uniqUnitCnt);
}

void process_recv_buff(int curd, vpt *gm, setupH *sh, int unitSize, int maxunitID){
    int i,j,k, dimsize = gm->size[curd];
    int rank, nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int *uniqUnitCnt = calloc( dimsize, sizeof(*uniqUnitCnt));
    int *cur_rbuffSizes = sh[curd].rbuffSizes; int **cur_rbuff = sh[curd].rbuff;

#ifdef D_SPEC_STFW
    MPI_Barrier(MPI_COMM_WORLD);
    fprintf(sdbgfp, "In process_recv_buff: pt0\n");
    fflush(sdbgfp);
#endif
    /* go over units in recv buffer and allocate st_buff and other stuff */
    gm->sc[curd].st_np = 0;
    for (i = 0; i < dimsize; ++i) {
        int j = 0; int *tp = cur_rbuff[i];
        if(cur_rbuffSizes[i] > 0) gm->sc[curd].st_np++;
        while(j < cur_rbuffSizes[i]){
#ifdef D_SPEC_STFW
            fprintf(sdbgfp, "process_recv_buff: rbuffSizes[%d]=%d *tp+1=%d\n",i, cur_rbuffSizes[i], *(tp+1));
            fflush(sdbgfp);
#endif
            uniqUnitCnt[i]++; j += ((*(tp+1)) + 3); tp += ((*(tp+1)) + 3); }
    }
#ifdef D_SPEC_STFW
    MPI_Barrier(MPI_COMM_WORLD);
    fprintf(sdbgfp, "process_recv_buff: pt1\n");
    fflush(sdbgfp);
#endif

    /* allocate store info for current dim */
    stfw_comm *cm = &gm->sc[curd];
    cm->st_buf = NULL;
    cm->st_procs = NULL;
    cm->st_locs = NULL;
    cm->st_reqs = NULL;
    cm->st_stts = NULL;
    cm->unit_RepCnt = NULL;
    cm->unit_ptr = NULL;
    int total_st_unit_cnt = 0;
    for(i = 0; i < dimsize; ++i) total_st_unit_cnt += (uniqUnitCnt[i]);
    if(total_st_unit_cnt){ cm->st_buf = malloc(total_st_unit_cnt * unitSize * sizeof(*cm->st_buf));
        cm->unit_RepCnt = malloc((1+total_st_unit_cnt) * sizeof(*cm->unit_RepCnt));
        for (i = 0; i <= total_st_unit_cnt; ++i) cm->unit_RepCnt[i] =0; }
    if(cm->st_np){ cm->st_procs = malloc(sizeof(*cm->st_procs) * cm->st_np);
        cm->st_stts = malloc(sizeof(*cm->st_stts) * cm->st_np);
        cm->st_reqs = malloc(sizeof(*cm->st_reqs) * cm->st_np); }
    cm->st_locs = malloc(sizeof(*cm->st_locs) * (cm->st_np+1)); cm->st_locs[0]=0; 
    int tcnt = 0; for(i = 0; i < dimsize; ++i) if(gm->crd[curd] != i && uniqUnitCnt[i] > 0)
    {cm->st_procs[tcnt] = get_nghbr(gm, curd, i); cm->st_locs[tcnt+1] = uniqUnitCnt[i]*unitSize; ++tcnt;} 
    for (i = 1; i <= cm->st_np; ++i) cm->st_locs[i]+= cm->st_locs[i-1];
#ifdef D_SPEC_STFW
    MPI_Barrier(MPI_COMM_WORLD);
    fprintf(sdbgfp, "process_recv_buff: pt2\n");
    fflush(sdbgfp);
#endif

    /* count repFactors and allocate pointers array */
    tcnt = 0;
    for (i = 0; i < dimsize; ++i) {
        int j = 0; int *tp = cur_rbuff[i];
        while(j < cur_rbuffSizes[i]){
            int id=*(tp++); int repFactor = *(tp++); int src=*(tp++);

#ifdef D_SPEC_STFW
            na_log(sdbgfp, "id=%d repFactor=%d src=%d\n", id, repFactor, src);
#endif
            for(k=0; k<repFactor; ++k){
                int dst = *(tp++);
                if(dst != gm->map[rank]){
                    int d = gm->comm_dim[dst];
#ifdef D_SPEC_STFW
                    na_log(sdbgfp, "dst=%d dim=%d\n", dst, d);
#endif
                    int proc_crd = get_crd(gm, dst, d);
                    add_to_fw_list(&sh[d].fw_list[proc_crd], src, dst, id);
                }
            }
            cm->unit_RepCnt[tcnt+1] = repFactor; tcnt++;
            j += (repFactor + 3); 
        }

    }
#ifdef D_SPEC_STFW
    MPI_Barrier(MPI_COMM_WORLD);
    fprintf(sdbgfp, "process_recv_buff: pt3\n");
    fflush(sdbgfp);
#endif

    for (i = 1; i <= total_st_unit_cnt; ++i) {cm->unit_RepCnt[i]+= cm->unit_RepCnt[i-1];}
    if(total_st_unit_cnt) cm->unit_ptr = malloc(sizeof(*cm->unit_ptr) * cm->unit_RepCnt[total_st_unit_cnt]);
    free(uniqUnitCnt);

}

void setup_commSH(vpt *gm, Comm *comm, commSH *csh ){
    int i,j,k, *perUnitLoc;
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
#ifdef D_SPEC_STFW
    MPI_Barrier(MPI_COMM_WORLD);
    na_log(sdbgfp," > In setup_commSH\n"); 
#endif
    perUnitLoc = malloc(sizeof(*perUnitLoc) * comm->maxunitID);
    csh->inv_recv_map = malloc(sizeof(*csh->inv_recv_map) * nprocs);
    //csh->inv_send_map = malloc(sizeof(*csh->inv_send_map) * nprocs);
    if(comm->rno){ csh->recvindsLocs = malloc(sizeof(*csh->recvindsLocs) * comm->rno);
        //csh->sendindsLocs = malloc(sizeof(*csh->sendindsLocs) * comm->sno);
    }
    if(comm->rno){ csh->recvindsSorted = malloc(sizeof(*csh->recvindsSorted) * comm->rno);
        //csh->sendindsSorted = malloc(sizeof(*csh->sendindsSorted) * comm->sno);
    }
#ifdef D_SPEC_STFW
    na_log(sdbgfp," > setup_commSH initial allocations done\n"); 
#endif
    /* fill inv_recv_map */
    for (i = 0; i < comm->rno; ++i){
        csh->inv_recv_map[gm->map[comm->rprocs[i]]] = i;
        int rsize = (comm->rsizes[i] / comm->unitSize);
#ifdef D_SPEC_STFW
        na_log(sdbgfp,"setup_commSH: rno=%d i=%d comm->rsizes[i]=%d unitSize=%d rsize =%d\n", comm->rno, i, comm->rsizes[i], comm->unitSize, rsize); 
#endif
        if(rsize > 0){
            csh->recvindsSorted[i] = malloc(sizeof(*csh->recvindsSorted[i]) * (rsize));
            csh->recvindsLocs[i] = malloc(sizeof(*csh->recvindsLocs[i]) * (rsize));
        }
    }


    /*     for (i = 0; i < comm->sno; ++i) {csh->inv_send_map[gm->map[comm->sprocs[i]]] = gm->map[i];
     *         int ssize = (comm->ssizes[i] / comm->unitSize);
     *         csh->sendindsSorted[i] = malloc(sizeof(*csh->sendindsSorted[i]) * (ssize));
     *         csh->sendindsLocs[i] = malloc(sizeof(*csh->sendindsLocs[i]) * (ssize)); }
     */
#ifdef D_SPEC_STFW
    na_log(sdbgfp," > setup_commSH pt 1.0\n"); 
#endif

    for (i = 0; i < comm->rno; ++i) {
        for (j = 0; j < comm->maxunitID; ++j) perUnitLoc[j] = -1; 
#ifdef D_SPEC_STFW
        fprintf(sdbgfp, "%d %d %d\n", comm->rno, comm->rsizes[i], comm->unitSize); fflush(sdbgfp);
#endif
        for (j = 0; j < comm->rsizes[i]/comm->unitSize; ++j) perUnitLoc[comm->LGM[comm->recvinds[i][j]]] = j; 
        int tcnt = 0;
        for (j = 0; j < comm->maxunitID; ++j) if(perUnitLoc[j] > -1){ csh->recvindsLocs[i][tcnt]=perUnitLoc[j];  csh->recvindsSorted[i][tcnt++] = j; }
    }

    free(perUnitLoc);
}

void init_stfw_nh_nbx_ExpRed(Comm *comm, vpt *gm) {
#define TAG_HDR 7
    int i, j, k, rank, nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* data sizes communicated with each proc in each dimension  */

    int max_buf_sz = -1;

    /* setup Helper */
    setupH *sh = malloc(sizeof(*sh) * gm->nd);
    for (i = 0; i < gm->nd; ++i) {
        int dim_size = gm->size[i];
        init_SH(&sh[i]);
        sh[i].fw_list = malloc(dim_size * sizeof(*sh[i].fw_list));
        for (j = 0; j < dim_size; ++j) {
            sh[i].fw_list[j] = NULL;
        }
        sh[i].reqs = malloc(sizeof(*sh[i].reqs) * dim_size);

    }

#ifdef D_SPEC_STFW
    fprintf(sdbgfp, "> init_stfw_gmesh\n");
    fprintf(sdbgfp, "number of processors to communicate = %d\n", comm->sno);
    for (i = 0; i < comm->sno; ++i) fprintf(sdbgfp, "  P %d, send size = %d\n", comm->sprocs[i], comm->ssizes[i]);
    for (i = 0; i < comm->rno; ++i) fprintf(sdbgfp, "  P %d, recv size = %d\n", comm->rprocs[i], comm->rsizes[i]);
    fprintf(sdbgfp, "\n");
    fflush(sdbgfp);
#endif


    commSH csh;
    init_CSH(&csh);
    setup_commSH(gm, comm, &csh);
    /* process my send list and fill in the communicators in the respective
     * dimension */
    for (i = 0; i < comm->sno; ++i) {
        int proc = gm->map[comm->sprocs[i]];
        int d = gm->comm_dim[proc];

        //gm->inv_recv_map[comm->sprocs[i]] = i;

#ifdef D_SPEC_STFW
        fprintf(sdbgfp, "will send data to P%d in comm dim %d, ", proc, d);
        fflush(sdbgfp);
#endif
        int proc_crd = get_crd(gm, proc, d);
        for(j = 0; j < comm->ssizes[i]/comm->unitSize; ++j)
            add_to_fw_list(&sh[d].fw_list[proc_crd], gm->map[rank], proc, comm->LGM[comm->sendinds[i][j]]);
    }

#ifdef D_SPEC_STFW
    //print_hdr(gm, hdr_send, hdr_recv);
#endif

    for (i = 0; i < gm->nd; ++i) {
        int curd = gm->order[i];
        int my_crd = gm->crd[curd];
        unit* *cur_fwList = sh[curd].fw_list;
        MPI_Request *cur_hreqs = sh[curd].reqs;
        /*  process cur_hsend  and allocate send buffers*/
        process_fw_list(curd, gm, &sh[curd], comm->unitSize, comm->maxunitID);

#ifdef D_SPEC_STFW
        MPI_Barrier(MPI_COMM_WORLD);
        fprintf(sdbgfp, "\nissuing sends...\n");
        fflush(sdbgfp);
#endif

        sh[curd].rbuffSizes = calloc(gm->size[curd], sizeof(*sh[curd].rbuffSizes) );
        sh[curd].rbuff = malloc(sizeof(*sh[curd].rbuff) * gm->size[curd]);
        for (j = 0; j < gm->size[curd]; ++j) sh[curd].rbuff[j] = NULL; 

        int nsendto=0;
        for (j = 0; j < gm->size[curd]; ++j) {
            if (j != my_crd && sh[curd].sbuffSizes[j] > 0) {
#ifdef D_SPEC_STFW
                fprintf(sdbgfp, "  sending to P%d (crd %d at dim %d), send size = %d\n", get_nghbr(gm, curd, j), j, curd, sh[curd].sbuffSizes[j]);
                fflush(sdbgfp);
#endif

                MPI_Issend(sh[curd].sbuff[j], sh[curd].sbuffSizes[j], MPI_INT, gm->inv_map[get_nghbr(gm, curd, j)], TAG_HDR, MPI_COMM_WORLD, &cur_hreqs[nsendto]);
                ++nsendto;
            }
        }

#ifdef D_SPEC_STFW
        MPI_Barrier(MPI_COMM_WORLD);
        fprintf(sdbgfp, "\n\n"); fprintf(sdbgfp, "comm operations turn %d, dimension %d, my_crd = %d, ", i, curd, my_crd);
        fflush(sdbgfp);
#endif


        /* process hdr_recv of the prev step to form send data of this step */
        int done=0, barrier_act=0, flag, tot_pck_cnt = 0;
        MPI_Status stts;
        MPI_Request ibrrqst; 
        while(!done){
            MPI_Iprobe(MPI_ANY_SOURCE, TAG_HDR, MPI_COMM_WORLD, &flag, &stts);
            if(flag){
                int source = stts.MPI_SOURCE, nrcvdhdrs; MPI_Get_count(&stts, MPI_INT, &nrcvdhdrs);
                int rcrd=get_crd(gm, gm->map[source], curd);
#ifdef D_SPEC_STFW
                fprintf(sdbgfp, "recvd data: turn %d, dimension %d, crd = %d, source=%d size=%d\n", i, curd, rcrd, source, nrcvdhdrs);
                fflush(sdbgfp);
#endif
                sh[curd].rbuffSizes[rcrd] = nrcvdhdrs; tot_pck_cnt+= nrcvdhdrs;
                sh[curd].rbuff[rcrd] = malloc(sizeof(*sh[curd].rbuff[rcrd]) * nrcvdhdrs);
                MPI_Recv(sh[curd].rbuff[rcrd], nrcvdhdrs, MPI_INT, source, TAG_HDR, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                /* update store communicators */
            }
            if(barrier_act){ int ibrflag; MPI_Test(&ibrrqst, &ibrflag, MPI_STATUS_IGNORE); if(ibrflag) done=1; }
            else{
                int testflag=0;
                MPI_Testall(nsendto, cur_hreqs, &testflag, MPI_STATUSES_IGNORE);
                if(testflag) {MPI_Ibarrier(MPI_COMM_WORLD, &ibrrqst); barrier_act=1;}
            }
        }

#ifdef D_SPEC_STFW
        MPI_Barrier(MPI_COMM_WORLD);
        fprintf(sdbgfp, "  total recved size %d\n", tot_pck_cnt);
        fflush(sdbgfp);
#endif


        /* free sendbuff of this dim */
        int dimsize = gm->size[curd];
        if(sh[curd].sbuff){
            for(j = 0; j < dimsize; ++j) if(sh[curd].sbuff[j]) free(sh[curd].sbuff[j]);
            free(sh[curd].sbuff);
        }
        if(sh[curd].sbuffSizes) free(sh[curd].sbuffSizes); 
#ifdef D_SPEC_STFW
        MPI_Barrier(MPI_COMM_WORLD);
        fprintf(sdbgfp, "free sbuffs done\n");
        fflush(sdbgfp);
#endif
        /* process received data */
        process_recv_buff(curd, gm, sh, comm->unitSize, comm->maxunitID);
#ifdef D_SPEC_STFW
        MPI_Barrier(MPI_COMM_WORLD);
        fprintf(sdbgfp, "  process_recv_buff done \n");
        fflush(sdbgfp);
#endif

    }




    /* finalize copy data */
#ifdef D_SPEC_STFW
    fprintf(sdbgfp, "\n\nlink data:\n");
    fflush(sdbgfp);
#endif


    /*  go over send list again */
    int tidx = 0;
    for (i = 0; i < comm->sno; ++i) for (j = 0 ; j < comm->ssizes[i]/comm->unitSize ; ++j) tidx++; 
    comm->send_unit_ptr = malloc(sizeof(*comm->send_unit_ptr)*tidx); tidx=0;
    for (i = 0; i < comm->sno; ++i) {
        int dst = gm->map[comm->sprocs[i]];
        int d = gm->comm_dim[dst];
        int proc_crd = get_crd(gm, dst, d);
        for (j = 0 ; j < comm->ssizes[i]/comm->unitSize ; ++j) {
            int id = comm->LGM[comm->sendinds[i][j]];
            int *base = sh[d].fwbuff_gids[proc_crd];
            int size = gm->sc[d].fw_locs[sh[d].fw_GtoL[proc_crd] + 1] - gm->sc[d].fw_locs[sh[d].fw_GtoL[proc_crd]];
#ifdef D_SPEC_STFW
            //na_log(sdbgfp, "size=%d base+1-base = %d id=%d dst=%d d=%d proc_crd=%d fw_GtoL=%d fw_locs=%d\n", size, *(&base+1)-base, id, dst, d, proc_crd,sh[d].fw_GtoL[proc_crd], gm->sc[d].fw_locs[sh[d].fw_GtoL[proc_crd]+1]);
#endif
            size /= comm->unitSize;
            int *res = (int*) bsearch(&id, base, size, sizeof(int), compareints) ;
            int loc = res - base;
#ifdef D_SPEC_STFW
            na_log(sdbgfp, "[IN dim %d] loc of %d src=%d dst=%d repFactor=%d is %d and w.r.t fw_buf %d \n ",d, id, gm->map[rank], dst, 1, loc, gm->sc[d].fw_locs[sh[d].fw_GtoL[proc_crd]] + (loc * comm->unitSize));
#endif
            comm->send_unit_ptr[tidx++] = &(gm->sc[d].fw_buf[gm->sc[d].fw_locs[sh[d].fw_GtoL[proc_crd]] + (loc * comm->unitSize)]);
        }
    }

    /* go over recv buffers again */
    for (i = 0; i < gm->nd; ++i) {
        stfw_comm *cm = &gm->sc[i];
        setupH *curr_sh = &sh[i];
        int dimsize = gm->size[i];
        int tidx = 0;
        for (j = 0; j < dimsize; ++j) {
            int z = 0; int *tp = curr_sh->rbuff[j];
            while(z < curr_sh->rbuffSizes[j]){
                int id=*(tp++); int repFactor = *(tp++); int src=*(tp++);
                for(k=0; k<repFactor; ++k){
                    int dst = *(tp++);
                    /* Important: if dst is the same as already existing dst, then
                     * this is reduce, and no need for repFactor to be more than 1
                     * THIS DOES NOT HAPPEN IN EXPAND*/
                    if(k==0 && repFactor > 1 && dst == (*tp)){ k+= repFactor; tp+=(repFactor-1);}
                    if(dst != gm->map[rank]){
                        int d = gm->comm_dim[dst];
                        int proc_crd = get_crd(gm, dst, d);
                        int *base = sh[d].fwbuff_gids[proc_crd];
                        int size = gm->sc[d].fw_locs[sh[d].fw_GtoL[proc_crd] + 1] - gm->sc[d].fw_locs[sh[d].fw_GtoL[proc_crd]];
                        size /= comm->unitSize;
                        int *res = (int *) bsearch(&id, base, size, sizeof(int), compareints) ;
                        int loc = res - base;
#ifdef D_SPEC_STFW
                        na_log(sdbgfp, "[FW dim %d] loc of %d src=%d dst=%d repFactor=%d is %d and w.r.t fw_buf %d \n ",d, id, src, dst, repFactor, loc, gm->sc[d].fw_locs[sh[d].fw_GtoL[proc_crd]] + (loc * comm->unitSize));
#endif
                        cm->unit_ptr[tidx++] = &(gm->sc[d].fw_buf[gm->sc[d].fw_locs[sh[d].fw_GtoL[proc_crd]] + (loc * comm->unitSize)]);
                    }
                    else{ /* I'm the receiver */
                        int lidx = csh.inv_recv_map[src];
                        int *base = csh.recvindsSorted[lidx];
                        int nmem = comm->rsizes[lidx]/comm->unitSize;
                        int *res = (int *) bsearch(&id, base, nmem, sizeof(int), compareints);
                        int loc = csh.recvindsLocs[lidx][res-base];
#ifdef D_SPEC_STFW
                        na_log(sdbgfp, "[RV dim %d] loc of %d src=%d dst=%d repFactor=%d is %d and w.r.t fw_buf %d \n ",i, id, src, dst, repFactor, loc, loc*comm->unitSize);
#endif
                        cm->unit_ptr[tidx++] = &(comm->recv[lidx][loc*comm->unitSize]);
                    }
                }
                z += (repFactor + 3); 
            }

        }
    }
#ifdef D_SPEC_STFW
    fprintf(sdbgfp, "\n\nlink data COMPLETED\n");
    fflush(sdbgfp);
#endif
    /* cleanup */
    for (i = 0; i < gm->nd; ++i) {
        int dimsize = gm->size[i];
        free_SH(&sh[i], dimsize);
    }
    free_CSH(&csh,comm); free(sh);

    //gm_cstats(gm);

    /*     int max_buf_sz_reduced;
     *     MPI_Reduce(&max_buf_sz, &max_buf_sz_reduced, 1, MPI_INT, MPI_MAX, 0,
     *             MPI_COMM_WORLD);
     *     if (rank == 0)
     *         fprintf(stdout, "max buf size across all dims and procs = %d\n",
     *                 max_buf_sz_reduced);
     */

    return;
}

void init_stfw_ExpRed(Comm *comm, vpt *gm) {
    int STFW_INIT_STRATEGY = DRYRUN_SPARSE_P2P;
    switch (STFW_INIT_STRATEGY) {
        case DRYRUN_HEAVY_P2P:
            //init_stfw_nh_heavy(comm, gm);
            break;
        case DRYRUN_SPARSE_P2P:
            init_stfw_nh_nbx_ExpRed(comm, gm);
    }
}



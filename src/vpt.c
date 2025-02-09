/**
 * @author      : Nabil Abubaker (nabil.abubaker@bilkent.edu.tr)
 * @file        : vpt
 * @created     : Monday Sep 05, 2022 12:24:22 +03
 */

#include "vpt.h"
#include "mpi.h"
#include "stfw.h"
#include <stdio.h>
#include <stdlib.h>

void vpt_check_map(int *map, int *helper, int n){
    int i; for (i = 0; i < n; ++i) helper[i] = -1;
    for (i = 0; i < n; ++i) helper[map[i]] += 1; for (i = 0; i < n; ++i) 
        if(helper[i] != 0 ){ fprintf(stderr, "SpComm Error: the provided vpt-aware mapping is incorrect\n");
            exit(1); }
}

void check_vpt_dims(int *dims, int n){
    int i, nprocs, mult=1;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    for (i = 0; i < n; ++i) {
        if(dims[i] <= 0){fprintf(stderr, "ERROR: invalid Dimension Vector: dim %d is %d!\n", i, dims[i]); exit(EXIT_FAILURE);}
        mult *= dims[i]; 
    }
    if(mult != nprocs) {fprintf(stderr, "ERROR: invalid Dimension Vector: mult of dims (%d) not equal to nprocs(%d)!\n", mult, nprocs); exit(EXIT_FAILURE);}
}

/*
 * DESCRIPTION
 *	Initializes all necessary communication structures by setting their
 *	sizes and allocating them. These include:
 *  comm->sno		: Number of procs that will be communicated
 *  comm->sprocs		: The ids of these procs
 *  comm-> invprocs	: For O(1) access to procs array. Nonexistent entries =
 *-1 comm->rsizes	: Recv sizes to each proc that will be communicated
 *  comm->recv		: The buffer to get the values from other procs
 *  comm->ssizes	: Send sizes to each proc that will be communicated
 *  comm->send		: The buffer to send values to other procs
 *  comm->send_map	: The ids of the rows for the send values (one to one
 *  				  correspondance between comm->send)
 */


void gm_cstats_na(vpt *gm, int *maxSendVol, int *maxRecvVol, int *totalVol , int *maxSendMsgs, int *maxRecvMsgs, int *totalMsgs, int factor)
{
    int i;
    int rank;
    int nprocs;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

#ifdef D_SPEC_STFW
    MPI_Barrier(MPI_COMM_WORLD);
    fprintf(stderr, "Hello from gm_cstats_na\n");
#endif
    /* these include header */
    int *vstot = malloc(sizeof(*vstot) * gm->nd);
    int *vrtot = malloc(sizeof(*vrtot) * gm->nd);
    int *vsmax = malloc(sizeof(*vsmax) * gm->nd);
    int *vrmax = malloc(sizeof(*vrmax) * gm->nd);
    int *mstot = malloc(sizeof(*mstot) * gm->nd);
    int *mrtot = malloc(sizeof(*mrtot) * gm->nd);
    int *msmax = malloc(sizeof(*msmax) * gm->nd);
    int *mrmax = malloc(sizeof(*mrmax) * gm->nd);

    int vsend, vrecv, msend, mrecv;
    for (i = 0; i < gm->nd; ++i) {
        int curd = gm->order[i];
        stfw_comm *cm = &gm->sc[curd];

        vsend = cm->fw_locs[cm->fw_np]/factor;
        vrecv = cm->st_locs[cm->st_np]/factor;
        msend = cm->fw_np;
        mrecv = cm->st_np;

        MPI_Reduce(&vsend, &vstot[curd], 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&vsend, &vsmax[curd], 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&vrecv, &vrtot[curd], 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&vrecv, &vrmax[curd], 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

        MPI_Reduce(&msend, &mstot[curd], 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&msend, &msmax[curd], 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&mrecv, &mrtot[curd], 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&mrecv, &mrmax[curd], 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

    }
#ifdef D_SPEC_STFW
    MPI_Barrier(MPI_COMM_WORLD);
    fprintf(stderr, "gm_cstats_na: phase 1 done\n");
    fflush(sdbgfp);
#endif

    *totalVol = 0;  *maxSendVol= 0; *maxRecvVol = 0;
    *totalMsgs = 0; *maxSendMsgs = 0; *maxRecvMsgs= 0;

    for (i = 0; i < gm->nd; ++i) {
        *totalVol += vstot[i];
        *maxSendVol += vsmax[i];
        *maxRecvVol += vrmax[i];
        *totalMsgs += mstot[i];
        *maxSendMsgs += msmax[i];
        *maxRecvMsgs += mrmax[i];
    }

#ifdef D_SPEC_STFW
    MPI_Barrier(MPI_COMM_WORLD);
    fprintf(stderr, "gm_cstats_na: phase 2 done\n");
    fflush(sdbgfp);
#endif
    int use_maxofsum = 1;


    if (use_maxofsum) {
        vsend = vrecv = msend = mrecv = 0;
        for (i = 0; i < gm->nd; ++i) {
            stfw_comm *cm = &gm->sc[i];
            vsend += cm->fw_locs[cm->fw_np]/factor;
            vrecv += cm->st_locs[cm->st_np]/factor;
            msend += cm->fw_np;
            mrecv += cm->st_np;
        }

        int vsend_max, vrecv_max, msend_max, mrecv_max;
        MPI_Reduce(&vsend, &vsend_max, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&vrecv, &vrecv_max, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&msend, &msend_max, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&mrecv, &mrecv_max, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

        *maxSendVol = vsend_max;
        *maxRecvVol = vrecv_max; 
        *maxSendMsgs = msend_max;
        *maxRecvMsgs = mrecv_max;
    }
#ifdef D_SPEC_STFW
    MPI_Barrier(MPI_COMM_WORLD);
    fprintf(stderr, "gm_cstats_na: phase 3 done\n");
#endif

    free(vstot);
    free(vrtot);
    free(vsmax);
    free(vrmax);
    free(mstot);
    free(mrtot);
    free(msmax);
    free(mrmax);

}

void gm_cstats(vpt *gm) {
    int i;
    int rank;
    int nprocs;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    FILE *fp = NULL;
    if (rank == 0) {
        char *statfile = calloc(1024, sizeof(*statfile));
        //sprintf(statfile, "%s.stfw.cstats.%d.d%d.", inp_hg_name, nprocs, gm->nd);
        for (i = 0; i < gm->nd; ++i)
            sprintf(&statfile[strlen(statfile)], "%dx", gm->size[i]);
        statfile[strlen(statfile) - 1] = '\0';

        fp = fopen(statfile, "w");
        free(statfile);

        fprintf(fp, "dimension order: ");
        for (i = 0; i < gm->nd; ++i)
            fprintf(fp, "%d ", gm->order[i]);
        fprintf(fp, "\n");

        fprintf(fp, "max method: max of each dim, then sum\n");
    }

    /* these include header */
    int *vstot = malloc(sizeof(*vstot) * gm->nd);
    int *vrtot = malloc(sizeof(*vrtot) * gm->nd);
    int *vsmax = malloc(sizeof(*vsmax) * gm->nd);
    int *vrmax = malloc(sizeof(*vrmax) * gm->nd);
    int *mstot = malloc(sizeof(*mstot) * gm->nd);
    int *mrtot = malloc(sizeof(*mrtot) * gm->nd);
    int *msmax = malloc(sizeof(*msmax) * gm->nd);
    int *mrmax = malloc(sizeof(*mrmax) * gm->nd);

    int vsend, vrecv, msend, mrecv;
    for (i = 0; i < gm->nd; ++i) {
        int curd = gm->order[i];
        stfw_comm *cm = &gm->sc[curd];

        vsend = cm->fw_locs[cm->fw_np];
        vrecv = cm->st_locs[cm->st_np];
        msend = cm->fw_np;
        mrecv = cm->st_np;

        MPI_Reduce(&vsend, &vstot[curd], 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&vsend, &vsmax[curd], 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&vrecv, &vrtot[curd], 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&vrecv, &vrmax[curd], 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

        MPI_Reduce(&msend, &mstot[curd], 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&msend, &msmax[curd], 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&mrecv, &mrtot[curd], 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&mrecv, &mrmax[curd], 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

        if (rank == 0) {
            fprintf(fp, ">>> comm stats at dim %d size %d <<<\n", curd,
                    gm->size[curd]);
            fprintf(fp, "vstot %6d\n", vstot[curd]);
            fprintf(fp, "vsmax %6d\n", vsmax[curd]);
            fprintf(fp, "vrtot %6d\n", vrtot[curd]);
            fprintf(fp, "vrmax %6d\n", vrmax[curd]);
            fprintf(fp, "mstot %6d\n", mstot[curd]);
            fprintf(fp, "msmax %6d\n", msmax[curd]);
            fprintf(fp, "mrtot %6d\n", mrtot[curd]);
            fprintf(fp, "mrmax %6d\n", mrmax[curd]);
            fprintf(fp, "\n");
            fflush(fp);
        }
    }

    if (rank == 0) {
        fprintf(fp, ">>> total (all dimensions) <<<\n");
        int vstot_all = 0, vsmax_all = 0, vrtot_all = 0, vrmax_all = 0;
        int mstot_all = 0, msmax_all = 0, mrtot_all = 0, mrmax_all = 0;

        for (i = 0; i < gm->nd; ++i) {
            vstot_all += vstot[i];
            vsmax_all += vsmax[i];
            vrtot_all += vrtot[i];
            vrmax_all += vrmax[i];
            mstot_all += mstot[i];
            msmax_all += msmax[i];
            mrtot_all += mrtot[i];
            mrmax_all += mrmax[i];
        }

        fprintf(fp, "vstot %6d\n", vstot_all);
        fprintf(fp, "vsmax %6d\n", vsmax_all);
        fprintf(fp, "vrtot %6d\n", vrtot_all);
        fprintf(fp, "vrmax %6d\n", vrmax_all);
        fprintf(fp, "mstot %6d\n", mstot_all);
        fprintf(fp, "msmax %6d\n", msmax_all);
        fprintf(fp, "mrtot %6d\n", mrtot_all);
        fprintf(fp, "mrmax %6d\n", mrmax_all);
        fflush(fp);
        /* fclose(fp); */
    }

    if (rank == 0)
        fprintf(fp, "\n\nmax method: sum each dim, then max\n");

    vsend = vrecv = msend = mrecv = 0;
    for (i = 0; i < gm->nd; ++i) {
        stfw_comm *cm = &gm->sc[i];
        vsend += cm->fw_locs[cm->fw_np];
        vrecv += cm->st_locs[cm->st_np];
        msend += cm->fw_np;
        mrecv += cm->st_np;
    }

    int vsend_max, vrecv_max, msend_max, mrecv_max;
    MPI_Reduce(&vsend, &vsend_max, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&vrecv, &vrecv_max, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&msend, &msend_max, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&mrecv, &mrecv_max, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        fprintf(fp, ">>> all dimensions max <<<\n");
        fprintf(fp, "vsmax %6d\n", vsend_max);
        fprintf(fp, "vrmax %6d\n", vrecv_max);
        fprintf(fp, "msmax %6d\n", msend_max);
        fprintf(fp, "mrmax %6d\n", mrecv_max);
        fflush(fp);
        fclose(fp);
    }

    free(vstot);
    free(vrtot);
    free(vsmax);
    free(vrmax);
    free(mstot);
    free(mrtot);
    free(msmax);
    free(mrmax);
}

void assign_crd(vpt *gm) {
    int i, tmp, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    for (i = 0, tmp = 1; i < gm->nd; tmp *= gm->size[i], ++i)
        gm->crd[i] = (rank % (tmp * gm->size[i])) / tmp;
    return;
}

void print_crd(vpt *gm, int *crd) {
    int i;
    for (i = 0; i < gm->nd; ++i) {
        fprintf(stderr, "%d", crd[i]);
        if (i != gm->nd - 1)
            fprintf(stderr, ", ");
    }
    fprintf(stderr, "\n");

    return;
}

/* Returns the caller's proc xth neighbor on dimension d */
int get_nghbr(vpt *gm, int dim, /* on dimension d */
        int x               /* xth order in d */
        ) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    rank = gm->map[rank];
    return (x << gm->cnbit[dim]) | ((~gm->mask[dim]) & rank);
}

void print_vpt(vpt *gm) {
    int i;

    fprintf(stderr, "# dims = %d\n", gm->nd);
    fprintf(stderr, "dim sizes = ");
    for (i = 0; i < gm->nd; ++i) {
        fprintf(stderr, "%d ", gm->size[i]);
        if (i != gm->nd - 1)
            fprintf(stderr, "x ");
    }
    fprintf(stderr, "\n");

    fprintf(stderr, "cumulative sizes ");
    for (i = 0; i < gm->nd + 1; ++i)
        fprintf(stderr, "%d ", gm->csize[i]);
    fprintf(stderr, "\n");

    fprintf(stderr, "proc coords = ");
    print_crd(gm, gm->crd);

    fprintf(stderr, "dim order = ");
    for (i = 0; i < gm->nd; ++i)
        fprintf(stderr, "%d ", gm->order[i]);
    fprintf(stderr, "\n");

    fprintf(stderr, "number of bits = ");
    for (i = 0; i < gm->nd; ++i)
        fprintf(stderr, "%d ", gm->nbit[i]);
    fprintf(stderr, "\n");

    fprintf(stderr, "cumulative number of bits = ");
    for (i = 0; i < gm->nd + 1; ++i)
        fprintf(stderr, "%d ", gm->cnbit[i]);
    fprintf(stderr, "\n");

    fprintf(stderr, "masks\n");
    for (i = 0; i < gm->nd; ++i)
        fprintf(stderr, "%x\n", gm->mask[i]);

    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    fprintf(stderr, "comm dims:\n");
    for (i = 0; i < nprocs; ++i)
        fprintf(stderr, "  dim of P%d = %d\n", i, gm->comm_dim[i]);

    fprintf(stderr, "\n");

    return;
}

void print_hdr(vpt *gm, int **hdr_send, int **hdr_recv) {
    int rank, nprocs;
    int i, d;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    fprintf(stderr, "\nsend headers of proc %d\n", rank);
    for (i = 0; i < gm->nd; ++i) {
        d = gm->order[i];
        int crd = gm->crd[d];
        int hlen = 3 * (nprocs / gm->size[d]) + 1;

        fprintf(stderr, "communication no %d, dimension %d, my crd = %d\n", i, d,
                crd);
        fprintf(stderr, "  number of max packets in this dimension = %d\n",
                nprocs / gm->size[d]);
        fprintf(stderr, "  number of max procs to communicate = %d\n",
                gm->size[d] - 1);
        fprintf(stderr, "  hlen = %d\n", hlen);

        int j;
        for (j = 0; j < gm->size[d]; ++j) {
            if (j == crd) /* skip proc itself */
                continue;

            int *ptr = &hdr_send[d][hlen * j];
            int npackets = *ptr;

            fprintf(stderr, "  proc crd %d on dim %d, P%d, npackets = %d\n", j, d,
                    get_nghbr(gm, d, j), npackets);

            int k;
            ptr += 1;
            for (k = 0; k < npackets; ++k, ptr += 3)
                fprintf(stderr, "    packet %d, content (P%d, P%d, %d)\n", k + 1, *ptr,
                        *(ptr + 1), *(ptr + 2));
        }
    }

    return;
}

/* @TODO Alter this to get coordinates after some communication order */
inline void get_all_crd(vpt *gm, int proc_rank, int *crd) {
    int i, tmp;
    /* for (i = 0, tmp = 1; i < gm->nd; tmp *= gm->size[i], ++i) */
    /* 	crd[i] = (proc_rank % (tmp * gm->size[i])) / tmp; */
    /* return; */

    for (i = 0, tmp = 0; i < gm->nd; tmp += gm->nbit[i++])
        crd[i] = (proc_rank & gm->mask[i]) >> tmp;

    return;
}

int get_crd(vpt *gm, int proc_rank, int d) {
    return (proc_rank & gm->mask[d]) >> gm->cnbit[d];
}

/* NA: this is for re-calculating comm_dim   */
void init_comm_dim(vpt *gm){
    int i, nprocs, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    rank = gm->map[rank];

    int *tmp_crds = malloc(sizeof(*tmp_crds) * gm->nd);
    for (i = 0; i < nprocs; ++i) {

        /* Determine proc comm dimensions */
        if (i == rank) {
            gm->comm_dim[i] = NE;
            continue;
        }

        get_all_crd(gm, i, tmp_crds);
        int d = 0;
        while (tmp_crds[gm->order[d]] == gm->crd[gm->order[d]])
            ++d;
        assert(d >= 0 && d < gm->nd); /* at least one dim should differ */
        gm->comm_dim[i] = gm->order[d];
    }

    free(tmp_crds);

}

void init_vpt(vpt *gm, dim_order ord, int *map) {
    /* va_list valist; */
    int i, j, tmp;
    int nprocs, rank;
    int nd = gm->nd;
    
    check_vpt_dims(gm->size, gm->nd);

    /* gm->nd	  = nd; */
    gm->order = malloc(sizeof(*gm->order) * nd);
    /* gm->size  = malloc(sizeof(*gm->size)*nd); */
    gm->csize = calloc(nd + 1, sizeof(*gm->csize));
    gm->nbit = calloc(nd, sizeof(*gm->nbit));
    gm->cnbit = calloc(nd + 1, sizeof(*gm->cnbit));
    gm->mask = malloc(sizeof(*gm->mask) * nd);
    gm->crd = malloc(sizeof(*gm->crd) * nd);
    gm->sc = malloc(sizeof(*gm->sc) * nd);
    /* va_start(valist, nd); */
    for (i = 0; i < nd; ++i) {
        /* gm->size[i]	   = va_arg(valist, int); */
        gm->csize[i + 1] = gm->csize[i] + gm->size[i];
    }

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    gm->map = malloc(sizeof(*gm->map)*nprocs);
    gm->inv_map = malloc(sizeof(*gm->inv_map)*nprocs);
    
    if(map != NULL){
       vpt_check_map(map, gm->map, nprocs);
       for (i = 0; i < nprocs; ++i){gm->map[i]=map[i]; gm->inv_map[map[i]]=i;}
    }    
    else
        for (i = 0; i < nprocs; ++i){gm->map[i]=i; gm->inv_map[i]=i;}

    tmp = 1;
    rank = gm->map[rank];
    int tmp2 = 0;
    for (i = 0; i < nd; ++i) {
        if(ord == D_LINEARREV) gm->order[i] = nd-i-1; 
        else  gm->order[i] = i; /* LINEAR here, will be changed if RANDOM */

        /* all dim sizes must be powers of 2 */
        assert((gm->size[i] & (gm->size[i] - 1)) == 0);
        tmp *= gm->size[i];
        j = gm->size[i];
        while (j >>= 1)
            ++(gm->nbit[i]);
        gm->cnbit[i + 1] = gm->cnbit[i] + gm->nbit[i];

        gm->mask[i] = (gm->size[i] - 1) << tmp2;
        gm->crd[i] = (rank & gm->mask[i]) >> tmp2;
        tmp2 += gm->nbit[i];
        init_cm_p(&(gm->sc[i]));      
        gm->sc[i].fw_np = 0;
        gm->sc[i].st_np = 0;
    }

    assert(tmp == nprocs);

    if (ord == D_RANDOM) {
        srand(17); /* must be same for all processes, don't use
                    * time here */
        for (i = 0; i < nd - 1; ++i) {
            int r = rand() % (nd - i);
            tmp = gm->order[r];
            gm->order[r] = gm->order[nd - i - 1];
            gm->order[nd - i - 1] = tmp;
        }
    }

    /* gm->order[0] = 1; */
    /* gm->order[1] = 0; */
    /* gm->order[2] = 2; */

    gm->comm_dim = malloc(sizeof(*gm->comm_dim) * nprocs);

    init_comm_dim(gm);
    MPI_Barrier(MPI_COMM_WORLD);

#ifdef D_SPEC_GMESH
    print_vpt(gm);
    /* for (i = 0; i < gm->nd; ++i) */
    /* 	fprintf(stderr, "my coord at %d: %d\n", i, get_crd(gm, rank, i)); */
#endif

    return;
}

void free_vpt(vpt *gm) {

#ifdef D_SPEC_STFW
    fprintf(stderr, "> inside free vpt\n\n"); 
    fflush(sdbgfp);
#endif
    int i;
    free(gm->order); free(gm->size); free(gm->csize); free(gm->nbit); free(gm->cnbit); free(gm->mask); free(gm->crd); free(gm->map); free(gm->inv_map);

#ifdef D_SPEC_STFW
    fprintf(stderr, "> free vpt arrays done\n"); 
    fflush(sdbgfp);
#endif
    for (i = 0; i < gm->nd; ++i) {
        if (gm->sc[i].fw_buf) free(gm->sc[i].fw_buf);
        if (gm->sc[i].fw_procs) free(gm->sc[i].fw_procs);
        if (gm->sc[i].fw_offset) free(gm->sc[i].fw_offset);
        if(gm->sc[i].fw_locs) free(gm->sc[i].fw_locs);
        if(gm->sc[i].fw_imap) free(gm->sc[i].fw_imap);

#ifdef D_SPEC_STFW
        fprintf(stderr, "> freeing fw buffers of dim=%d done\n", i); 
        fflush(sdbgfp);
#endif
        if (gm->sc[i].st_buf) free(gm->sc[i].st_buf);
        if (gm->sc[i].st_procs) free(gm->sc[i].st_procs);
        if (gm->sc[i].st_reqs) free(gm->sc[i].st_reqs);
        if (gm->sc[i].st_stts) free(gm->sc[i].st_stts);
        if(gm->sc[i].st_locs) free(gm->sc[i].st_locs);
#ifdef D_SPEC_STFW
        fprintf(stderr, "> freeing st buffers of dim=%d done\n", i); 
        fflush(sdbgfp);
#endif
        if (_DE) {
            if(gm->sc[i].fw_imap_d) free(gm->sc[i].fw_imap_d);
            if(gm->sc[i].fw_offset_d) free(gm->sc[i].fw_offset_d);
            if(gm->sc[i].st_reqs_d) free(gm->sc[i].st_reqs_d);
            if(gm->sc[i].st_stts_d) free(gm->sc[i].st_stts_d);
#ifdef D_SPEC_STFW
            fprintf(stderr, "> freeing additional dual buffers of dim=%d done\n", i); 
            fflush(sdbgfp);
#endif
        }
        if(gm->sc[i].pck_cnt) free(gm->sc[i].pck_cnt);
        if (gm->sc[i].pck_sizes) free(gm->sc[i].pck_sizes);
        if (gm->sc[i].pck_cp) free(gm->sc[i].pck_cp);
#ifdef D_SPEC_STFW
        fprintf(stderr, "> freeing pck buffers of dim=%d done\n", i); 
        fflush(sdbgfp);
#endif
#ifdef D_SPEC_STFW
        if (gm->sc[i].pck_fw_idx) free(gm->sc[i].pck_fw_idx);
        if (gm->sc[i].pck_src) free(gm->sc[i].pck_src);
        if (gm->sc[i].pck_dest) free(gm->sc[i].pck_dest);

        fprintf(stderr, "> freeing dbg buffers of dim=%d done\n", i); 
        fflush(sdbgfp);
#endif

        if(gm->sc[i].unit_RepCnt) free(gm->sc[i].unit_RepCnt);
        if(gm->sc[i].unit_ptr) free(gm->sc[i].unit_ptr);
    }
    free(gm->sc);

    free(gm->comm_dim);

#ifdef D_SPEC_STFW
    fprintf(stderr, "> all free vpt done\n"); 
#endif

}


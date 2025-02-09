//include "dbg_print.h"
#include "vptcomm.h"
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

void arrange_comm_nh(vpt *gm, int dim, int *fw_sz, int fw_tot_sz,
        int *st_sz, int st_tot_sz, int *pck_cnt) {
    int i;
    stfw_comm *cm = &gm->sc[dim];

    cm->fw_buf = NULL;
    cm->fw_procs = NULL;
    cm->fw_locs = NULL;
    cm->fw_offset = NULL;
    assert(cm->fw_np >= 0);
    if (fw_tot_sz)
        cm->fw_buf = malloc(sizeof(*cm->fw_buf) * fw_tot_sz);
    if (cm->fw_np) {
        cm->fw_procs = malloc(sizeof(*cm->fw_procs) * cm->fw_np);
        cm->fw_offset = malloc(sizeof(*cm->fw_offset) * cm->fw_np);
    }
    cm->fw_locs = calloc(cm->fw_np + 1, sizeof(*cm->fw_locs));
    cm->fw_imap = malloc(sizeof(*cm->fw_imap) * gm->size[dim]);
    cm->fw_locs[0] = 0;

    cm->st_buf = NULL;
    cm->st_procs = NULL;
    cm->st_reqs = NULL;
    cm->st_stts = NULL;
    assert(cm->st_np >= 0);
    if (st_tot_sz)
        cm->st_buf = malloc(sizeof(*cm->st_buf) * st_tot_sz);
    if (cm->st_np) {
        cm->st_procs = malloc(sizeof(*cm->st_procs) * cm->st_np);
        cm->st_reqs = malloc(sizeof(*cm->st_reqs) * cm->st_np);
        cm->st_stts = malloc(sizeof(*cm->st_stts) * cm->st_np);
    }
    cm->st_locs = calloc((cm->st_np + 1), sizeof(*cm->st_locs));
    cm->st_locs[0] = 0;

    cm->pck_cnt = calloc(cm->st_np + 1, sizeof(*cm->pck_cnt));
    cm->pck_cnt[0] = 0;

    int fw_idx, st_idx;
    for (i = 0, fw_idx = 0, st_idx = 0; i < gm->size[dim]; ++i) {
        cm->fw_imap[i] = NE;

        if (gm->crd[dim] == i)
            continue;

        if (fw_sz[i]) {
            cm->fw_offset[fw_idx] = cm->fw_locs[fw_idx];
            //cm->fw_buf[cm->fw_locs[fw_idx]] = 0.0;
            cm->fw_locs[fw_idx + 1] = cm->fw_locs[fw_idx] + fw_sz[i];
            cm->fw_imap[i] = fw_idx;
            cm->fw_procs[fw_idx++] = get_nghbr(gm, dim, i);
        }

        if (st_sz[i] > 0) {
            cm->st_locs[st_idx + 1] = cm->st_locs[st_idx] + st_sz[i];
            assert(pck_cnt[i] != 0);
            cm->pck_cnt[st_idx + 1] = cm->pck_cnt[st_idx] + pck_cnt[i];
            cm->st_procs[st_idx++] = get_nghbr(gm, dim, i);
        }
    }

    //fprintf(stderr, "st_np=%d st_idx=%d\n", cm->st_np, st_idx);
    assert(cm->fw_np == fw_idx);
    assert(cm->st_np == st_idx);

#ifdef D_SPEC_STFW
    fprintf(sdbgfp, "\n\nfinalizing communicator at dim %d\n", dim);
    fprintf(sdbgfp, "forward:\n");
    fprintf(sdbgfp, "  total buffer size %d\n", fw_tot_sz);
    fprintf(sdbgfp, "  number of procs %d\n", cm->fw_np);
    for (i = 0; i < cm->fw_np; ++i)
        fprintf(sdbgfp,
                "    P%d, send amount = %d, locations = [%d, %d), "
                "offset = %d\n",
                cm->fw_procs[i], cm->fw_locs[i + 1] - cm->fw_locs[i],
                cm->fw_locs[i], cm->fw_locs[i + 1], cm->fw_offset[i]);
    fprintf(sdbgfp, "  inverse fwd map:\n");
    for (i = 0; i < gm->size[dim]; ++i)
        fprintf(sdbgfp, "    crd %d: %d\n", i, cm->fw_imap[i]);

    fprintf(sdbgfp, "store:\n");
    fprintf(sdbgfp, "  total buffer size %d\n", st_tot_sz);
    fprintf(sdbgfp, "  number of procs %d\n", cm->st_np);
    for (i = 0; i < cm->st_np; ++i)
        fprintf(sdbgfp, "    P%d, recv amount = %d, locations = [%d, %d)\n",
                cm->st_procs[i], cm->st_locs[i + 1] - cm->st_locs[i],
                cm->st_locs[i], cm->st_locs[i + 1]);

    fprintf(sdbgfp, "packet:\n");
    fprintf(sdbgfp, "  number of procs %d\n", cm->st_np);
    for (i = 0; i < cm->st_np; ++i)
        fprintf(sdbgfp, "    P%d, pck count %d, locations = [%d, %d)\n",
                cm->st_procs[i], cm->pck_cnt[i + 1] - cm->pck_cnt[i],
                cm->pck_cnt[i], cm->pck_cnt[i + 1]);

    fflush(sdbgfp);
#endif

    return;
}


/* nh: no headers version */
void init_stfw_nh_heavy(Comm *comm, vpt *gm) {
#define TAG_HDR 7

    int i, j, k;
    int rank, nprocs;
    int tmp;

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int *tmp_crds = malloc(sizeof(*tmp_crds) * gm->nd);

    /* data sizes communicated with each proc in each dimension  */
    int **fw_sz = malloc(sizeof(*fw_sz) * gm->nd);
    int **st_sz = malloc(sizeof(*st_sz) * gm->nd);
    int *fw_tot_sz = calloc(gm->nd, sizeof(*fw_tot_sz));
    int *st_tot_sz = calloc(gm->nd, sizeof(*st_tot_sz));

    int max_buf_sz = -1;

    /* headers */
    int **hdr_send = malloc(sizeof(*hdr_send) * gm->nd);
    int **hdr_recv = malloc(sizeof(*hdr_recv) * gm->nd);
    MPI_Request **hdr_reqs = malloc(sizeof(*hdr_reqs) * gm->nd);

    /* packets helpers */
    int **pck_cnt = malloc(sizeof(*pck_cnt) * gm->nd);
    struct pck_cp_data {
        /* if d = -1, packet is for me and dim_idx is the proc I recved data */
        /* if d != -1, use fw_idx, if it is -1 use ptr */
        int d;         /* which dim the packet will be forwarded in */
        int dim_idx;   /* index in the respective dim */
        int fw_idx;    /* forward index for fw_buf */
        vptcomm_real_t *ptr; /* direct address for vector p */
    };
    struct pck_cp_data **pck_cp = malloc(sizeof(*pck_cp) * gm->nd);

    for (i = 0; i < gm->nd; ++i) {
        int dim_size = gm->size[i];
        /* each packet: src dest size */
        int hlen = 3 * (nprocs / dim_size) + 1;
        hdr_send[i] = calloc(dim_size * hlen, sizeof(**hdr_send));
        hdr_recv[i] = calloc(dim_size * hlen, sizeof(**hdr_recv));
        hdr_reqs[i] = malloc(sizeof(**hdr_reqs) * dim_size);

        fw_sz[i] = calloc(dim_size, sizeof(**fw_sz));
        st_sz[i] = calloc(dim_size, sizeof(**st_sz));

        pck_cnt[i] = calloc(dim_size, sizeof(**pck_cnt));
    }

#ifdef D_SPEC_STFW
    fprintf(sdbgfp, "> init_stfw_gmesh\n");
    fprintf(sdbgfp, "number of processors to communicate = %d\n", comm->sno);
    for (i = 0; i < comm->sno; ++i)
        fprintf(sdbgfp, "  P %d, send size = %d, recv size = %d\n", gm->map[comm->sprocs[i]],
                comm->ssizes[i], comm->rsizes[i]);
    fprintf(sdbgfp, "\n");
    fflush(sdbgfp);
#endif

    /* fill inv_recv_map */
    int *inv_recv_map = malloc(sizeof(*inv_recv_map) * nprocs);
    for (i = 0; i < comm->rno; ++i) {
        inv_recv_map[gm->map[comm->rprocs[i]]] = gm->map[i];
    }
    /* process my send list and fill in the communicators in the respective
     * dimension */
    for (i = 0; i < comm->sno; ++i) {
        int proc = gm->map[comm->sprocs[i]];
        int d = gm->comm_dim[proc];

        /* NA: what the hell is this?? */
        //gm->inv_recv_map[comm->sprocs[i]] = i;

#ifdef D_SPEC_STFW
        fprintf(sdbgfp, "will send data to P%d in comm dim %d, ", proc, d);
        fflush(sdbgfp);
#endif

        /* add the header */
        int hlen = 3 * (nprocs / gm->size[d]) + 1;
        int proc_crd = get_crd(gm, proc, d);
        int *ptr = &hdr_send[d][hlen * proc_crd];
        int packet_id = *ptr;
        ptr += 3 * packet_id + 1;
        *ptr = gm->map[rank];                  /* src */
        *(ptr + 1) = proc;            /* dest */
        *(ptr + 2) = comm->ssizes[i]; /* size */
        hdr_send[d][hlen * proc_crd] += 1;
        assert(packet_id >= 0 && packet_id < (nprocs / gm->size[d]));

#ifdef D_SPEC_STFW
        fprintf(sdbgfp,
                "  adding to respective header, this is packet %d, "
                "size %d, proc crd in respective dim is %d, max packets %d\n",
                hdr_send[d][hlen * proc_crd], comm->ssizes[i], proc_crd,
                nprocs / gm->size[d]);
        fflush(sdbgfp);
#endif

        /* update forward communicators */
        if (fw_sz[d][proc_crd] == 0) /* first packet */
            ++(gm->sc[d].fw_np);
        fw_sz[d][proc_crd] += comm->ssizes[i];
        fw_tot_sz[d] += comm->ssizes[i];

#ifdef D_SPEC_STFW
        fprintf(sdbgfp,
                "  [comm fw] fw_sz and fw_np: fw_sz[%d][%d] = %d, "
                "gm->sc[%d].fw_np = %d, fw_tot_sz[%d] = %d\n",
                d, proc_crd, fw_sz[d][proc_crd], d, gm->sc[d].fw_np, d,
                fw_tot_sz[d]);
#endif
    }

#ifdef D_SPEC_STFW
    print_hdr(gm, hdr_send, hdr_recv);
#endif

    for (i = 0; i < gm->nd; ++i) {
        int curd = gm->order[i];
        int my_crd = gm->crd[curd];
        int cur_hlen = 3 * (nprocs / gm->size[curd]) + 1;
        int *cur_hsend = hdr_send[curd];
        int *cur_hrecv = hdr_recv[curd];
        MPI_Request *cur_hreqs = hdr_reqs[curd];

#ifdef D_SPEC_STFW
        fprintf(sdbgfp, "\n\n");
        fprintf(sdbgfp,
                "comm operations turn %d, dimension %d, my_crd = %d, "
                "hlen = %d, max #packets = %d\n",
                i, curd, my_crd, cur_hlen, nprocs / gm->size[curd]);
        fflush(sdbgfp);
#endif

        for (j = 0; j < gm->size[curd]; ++j) {
            if (j != my_crd)
                MPI_Irecv(&cur_hrecv[j * cur_hlen], cur_hlen, MPI_INT,
                        gm->inv_map[get_nghbr(gm, curd, j)], TAG_HDR, MPI_COMM_WORLD,
                        &cur_hreqs[j]);
        }

        /* process hdr_recv of the prev step to form send data of this step */
        if (i != 0) {
            int prev_dim = gm->order[i - 1];
            int prev_crd = gm->crd[prev_dim];
            int prev_hlen = 3 * (nprocs / gm->size[prev_dim]) + 1;
            int *prev_hrecv = hdr_recv[prev_dim];
            MPI_Request *prev_hreqs = hdr_reqs[prev_dim];

#ifdef D_SPEC_STFW
            fprintf(sdbgfp, "\nlisting recv data of prev turn...\n");
            fprintf(sdbgfp, "  prev dim = %d, my prev crd = %d\n", prev_dim, prev_crd);
            fflush(sdbgfp);
#endif

            /* wait for recvs and update packet-related data */
            int tot_pck_cnt = 0;
            for (j = 0; j < gm->size[prev_dim]; ++j) {
                if (j == prev_crd)
                    continue;

                MPI_Status stts;
                MPI_Wait(&prev_hreqs[j], &stts);

                int npackets = prev_hrecv[j * prev_hlen];
                assert(npackets >= 0 && npackets <= nprocs / gm->size[prev_dim]);

                /* update store communicators */
                if (npackets > 0 && st_sz[prev_dim][j] == 0) /* first packet */
                    ++(gm->sc[prev_dim].st_np);

                tot_pck_cnt += npackets;
                pck_cnt[prev_dim][j] = npackets;
            }

#ifdef D_SPEC_STFW
            fprintf(sdbgfp, "  total recved packet count %d\n", tot_pck_cnt);
            fflush(sdbgfp);
#endif

            stfw_comm *cm = &gm->sc[prev_dim];
            cm->pck_sizes = NULL;
            pck_cp[prev_dim] = NULL;
#ifdef D_SPEC_STFW
            cm->pck_dest = NULL;
            cm->pck_src = NULL;
#endif
            if (tot_pck_cnt) {
                cm->pck_sizes = malloc(sizeof(*cm->pck_sizes) * tot_pck_cnt);
                pck_cp[prev_dim] = malloc(sizeof(**pck_cp) * tot_pck_cnt);
#ifdef D_SPEC_STFW
                cm->pck_dest = malloc(sizeof(*cm->pck_dest) * tot_pck_cnt);
                cm->pck_src = malloc(sizeof(*cm->pck_src) * tot_pck_cnt);
#endif
            }

            /* process packets */
            int pck_idx = 0;
            for (j = 0; j < gm->size[prev_dim]; ++j) {
                if (j == prev_crd)
                    continue;

                int npackets = prev_hrecv[j * prev_hlen];
                int *ptr = NULL;
                if (npackets > 0) /* non-empty recv */
                    ptr = &prev_hrecv[j * prev_hlen + 1];

#ifdef D_SPEC_STFW
                fprintf(sdbgfp,
                        "  recv from P%d (crd %d at dim %d), "
                        "npackets = %d, contents:\n",
                        get_nghbr(gm, prev_dim, j), j, prev_dim, npackets);
                fflush(sdbgfp);
#endif

                for (k = 0; k < npackets; ++k, ptr += 3) {
                    int pck_src = *ptr;
                    int pck_dest = *(ptr + 1);
                    int pck_size = *(ptr + 2);
                    cm->pck_sizes[pck_idx] = pck_size;
#ifdef D_SPEC_STFW
                    cm->pck_dest[pck_idx] = pck_dest;
                    cm->pck_src[pck_idx] = pck_src;
#endif

#ifdef D_SPEC_STFW
                    fprintf(sdbgfp, "    P%d P%d %d => ", pck_src, pck_dest, pck_size);
                    fflush(sdbgfp);
#endif

                    /* update store communicators */
                    st_sz[prev_dim][j] += pck_size;
                    st_tot_sz[prev_dim] += pck_size;

                    /* packet is for me */
                    if (pck_dest == gm->map[rank]) {
#ifdef D_SPEC_STFW
                        assert(inv_recv_map[pck_src] != NE);
                        fprintf(sdbgfp,
                                "package is for me, pck_idx = %d, "
                                "tuple info d = -1, dim_dix = %d, ptr = %p\n",
                                pck_idx, pck_src, comm->recv[gm->inv_map[inv_recv_map[pck_src]]]);

                        /* also assert recved proc rank is in my recv list */
                        int l;
                        int rsize = -1;
                        for (l = 0; l < comm->rno; ++l) {
                            if (gm->map[comm->rprocs[l]] == pck_src) {
                                rsize = comm->rsizes[l];
                                break;
                            }
                        }
                        fprintf(sdbgfp, "rsize=%d pck_size=%d\n", rsize, pck_size);
                        fflush(sdbgfp);
                        assert(l >= 0 && l < comm->rno && rsize == pck_size);
                        fprintf(sdbgfp,
                                "    found packet in my recv list %d = P%d, "
                                "size = %d\n",
                                l, gm->map[comm->rprocs[l]], rsize);
                        fflush(sdbgfp);
#endif
                        pck_cp[prev_dim][pck_idx].d = -1;
                        pck_cp[prev_dim][pck_idx].dim_idx = pck_src;
                        pck_cp[prev_dim][pck_idx++].ptr =
                            comm->recv[gm->inv_map[inv_recv_map[pck_src]]];
                        continue;
                    }

                    /* packet is to be forwarded */
                    int d = gm->comm_dim[pck_dest];
                    int hlen = 3 * (nprocs / gm->size[d]) + 1;
                    int proc_crd = get_crd(gm, pck_dest, d);
                    int *ptr2 = &hdr_send[d][hlen * proc_crd];
                    int packet_id = *ptr2;
                    ptr2 += 3 * packet_id + 1;
                    *ptr2 = pck_src;
                    *(ptr2 + 1) = pck_dest;
                    *(ptr2 + 2) = pck_size;
                    hdr_send[d][hlen * proc_crd] += 1;
                    assert(packet_id >= 0 && packet_id < (nprocs / gm->size[d]));

#ifdef D_SPEC_STFW
                    fprintf(sdbgfp,
                            "adding to the comm in dim %d, proc crd %d, "
                            "packet id %d, size %d, max packets %d\n",
                            d, proc_crd, hdr_send[d][hlen * proc_crd], pck_size,
                            nprocs / gm->size[d]);
                    fprintf(sdbgfp,
                            "    package will be forwarded, "
                            "pck_idx = %d, tuple info d = %d, dim_dix = %d, "
                            "fw_idx = %d\n",
                            pck_idx, d, proc_crd, fw_sz[d][proc_crd]);
                    fflush(sdbgfp);
#endif

                    pck_cp[prev_dim][pck_idx].d = d;
                    pck_cp[prev_dim][pck_idx].dim_idx = proc_crd;
                    pck_cp[prev_dim][pck_idx++].fw_idx = fw_sz[d][proc_crd];

                    /* update forward communicators */
                    if (fw_sz[d][proc_crd] == 0) /* first packet */
                        ++(gm->sc[d].fw_np);
                    fw_sz[d][proc_crd] += pck_size;
                    fw_tot_sz[d] += pck_size;

#ifdef D_SPEC_STFW
                    fprintf(sdbgfp,
                            "    [comm fw] fw_sz and fw_np: "
                            "fw_sz[%d][%d] = %d, gm->sc[%d].fw_np = %d, "
                            "fw_tot_sz[%d] = %d\n",
                            d, proc_crd, fw_sz[d][proc_crd], d, gm->sc[d].fw_np, d,
                            fw_tot_sz[d]);
#endif
                }

#ifdef D_SPEC_STFW
                fprintf(sdbgfp,
                        "    [comm st] st_sz and st_np: "
                        "st_sz[%d][%d] = %d, gm->sc[%d].st_np = %d, "
                        "st_tot_sz[%d] = %d\n",
                        prev_dim, j, st_sz[prev_dim][j], prev_dim,
                        gm->sc[prev_dim].st_np, prev_dim, st_tot_sz[prev_dim]);
                fflush(sdbgfp);
#endif
            }

            /* now can finalize the previous dimension's communicators */
            arrange_comm_nh(gm, prev_dim, fw_sz[prev_dim], fw_tot_sz[prev_dim],
                    st_sz[prev_dim], st_tot_sz[prev_dim], pck_cnt[prev_dim]);
            if (fw_tot_sz[prev_dim] + st_tot_sz[prev_dim] > max_buf_sz)
                max_buf_sz = fw_tot_sz[prev_dim] + st_tot_sz[prev_dim];

#ifdef D_SPEC_STFW
            assert(pck_idx == tot_pck_cnt);
            fprintf(sdbgfp, "recv packet info:\n");
            for (j = 0; j < cm->st_np; ++j) {
                fprintf(sdbgfp, "  P%d pck count %d\n", cm->st_procs[j],
                        cm->pck_cnt[j + 1] - cm->pck_cnt[j]);
                int l;
                for (l = cm->pck_cnt[j]; l < cm->pck_cnt[j + 1]; ++l)
                    fprintf(sdbgfp, "    packet %d size %d, src %d, dest %d\n",
                            l - cm->pck_cnt[j] + 1, cm->pck_sizes[l], cm->pck_src[l],
                            cm->pck_dest[l]);
            }
#endif
        }

#ifdef D_SPEC_STFW
        fprintf(sdbgfp, "\nissuing sends...\n");
        fflush(sdbgfp);
#endif

        for (j = 0; j < gm->size[curd]; ++j) {
            if (j != my_crd) {
#ifdef D_SPEC_STFW
                fprintf(sdbgfp,
                        "  sending to P%d (crd %d at dim %d), "
                        "npackets = %d, send size = %d\n",
                        get_nghbr(gm, curd, j), j, curd, cur_hsend[j * cur_hlen],
                        cur_hsend[j * cur_hlen] * 3 + 1);
                fflush(sdbgfp);
#endif

                MPI_Send(&cur_hsend[j * cur_hlen], cur_hsend[j * cur_hlen] * 3 + 1,
                        MPI_INT, gm->inv_map[get_nghbr(gm, curd, j)], TAG_HDR, MPI_COMM_WORLD);
            }
        }
    }

    /* wait for the final recvs, these will not be used to form new forwarding,
     * they are final and will be written over to input vector */
    int prev_dim = gm->order[gm->nd - 1];
    int prev_crd = gm->crd[prev_dim];
    int prev_hlen = 3 * (nprocs / gm->size[prev_dim]) + 1;
    int *prev_hrecv = hdr_recv[prev_dim];
    MPI_Request *prev_hreqs = hdr_reqs[prev_dim];

#ifdef D_SPEC_STFW
    fprintf(sdbgfp, "\nlisting recv data of prev turn...\n");
    fprintf(sdbgfp, "  prev dim = %d, my prev crd = %d\n", prev_dim, prev_crd);
    fflush(sdbgfp);
#endif

    /* wait for recvs and update packet-related data */
    int tot_pck_cnt = 0;
    for (j = 0; j < gm->size[prev_dim]; ++j) {
        if (j == prev_crd)
            continue;

        MPI_Status stts;
        MPI_Wait(&prev_hreqs[j], &stts);

        int npackets = prev_hrecv[j * prev_hlen];
        assert(npackets >= 0 && npackets <= nprocs / gm->size[prev_dim]);

        /* update store communicators */
        if (npackets > 0 && st_sz[prev_dim][j] == 0) /* first packet */
            ++(gm->sc[prev_dim].st_np);

        tot_pck_cnt += npackets;
        pck_cnt[prev_dim][j] = npackets;
    }

#ifdef D_SPEC_STFW
    fprintf(sdbgfp, "  total recved packet count %d\n", tot_pck_cnt);
    fflush(sdbgfp);
#endif

    stfw_comm *cm = &gm->sc[prev_dim];
    cm->pck_sizes = NULL;
    pck_cp[prev_dim] = NULL;
#ifdef D_SPEC_STFW
    cm->pck_dest = NULL;
    cm->pck_src = NULL;
#endif
    if (tot_pck_cnt) {
        cm->pck_sizes = malloc(sizeof(*cm->pck_sizes) * tot_pck_cnt);
        pck_cp[prev_dim] = malloc(sizeof(**pck_cp) * tot_pck_cnt);
#ifdef D_SPEC_STFW
        cm->pck_dest = malloc(sizeof(*cm->pck_dest) * tot_pck_cnt);
        cm->pck_src = malloc(sizeof(*cm->pck_src) * tot_pck_cnt);
#endif
    }

    /* process packets */
    int pck_idx = 0;
    for (j = 0; j < gm->size[prev_dim]; ++j) {
        if (j == prev_crd)
            continue;

        int npackets = prev_hrecv[j * prev_hlen];
        int *ptr = NULL;
        if (npackets != 0) /* non-empty recv */
            ptr = &prev_hrecv[j * prev_hlen + 1];

#ifdef D_SPEC_STFW
        fprintf(sdbgfp,
                "  recv from P%d (crd %d at dim %d), "
                "npackets = %d, contents:\n",
                get_nghbr(gm, prev_dim, j), j, prev_dim, npackets);
        fflush(sdbgfp);
#endif

        for (k = 0; k < npackets; ++k, ptr += 3) {
            int pck_src = *ptr;
            int pck_dest = *(ptr + 1);
            int pck_size = *(ptr + 2);
            cm->pck_sizes[pck_idx] = pck_size;
#ifdef D_SPEC_STFW
            cm->pck_dest[pck_idx] = pck_dest;
            cm->pck_src[pck_idx] = pck_src;
#endif

#ifdef D_SPEC_STFW
            fprintf(sdbgfp, "    P%d P%d %d => ", pck_src, pck_dest, pck_size);
            fflush(sdbgfp);
#endif

            /* this is final packet, destination should be me */
            assert(pck_dest == gm->map[rank]);

            /* also assert recved proc rank is in my recv list */
#ifdef D_SPEC_STFW
            assert(inv_recv_map[pck_src] != NE);
            fprintf(sdbgfp,
                    "package is for me, pck_idx = %d, "
                    "tuple info d = -1, dim_dix = %d, ptr = %p\n",
                    pck_idx, pck_src, comm->recv[gm->inv_map[inv_recv_map[pck_src]]]);
            int l;
            int rsize = -1;
            for (l = 0; l < comm->rno; ++l) {
                if (gm->map[comm->rprocs[l]] == pck_src) {
                    rsize = comm->rsizes[l];
                    break;
                }
            }
            assert(l >= 0 && l < comm->rno && rsize == pck_size);
            fprintf(sdbgfp,
                    "    found packet in my recv list %d = P%d, "
                    "size = %d\n",
                    l, gm->map[comm->sprocs[l]], rsize);
            fflush(sdbgfp);
#endif

            /* update store communicators */
            st_sz[prev_dim][j] += pck_size;
            st_tot_sz[prev_dim] += pck_size;

            pck_cp[prev_dim][pck_idx].d = -1;
            pck_cp[prev_dim][pck_idx].dim_idx = pck_src;
            pck_cp[prev_dim][pck_idx++].ptr = comm->recv[gm->inv_map[inv_recv_map[pck_src]]];
        }

#ifdef D_SPEC_STFW
        fprintf(sdbgfp,
                "    [comm st] st_sz and st_np: "
                "st_sz[%d][%d] = %d, gm->sc[%d].st_np = %d, "
                "st_tot_sz[%d] = %d\n",
                prev_dim, j, st_sz[prev_dim][j], prev_dim, gm->sc[prev_dim].st_np,
                prev_dim, st_tot_sz[prev_dim]);
        fflush(sdbgfp);
#endif
    }

    /* now can finalize the previous dimension's communicators */
    arrange_comm_nh(gm, prev_dim, fw_sz[prev_dim], fw_tot_sz[prev_dim],
            st_sz[prev_dim], st_tot_sz[prev_dim], pck_cnt[prev_dim]);
    if (fw_tot_sz[prev_dim] + st_tot_sz[prev_dim] > max_buf_sz)
        max_buf_sz = fw_tot_sz[prev_dim] + st_tot_sz[prev_dim];

#ifdef D_SPEC_STFW
    assert(pck_idx == tot_pck_cnt);
    fprintf(sdbgfp, "recv packet info:\n");
    for (j = 0; j < cm->st_np; ++j) {
        fprintf(sdbgfp, "  P%d pck count %d\n", cm->st_procs[j],
                cm->pck_cnt[j + 1] - cm->pck_cnt[j]);
        int l;
        for (l = cm->pck_cnt[j]; l < cm->pck_cnt[j + 1]; ++l)
            fprintf(sdbgfp, "    packet %d size %d, src %d, dest %d\n",
                    l - cm->pck_cnt[j] + 1, cm->pck_sizes[l], cm->pck_src[l],
                    cm->pck_dest[l]);
    }
#endif

    /* finalize copy data */
#ifdef D_SPEC_STFW
    fprintf(sdbgfp, "\n\ncopy data:\n");
    fflush(sdbgfp);
#endif
    for (i = 0; i < gm->nd; ++i) {
        stfw_comm *cm = &gm->sc[i];
        int tot_pcks = cm->pck_cnt[cm->st_np];
        cm->pck_cp = NULL;
#ifdef D_SPEC_STFW
        cm->pck_fw_idx = NULL;
#endif
        if (tot_pcks) {
            cm->pck_cp = malloc(sizeof(*cm->pck_cp) * tot_pcks);
#ifdef D_SPEC_STFW
            cm->pck_fw_idx = malloc(sizeof(*cm->pck_fw_idx) * tot_pcks);
#endif
        }

#ifdef D_SPEC_STFW
        fprintf(sdbgfp, "d = %d, #packets = %d\n", i, tot_pcks);
        fflush(sdbgfp);
#endif

        for (j = 0; j < tot_pcks; ++j) {
            if (pck_cp[i][j].d == -1) /* copy to vector */
            {
                cm->pck_cp[j] = pck_cp[i][j].ptr;

#ifdef D_SPEC_STFW
                cm->pck_fw_idx[j] = -1; /* mine */
#endif

#ifdef D_SPEC_STFW
                fprintf(sdbgfp,
                        "  packet %d -> copy to vector, "
                        "addr %p, fw_idx %d\n",
                        j, cm->pck_cp[j], cm->pck_fw_idx[j]);
                fflush(sdbgfp);
#endif
            } else /* copy forward buffer */
            {
                int d = pck_cp[i][j].d;
                int dim_idx = pck_cp[i][j].dim_idx;
                int fw_idx = pck_cp[i][j].fw_idx;
                assert(gm->sc[d].fw_imap[dim_idx] != NE);

                int iproc = gm->sc[d].fw_imap[dim_idx];

#ifdef D_SPEC_STFW
                fprintf(sdbgfp,
                        "  packet %d -> fw, in dim %d, dim_idx %d, "
                        "iproc %d, fw_idx %d, ",
                        j, d, dim_idx, iproc, fw_idx);
                fflush(sdbgfp);
#endif
                assert(fw_idx + gm->sc[d].fw_locs[iproc] <
                        gm->sc[d].fw_locs[iproc + 1]);

                fw_idx += gm->sc[d].fw_locs[iproc];
                cm->pck_cp[j] = &gm->sc[d].fw_buf[fw_idx];

#ifdef D_SPEC_STFW
                cm->pck_fw_idx[j] = fw_idx;
#endif

#ifdef D_SPEC_STFW
                fprintf(sdbgfp, "fw_idx with offset %d, addr %p, copy size %d\n", fw_idx,
                        cm->pck_cp[j], cm->pck_sizes[j]);
                fflush(sdbgfp);
#endif
            }
        }
    }

    /* cleanup */
    free(tmp_crds);
    for (i = 0; i < gm->nd; ++i) {
        free(fw_sz[i]);
        free(st_sz[i]);
        free(hdr_send[i]);
        free(hdr_recv[i]);
        free(hdr_reqs[i]);
        free(pck_cnt[i]);
        if (pck_cp[i])
            free(pck_cp[i]);
    }
    free(fw_sz);
    free(st_sz);
    free(fw_tot_sz);
    free(st_tot_sz);
    free(hdr_send);
    free(hdr_recv);
    free(hdr_reqs);
    free(pck_cnt);
    free(pck_cp);

    //gm_cstats(gm);
    free(inv_recv_map);

    return;
}

void init_stfw_nh_nbx(Comm *comm, vpt *gm) {
#define TAG_HDR 7

    int i, j, k;
    int rank, nprocs;
    int tmp;

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int *tmp_crds = malloc(sizeof(*tmp_crds) * gm->nd);

    /* data sizes communicated with each proc in each dimension  */
    int **fw_sz = malloc(sizeof(*fw_sz) * gm->nd);
    int **st_sz = malloc(sizeof(*st_sz) * gm->nd);
    int *fw_tot_sz = calloc(gm->nd, sizeof(*fw_tot_sz));
    int *st_tot_sz = calloc(gm->nd, sizeof(*st_tot_sz));

    int max_buf_sz = -1;

    /* headers */
    int **hdr_send = malloc(sizeof(*hdr_send) * gm->nd);
    int **hdr_recv = malloc(sizeof(*hdr_recv) * gm->nd);
    MPI_Request **hdr_reqs = malloc(sizeof(*hdr_reqs) * gm->nd);

    /* packets helpers */
    int **pck_cnt = malloc(sizeof(*pck_cnt) * gm->nd);
    struct pck_cp_data {
        /* if d = -1, packet is for me and dim_idx is the proc I recved data */
        /* if d != -1, use fw_idx, if it is -1 use ptr */
        int d;         /* which dim the packet will be forwarded in */
        int dim_idx;   /* index in the respective dim */
        int fw_idx;    /* forward index for fw_buf */
        vptcomm_real_t *ptr; /* direct address for vector p */
    };
    struct pck_cp_data **pck_cp = malloc(sizeof(*pck_cp) * gm->nd);

    for (i = 0; i < gm->nd; ++i) {
        int dim_size = gm->size[i];
        /* each packet: src dest size */
        int hlen = 3 * (nprocs / dim_size) + 1;
        hdr_send[i] = calloc(dim_size * hlen, sizeof(**hdr_send));
        hdr_recv[i] = calloc(dim_size * hlen, sizeof(**hdr_recv));
        hdr_reqs[i] = malloc(sizeof(**hdr_reqs) * dim_size);

        fw_sz[i] = calloc(dim_size, sizeof(**fw_sz));
        st_sz[i] = calloc(dim_size, sizeof(**st_sz));

        pck_cnt[i] = calloc(dim_size, sizeof(**pck_cnt));
    }

#ifdef D_SPEC_STFW
    fprintf(sdbgfp, "> init_stfw_gmesh\n");
    fprintf(sdbgfp, "number of processors to communicate = %d\n", comm->sno);
    for (i = 0; i < comm->sno; ++i)
        fprintf(sdbgfp, "  P %d, send size = %d, recv size = %d\n", gm->map[comm->sprocs[i]],
                comm->ssizes[i], comm->rsizes[i]);
    fprintf(sdbgfp, "\n");
    fflush(sdbgfp);
#endif

    /* fill inv_recv_map */
    int *inv_recv_map = malloc(sizeof(*inv_recv_map) * nprocs);
    for (i = 0; i < comm->rno; ++i) {
        inv_recv_map[gm->map[comm->rprocs[i]]] = gm->map[i];
    }
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

        /* add the header */
        int hlen = 3 * (nprocs / gm->size[d]) + 1;
        int proc_crd = get_crd(gm, proc, d);
        int *ptr = &hdr_send[d][hlen * proc_crd];
        int packet_id = *ptr;
        ptr += 3 * packet_id + 1;
        *ptr = gm->map[rank];                  /* src */
        *(ptr + 1) = proc;            /* dest */
        assert(comm->ssizes[i] > 0);
        *(ptr + 2) = comm->ssizes[i]; /* size */
        hdr_send[d][hlen * proc_crd] += 1;
        assert(packet_id >= 0 && packet_id < (nprocs / gm->size[d]));

#ifdef D_SPEC_STFW
        fprintf(sdbgfp,
                "  adding to respective header, this is packet %d, "
                "size %d, proc crd in respective dim is %d, max packets %d\n",
                hdr_send[d][hlen * proc_crd], comm->ssizes[i], proc_crd,
                nprocs / gm->size[d]);
        fflush(sdbgfp);
#endif

        /* update forward communicators */
        if (fw_sz[d][proc_crd] == 0) /* first packet */
            ++(gm->sc[d].fw_np);
        fw_sz[d][proc_crd] += comm->ssizes[i];
        fw_tot_sz[d] += comm->ssizes[i];

#ifdef D_SPEC_STFW
        fprintf(sdbgfp,
                "  [comm fw] fw_sz and fw_np: fw_sz[%d][%d] = %d, "
                "gm->sc[%d].fw_np = %d, fw_tot_sz[%d] = %d\n",
                d, proc_crd, fw_sz[d][proc_crd], d, gm->sc[d].fw_np, d,
                fw_tot_sz[d]);
#endif
    }

#ifdef D_SPEC_STFW
    print_hdr(gm, hdr_send, hdr_recv);
#endif

    for (i = 0; i < gm->nd; ++i) {
        int curd = gm->order[i];
        int my_crd = gm->crd[curd];
        int cur_hlen = 3 * (nprocs / gm->size[curd]) + 1;
        int *cur_hsend = hdr_send[curd];
        int *cur_hrecv = hdr_recv[curd];
        MPI_Request *cur_hreqs = hdr_reqs[curd];


#ifdef D_SPEC_STFW
        fprintf(sdbgfp, "\nissuing sends...\n");
        fflush(sdbgfp);
#endif

        int nsendto=0;
        for (j = 0; j < gm->size[curd]; ++j) {
            if (j != my_crd && cur_hsend[j*cur_hlen] > 0) {
#ifdef D_SPEC_STFW
                fprintf(sdbgfp,
                        "  sending to P%d (crd %d at dim %d), "
                        "npackets = %d, send size = %d\n",
                        get_nghbr(gm, curd, j), j, curd, cur_hsend[j * cur_hlen],
                        cur_hsend[j * cur_hlen] * 3 + 1);
                fflush(sdbgfp);
#endif

                MPI_Issend(&cur_hsend[j * cur_hlen], cur_hsend[j * cur_hlen] * 3 + 1, MPI_INT, gm->inv_map[get_nghbr(gm, curd, j)], TAG_HDR, MPI_COMM_WORLD, &cur_hreqs[nsendto]);
                ++nsendto;
            }
        }

#ifdef D_SPEC_STFW
        fprintf(sdbgfp, "\n\n");
        fprintf(sdbgfp,
                "comm operations turn %d, dimension %d, my_crd = %d, "
                "hlen = %d, max #packets = %d\n",
                i, curd, my_crd, cur_hlen, nprocs / gm->size[curd]);
        fflush(sdbgfp);
#endif


        /* process hdr_recv of the prev step to form send data of this step */
            int done=0, barrier_act=0, flag, tot_pck_cnt = 0;
            MPI_Status stts;
            MPI_Request ibrrqst; 
            int *inds = malloc(sizeof(*inds) * gm->size[curd]);
            /* make sure cur_hrecv counts are zero before recv */
            memset(cur_hrecv, 0, sizeof(int) * gm->size[curd]*cur_hlen);
            while(!done){
                MPI_Iprobe(MPI_ANY_SOURCE, TAG_HDR, MPI_COMM_WORLD, &flag, &stts);
                if(flag){
                    int source = stts.MPI_SOURCE, nrcvdhdrs; MPI_Get_count(&stts, MPI_INT, &nrcvdhdrs);
                    int rcrd=get_crd(gm, gm->map[source], curd);
                    MPI_Recv(&cur_hrecv[rcrd * cur_hlen], nrcvdhdrs, MPI_INT, source, TAG_HDR, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    int npackets = cur_hrecv[rcrd * cur_hlen];
                    assert(npackets > 0 && npackets <= nprocs / gm->size[curd]);


                    tot_pck_cnt += npackets;
                    pck_cnt[curd][rcrd] = npackets;
                }
                if(barrier_act){ int ibrflag; MPI_Test(&ibrrqst, &ibrflag, MPI_STATUS_IGNORE); if(ibrflag) done=1; }
                else{
                    int testflag=0;
                    MPI_Testall(nsendto, cur_hreqs, &testflag, MPI_STATUSES_IGNORE);
                    if(testflag) {MPI_Ibarrier(MPI_COMM_WORLD, &ibrrqst); barrier_act=1;}
                }
            }
            free(inds);

#ifdef D_SPEC_STFW
            fprintf(sdbgfp, "  total recved packet count %d\n", tot_pck_cnt);
            fflush(sdbgfp);
#endif

            stfw_comm *cm = &gm->sc[curd];
            cm->pck_sizes = NULL;
            pck_cp[curd] = NULL;
#ifdef D_SPEC_STFW
            cm->pck_dest = NULL;
            cm->pck_src = NULL;
#endif
            if (tot_pck_cnt) {
                cm->pck_sizes = malloc(sizeof(*cm->pck_sizes) * tot_pck_cnt);
                pck_cp[curd] = malloc(sizeof(**pck_cp) * tot_pck_cnt);
#ifdef D_SPEC_STFW
                cm->pck_dest = malloc(sizeof(*cm->pck_dest) * tot_pck_cnt);
                cm->pck_src = malloc(sizeof(*cm->pck_src) * tot_pck_cnt);
#endif
            }

            /* process packets */
            int pck_idx = 0;
            for (j = 0; j < gm->size[curd]; ++j) {
                if (j == my_crd)
                    continue;

                int npackets = cur_hrecv[j * cur_hlen];
                int *ptr = NULL;
                if (npackets > 0) /* non-empty recv */
                    ptr = &cur_hrecv[j *cur_hlen + 1];
                
                /* update store communicators */
                if (npackets > 0 && st_sz[curd][j] == 0) /* first packet */
                    ++(gm->sc[curd].st_np);

#ifdef D_SPEC_STFW
                fprintf(sdbgfp,
                        "  recv from P%d (crd %d at dim %d), "
                        "npackets = %d, contents:\n",
                        get_nghbr(gm, curd, j), j, curd, npackets);
                fflush(sdbgfp);
#endif

                for (k = 0; k < npackets; ++k, ptr += 3) {
                    int pck_src = *ptr;
                    int pck_dest = *(ptr + 1);
                    int pck_size = *(ptr + 2);
                    assert(pck_size > 0);
                    cm->pck_sizes[pck_idx] = pck_size;
#ifdef D_SPEC_STFW
                    cm->pck_dest[pck_idx] = pck_dest;
                    cm->pck_src[pck_idx] = pck_src;
#endif

#ifdef D_SPEC_STFW
                    fprintf(sdbgfp, "    P%d P%d %d => ", pck_src, pck_dest, pck_size);
                    fflush(sdbgfp);
#endif

                    /* update store communicators */
                    st_sz[curd][j] += pck_size;
                    st_tot_sz[curd] += pck_size;

                    /* packet is for me */
                    if (pck_dest == gm->map[rank]) {
#ifdef D_SPEC_STFW
                        assert(inv_recv_map[pck_src] != NE);
                        fprintf(sdbgfp,
                                "package is for me, pck_idx = %d, "
                                "tuple info d = -1, dim_dix = %d, ptr = %p\n",
                                pck_idx, pck_src, comm->recv[gm->inv_map[inv_recv_map[pck_src]]]);

                        /* also assert recved proc rank is in my recv list */
                        int l;
                        int rsize = -1;
                        for (l = 0; l < comm->rno; ++l) {
                            if (comm->rprocs[l] == gm->inv_map[pck_src]) {
                                rsize = comm->rsizes[l];
                                break;
                            }
                        }
                        fprintf(sdbgfp, "rsize=%d pck_size=%d\n", rsize, pck_size);
                        fflush(sdbgfp);
                        assert(l >= 0 && l < comm->rno && rsize == pck_size);
                        fprintf(sdbgfp,
                                "    found packet in my recv list %d = P%d, "
                                "size = %d\n",
                                l, comm->rprocs[l], rsize);
                        fflush(sdbgfp);
#endif
                        pck_cp[curd][pck_idx].d = -1;
                        pck_cp[curd][pck_idx].dim_idx = pck_src;
                        pck_cp[curd][pck_idx++].ptr =
                            comm->recv[gm->inv_map[inv_recv_map[pck_src]]];
                        continue;
                    }

                    /* packet is to be forwarded */
                    int d = gm->comm_dim[pck_dest];
                    int hlen = 3 * (nprocs / gm->size[d]) + 1;
                    int proc_crd = get_crd(gm, pck_dest, d);
                    int *ptr2 = &hdr_send[d][hlen * proc_crd];
                    int packet_id = *ptr2;
                    ptr2 += 3 * packet_id + 1;
                    *ptr2 = pck_src;
                    *(ptr2 + 1) = pck_dest;
                    *(ptr2 + 2) = pck_size;
                    hdr_send[d][hlen * proc_crd] += 1;
                    assert(packet_id >= 0 && packet_id < (nprocs / gm->size[d]));

#ifdef D_SPEC_STFW
                    fprintf(sdbgfp,
                            "adding to the comm in dim %d, proc crd %d, "
                            "packet id %d, size %d, max packets %d\n",
                            d, proc_crd, hdr_send[d][hlen * proc_crd], pck_size,
                            nprocs / gm->size[d]);
                    fprintf(sdbgfp,
                            "    package will be forwarded, "
                            "pck_idx = %d, tuple info d = %d, dim_dix = %d, "
                            "fw_idx = %d\n",
                            pck_idx, d, proc_crd, fw_sz[d][proc_crd]);
                    fflush(sdbgfp);
#endif

                    pck_cp[curd][pck_idx].d = d;
                    pck_cp[curd][pck_idx].dim_idx = proc_crd;
                    pck_cp[curd][pck_idx++].fw_idx = fw_sz[d][proc_crd];

                    /* update forward communicators */
                    if (fw_sz[d][proc_crd] == 0) /* first packet */
                        ++(gm->sc[d].fw_np);
                    fw_sz[d][proc_crd] += pck_size;
                    fw_tot_sz[d] += pck_size;

#ifdef D_SPEC_STFW
                    fprintf(sdbgfp,
                            "    [comm fw] fw_sz and fw_np: "
                            "fw_sz[%d][%d] = %d, gm->sc[%d].fw_np = %d, "
                            "fw_tot_sz[%d] = %d\n",
                            d, proc_crd, fw_sz[d][proc_crd], d, gm->sc[d].fw_np, d,
                            fw_tot_sz[d]);
#endif
                }

#ifdef D_SPEC_STFW
                fprintf(sdbgfp,
                        "    [comm st] st_sz and st_np: "
                        "st_sz[%d][%d] = %d, gm->sc[%d].st_np = %d, "
                        "st_tot_sz[%d] = %d\n",
                        curd, j, st_sz[curd][j], curd,
                        gm->sc[curd].st_np, curd, st_tot_sz[curd]);
                fflush(sdbgfp);
#endif
            }

            /* now can finalize the previous dimension's communicators */
            arrange_comm_nh(gm, curd, fw_sz[curd], fw_tot_sz[curd],
                    st_sz[curd], st_tot_sz[curd], pck_cnt[curd]);
            if (fw_tot_sz[curd] + st_tot_sz[curd] > max_buf_sz)
                max_buf_sz = fw_tot_sz[curd] + st_tot_sz[curd];

#ifdef D_SPEC_STFW
            assert(pck_idx == tot_pck_cnt);
            fprintf(sdbgfp, "recv packet info:\n");
            for (j = 0; j < cm->st_np; ++j) {
                fprintf(sdbgfp, "  P%d pck count %d\n", cm->st_procs[j],
                        cm->pck_cnt[j + 1] - cm->pck_cnt[j]);
                int l;
                for (l = cm->pck_cnt[j]; l < cm->pck_cnt[j + 1]; ++l)
                    fprintf(sdbgfp, "    packet %d size %d, src %d, dest %d\n",
                            l - cm->pck_cnt[j] + 1, cm->pck_sizes[l], cm->pck_src[l],
                            cm->pck_dest[l]);
            }
#endif

    }




    /* finalize copy data */
#ifdef D_SPEC_STFW
    fprintf(sdbgfp, "\n\ncopy data:\n");
    fflush(sdbgfp);
#endif
    for (i = 0; i < gm->nd; ++i) {
        stfw_comm *cm = &gm->sc[i];
        int tot_pcks = cm->pck_cnt[cm->st_np];
        cm->pck_cp = NULL;
#ifdef D_SPEC_STFW
        cm->pck_fw_idx = NULL;
#endif
        if (tot_pcks) {
            cm->pck_cp = malloc(sizeof(*cm->pck_cp) * tot_pcks);
#ifdef D_SPEC_STFW
            cm->pck_fw_idx = malloc(sizeof(*cm->pck_fw_idx) * tot_pcks);
#endif
        }

#ifdef D_SPEC_STFW
        fprintf(sdbgfp, "d = %d, #packets = %d\n", i, tot_pcks);
        fflush(sdbgfp);
#endif

        for (j = 0; j < tot_pcks; ++j) {
            if (pck_cp[i][j].d == -1) /* copy to vector */
            {
                cm->pck_cp[j] = pck_cp[i][j].ptr;

#ifdef D_SPEC_STFW
                cm->pck_fw_idx[j] = -1; /* mine */
#endif

#ifdef D_SPEC_STFW
                fprintf(sdbgfp,
                        "  packet %d -> copy to vector, "
                        "addr %p, fw_idx %d\n",
                        j, cm->pck_cp[j], cm->pck_fw_idx[j]);
                fflush(sdbgfp);
#endif
            } else /* copy forward buffer */
            {
                int d = pck_cp[i][j].d;
                int dim_idx = pck_cp[i][j].dim_idx;
                int fw_idx = pck_cp[i][j].fw_idx;
                assert(gm->sc[d].fw_imap[dim_idx] != NE);

                int iproc = gm->sc[d].fw_imap[dim_idx];

#ifdef D_SPEC_STFW
                fprintf(sdbgfp,
                        "  packet %d -> fw, in dim %d, dim_idx %d, "
                        "iproc %d, fw_idx %d, ",
                        j, d, dim_idx, iproc, fw_idx);
                fflush(sdbgfp);
#endif
                assert(fw_idx + gm->sc[d].fw_locs[iproc] <
                        gm->sc[d].fw_locs[iproc + 1]);

                fw_idx += gm->sc[d].fw_locs[iproc];
                cm->pck_cp[j] = &gm->sc[d].fw_buf[fw_idx];

#ifdef D_SPEC_STFW
                cm->pck_fw_idx[j] = fw_idx;
#endif

#ifdef D_SPEC_STFW
                fprintf(sdbgfp, "fw_idx with offset %d, addr %p, copy size %d\n", fw_idx,
                        cm->pck_cp[j], cm->pck_sizes[j]);
                fflush(sdbgfp);
#endif
            }
        }
    }

    /* cleanup */
    free(tmp_crds);
    for (i = 0; i < gm->nd; ++i) {
        free(fw_sz[i]);
        free(st_sz[i]);
        free(hdr_send[i]);
        free(hdr_recv[i]);
        free(hdr_reqs[i]);
        free(pck_cnt[i]);
        if (pck_cp[i])
            free(pck_cp[i]);
    }
    free(fw_sz);
    free(st_sz);
    free(fw_tot_sz);
    free(st_tot_sz);
    free(hdr_send);
    free(hdr_recv);
    free(hdr_reqs);
    free(pck_cnt);
    free(pck_cp);

    //gm_cstats(gm);
    free(inv_recv_map);

    return;
}

void init_stfw(Comm *comm, vpt *gm, int STFW_INIT_STRATEGY) {
    switch (STFW_INIT_STRATEGY) {
        case DRYRUN_HEAVY_P2P:
            init_stfw_nh_heavy(comm, gm);
            break;
        case DRYRUN_SPARSE_P2P:
            init_stfw_nh_nbx(comm, gm);
    }
}



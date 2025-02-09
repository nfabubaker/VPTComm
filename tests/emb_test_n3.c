/**
 * @author      : nabeelooo (nabeelooo@$HOSTNAME)
 * @file        : emb_test
 * @created     : Saturday Jun 19, 2021 11:00:48 +03
 */

#include "../include/vptcomm.h"
#include "../src/util.h"
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <sys/stat.h>
#include <string.h>

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    STFW_init(1, 0);

    int i, mypid, npes, nsendwho,  *sendwho,  *xsendind,  *sendind,  nrecvwho,  *recvwho,  *xrecvind,  *recvind,  ndims, *indsmap,  embDataUnitSize;  
    vptcomm_real_t *data = malloc(sizeof(*data) * 4 * 100);
    
    MPI_Comm_rank(MPI_COMM_WORLD, &mypid);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);

    ndims = log2(npes);
    for (i = 0; i < 400; ++i) {
        data[i] = ((int)((i/40)) == mypid ? (vptcomm_real_t)mypid/10:i) ;
    }
#ifdef NA_DBG
    struct stat st = {0};

    if (stat("./dbg_logs", &st) == -1) {
        mkdir("./dbg_logs", 0700);
    }
    sprintf(dbg_fn, "./dbg_logs/outfile-%d-%d", mypid, npes);
    dbgfp = fopen(dbg_fn, "w");
#endif 

    switch (mypid) {
        case 0:
            nsendwho = 2; nrecvwho = 2;
            sendwho = malloc(sizeof(int) * nsendwho); sendwho[0] = 2; sendwho[1] = 7;
            recvwho = malloc(sizeof(int) * nrecvwho); recvwho[0] = 1; recvwho[1] = 3;
            xsendind = calloc(npes+1,sizeof(int)); xsendind[3] = 3; xsendind[8] = 3; 
            xrecvind = calloc(npes+1, sizeof(int)); xrecvind[2] = 3; xrecvind[4] = 3;
            sendind = malloc(sizeof(int) * 3); sendind[0] = 1; sendind[1] = 2; sendind[2] = 3;  sendind[3] = 4; sendind[4] = 5; sendind[5] = 6;
            recvind = malloc(sizeof(int) * 6); recvind[0] = 11; recvind[1] = 12; recvind[2] = 13; recvind[3] = 31; recvind[4] = 32; recvind[5] = 33;
            break;

        case 1:
            nrecvwho = 1; nsendwho = 3;
            recvwho = malloc(sizeof(int) * nrecvwho); recvwho[0] = 2;
            sendwho = malloc(sizeof(int) * nsendwho); sendwho[0] = 0; sendwho[1] = 2; sendwho[2] = 3;
            xsendind = calloc(npes+1,sizeof(int)); xsendind[1] = 3; xsendind[3] = 3, xsendind[4] = 3;
            xrecvind = calloc(npes+1, sizeof(int)); xrecvind[3] = 3;
            recvind = malloc(sizeof(int) * 3); recvind[0] = 21; recvind[1] = 22; recvind[2] = 23; 
            sendind = malloc(sizeof(int) * 9); sendind[0] = 11; sendind[1] = 12; sendind[2] = 13; sendind[3] = 14; sendind[4] = 15; sendind[5] = 16;sendind[6] = 17; sendind[7] = 18;sendind[8] = 19;

            break;
        case 2: 
            nsendwho = 1; nrecvwho = 2;
            sendwho = malloc(sizeof(int) * nsendwho); sendwho[0] = 1;
            recvwho = malloc(sizeof(int) * nrecvwho); recvwho[0] = 0; recvwho[1] = 1;
            xsendind = calloc(npes+1,sizeof(int)); xsendind[2] = 3; 
            xrecvind = calloc(npes+1, sizeof(int)); xrecvind[1] = 3; xrecvind[2] = 3;
            sendind = malloc(sizeof(int) * 3); sendind[0] = 21; sendind[1] = 22; sendind[2] = 23; 
            recvind = malloc(sizeof(int) * 6); recvind[0] = 1; recvind[1] = 2; recvind[2] = 3; recvind[3] = 14; recvind[4] = 15; recvind[5] = 16;
            break;
        case 3:
            nsendwho = 1; nrecvwho = 1;
            sendwho = malloc(sizeof(int) * nsendwho); sendwho[0] = 0;
            recvwho = malloc(sizeof(int) * nrecvwho); recvwho[0] = 1;
            xsendind = calloc(npes+1,sizeof(int)); xsendind[1] = 3; 
            xrecvind = calloc(npes+1, sizeof(int)); xrecvind[2] = 3; 
            sendind = malloc(sizeof(int) * 3); sendind[0] = 31; sendind[1] = 32; sendind[2] = 33; 
            recvind = malloc(sizeof(int) * 6); recvind[0] = 17; recvind[1] = 18; recvind[2] = 19;
            break;
        case 4:
            nsendwho = 1; nrecvwho = 2;
            sendwho = malloc(sizeof(int) * nsendwho); sendwho[0] = 6;
            recvwho = malloc(sizeof(int) * nrecvwho); recvwho[0] = 5; recvwho[1] = 7;
            xsendind = calloc(npes+1,sizeof(int)); xsendind[7] = 3; 
            xrecvind = calloc(npes+1, sizeof(int)); xrecvind[6] = 3; xrecvind[8] = 3;
            sendind = malloc(sizeof(int) * 3); sendind[0] = 41; sendind[1] = 42; sendind[2] = 43; 
            recvind = malloc(sizeof(int) * 6); recvind[0] = 51; recvind[1] = 52; recvind[2] = 53; recvind[3] = 71; recvind[4] = 72; recvind[5] = 73;
            break;

        case 5:
            nrecvwho = 1; nsendwho = 2;
            recvwho = malloc(sizeof(int) * nrecvwho); recvwho[0] = 6;
            sendwho = malloc(sizeof(int) * nsendwho); sendwho[0] = 4; sendwho[1] = 6;
            xsendind = calloc(npes+1,sizeof(int)); xsendind[5] = 3; xsendind[7]=3; 
            xrecvind = calloc(npes+1, sizeof(int)); xrecvind[7] = 3;
            recvind = malloc(sizeof(int) * 3); recvind[0] = 61; recvind[1] = 62; recvind[2] = 63; 
            sendind = malloc(sizeof(int) * 9); sendind[0] = 51; sendind[1] = 52; sendind[2] = 53; sendind[3] = 54; sendind[4] = 55; sendind[5] = 56;
            break;
        case 6: 
            nsendwho = 1; nrecvwho = 2;
            sendwho = malloc(sizeof(int) * nsendwho); sendwho[0] = 5;
            recvwho = malloc(sizeof(int) * nrecvwho); recvwho[0] = 4; recvwho[1] = 5;
            xsendind = calloc(npes+1,sizeof(int)); xsendind[6] = 3; 
            xrecvind = calloc(npes+1, sizeof(int)); xrecvind[5] = 3; xrecvind[6] = 3;
            sendind = malloc(sizeof(int) * 3); sendind[0] = 61; sendind[1] = 62; sendind[2] = 63; 
            recvind = malloc(sizeof(int) * 6); recvind[0] = 41; recvind[1] = 42; recvind[2] = 43; recvind[3] = 54; recvind[4] = 55; recvind[5] = 56;
            break;
        case 7:
            nsendwho = 1; nrecvwho = 1;
            sendwho = malloc(sizeof(int) * nsendwho); sendwho[0] = 4;
            recvwho = malloc(sizeof(int) * nrecvwho); recvwho[0] = 0;
            xsendind = calloc(npes+1,sizeof(int)); xsendind[5] = 3; 
            xrecvind = calloc(npes+1, sizeof(int)); xrecvind[1] = 3; 
            sendind = malloc(sizeof(int) * 3); sendind[0] = 71; sendind[1] = 72; sendind[2] = 73; 
            recvind = malloc(sizeof(int) * 6); recvind[0] = 4; recvind[1] = 5; recvind[2] = 6;
            break;
        default:
            break;

    }
    xsendind[0] = 0;
    xrecvind[0] = 0;
    for (i = 1; i < npes+1; ++i) {
        xsendind[i] += xsendind[i-1];   
        xrecvind[i] += xrecvind[i-1];
    }
    printf( "%s\n", "done init");
    embDataUnitSize = 4; 
    indsmap = malloc(sizeof(int) * 100);
    for (i = 0; i < 100; ++i) {
        indsmap[i] = i; /* 1-1 mappping */
    }
    int * HI = malloc(sizeof(*HI) * npes);
    for (i = 0; i < npes; ++i) {
       //HI[i] = npes-i; 
       HI[i] = i;
    }
    
    vptcomm_real_t *sbuff, *recvbuff, **sendp, **recvp;
    sbuff = malloc(sizeof(vptcomm_real_t) * xsendind[npes] * embDataUnitSize);
    recvbuff = malloc(sizeof(vptcomm_real_t) * xrecvind[npes] * embDataUnitSize);
    int j;
    vptcomm_real_t *ptr = sbuff;
    for (i = 0; i < nsendwho; ++i){
        int sid=sendwho[i]; 
        for (j = xsendind[sid]; j < xsendind[sid+1]; ++j) {
            memcpy(ptr, &data[sendind[j] * embDataUnitSize], sizeof(*data)*embDataUnitSize);
            ptr+= embDataUnitSize;
        }
    }
    sendp = malloc(sizeof(*sendp) * nsendwho); 
    recvp = malloc(sizeof(*recvp) * nrecvwho); 
    int *ssend, *srecv;
    ssend = malloc(sizeof(*ssend) * nsendwho);
    srecv = malloc(sizeof(*srecv) * nrecvwho);
    for (i = 0; i < nsendwho; ++i) {ssend[i]=(xsendind[sendwho[i]+1]-xsendind[sendwho[i]])*embDataUnitSize; sendp[i] = &sbuff[xsendind[sendwho[i]]*embDataUnitSize]; }
    for (i = 0; i < nrecvwho; ++i) { srecv[i]=(xrecvind[recvwho[i]+1]-xrecvind[recvwho[i]])*embDataUnitSize; recvp[i] = &recvbuff[xrecvind[recvwho[i]]*embDataUnitSize]; }
    printf( "%s\n", "done setup");


    //HI[0] = 2; HI[1] = 7; HI[2] = 4; HI[3] = 5; HI[4] = 0; HI[5] = 6; HI[6] = 1; HI[7] = 3;
    int dsizes[2]={2,4};
    STFW_init_instance(0, 2, dsizes, nsendwho, sendwho, nrecvwho, recvwho, ssend, sendp, srecv,recvp, HI);
    printf( "%s\n", "done stfw init");

    MPI_Barrier(MPI_COMM_WORLD);
/*     {
 *         volatile int tt = 0;
 *         printf("PID %d on %d ready for attach\n", mypid,  getpid());
 *         fflush(stdout);
 *         while (0 == tt)
 *             sleep(5);
 *     }
 */
/*     int dim, ct_cnt;
 *     for (dim = 0; dim < ec->ndims; ++dim) {
 *         for (ct_cnt = 1; ct_cnt < CT_CNT; ++ct_cnt) {
 *             for (i = ec->xsendptrs[dim][ct_cnt]; i < ec->xsendptrs[dim][ct_cnt+1]; ++i) {
 *                 printf( "\t\t*sendptr[%d][%d] = %f", dim, i , *ec->sendptrs[dim][i]);	
 *             }	
 *         }
 *     }
 *     for (dim = ndims-1; dim >=0; --dim) {
 *         for (ct_cnt = 1; ct_cnt < CT_CNT; ++ct_cnt) {
 *             for (i = ec_d->xsendptrs[dim][ct_cnt]; i < ec_d->xsendptrs[dim][ct_cnt+1]; ++i) {
 *                 printf( "\t\t[dual]*sendptr[%d][%d] = %f", dim, i , *ec_d->sendptrs[dim][i]);	
 *             }	
 *         }
 *     }
 */

    STFW_Comm(0);
    MPI_Barrier(MPI_COMM_WORLD);
    printf( "STFW_comm done\n");


    ptr = recvbuff;
    for (i = 0; i < nrecvwho; ++i){
        int sid=recvwho[i]; 
        for (j = xrecvind[sid]; j < xrecvind[sid+1]; ++j) {
            memcpy(&data[recvind[j]*embDataUnitSize], ptr, sizeof(*data)*embDataUnitSize);
            ptr+= embDataUnitSize;
        }
    }

    for (i = 0; i < 100; ++i) {
        printf( "[%d] ", i);
        for (j = 0; j < 4; ++j) {
            printf( "%0.2f ", data[(i*4)+j]);
        }
        printf( "\n");
    }
    printf( "%s\n", "done comm");
    STFW_finalize();
    MPI_Finalize();
    return 0;
}

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
#include <assert.h>

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    STFW_init(1, 0);

    int i, mypid, npes, nsendwho,  *sendwho,  *ssize,  **sendind,  nrecvwho,  *recvwho,  *rsize,  **recvind,  ndims, *indsmap,  embDataUnitSize;  
    
    MPI_Comm_rank(MPI_COMM_WORLD, &mypid);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);

    ndims = log2(npes);

    embDataUnitSize = 4; 
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
            ssize = calloc(nsendwho,sizeof(int)); rsize = calloc(nrecvwho,sizeof(int));
            sendind = malloc(sizeof(int*) * nsendwho); recvind = malloc(sizeof(int*) * nrecvwho);  
            for(i=0; i < nsendwho; ++i) {sendind[i] = malloc(sizeof(int) * 3); ssize[i] = 3*embDataUnitSize;}
            for(i=0; i < nrecvwho; ++i) {recvind[i] = malloc(sizeof(int) * 3); rsize[i] = 3*embDataUnitSize;}
            sendind[0][0] = 1; sendind[0][1] = 2; sendind[0][2] = 3;  sendind[1][0] = 1; sendind[1][1] = 2; sendind[1][2] = 3;
            recvind[0][0] = 11; recvind[0][1] = 12; recvind[0][2] = 13; recvind[1][0] = 31; recvind[1][1] = 32; recvind[1][2] = 33;
            break;

        case 1:
            nrecvwho = 1; nsendwho = 3;
            recvwho = malloc(sizeof(int) * nrecvwho); recvwho[0] = 2;
            sendwho = malloc(sizeof(int) * nsendwho); sendwho[0] = 0; sendwho[1] = 2; sendwho[2] = 3;
            ssize = calloc(nsendwho,sizeof(int)); rsize = calloc(nrecvwho,sizeof(int));
            sendind = malloc(sizeof(int*) * nsendwho); recvind = malloc(sizeof(int*) * nrecvwho);  
            for(i=0; i < nsendwho; ++i) {sendind[i] = malloc(sizeof(int) * 3); ssize[i] = 3*embDataUnitSize;}
            for(i=0; i < nrecvwho; ++i) {recvind[i] = malloc(sizeof(int) * 3); rsize[i] = 3*embDataUnitSize;}
            sendind[0][0] = 11; sendind[0][1] = 12; sendind[0][2] = 13;  sendind[1][0] = 14; sendind[1][1] = 15; sendind[1][2] = 16;sendind[2][0] = 17; sendind[2][1] = 18; sendind[2][2] = 19;
            recvind[0][0] = 21; recvind[0][1] = 22; recvind[0][2] = 23; 

            break;
        case 2: 
            nsendwho = 1; nrecvwho = 2;
            sendwho = malloc(sizeof(int) * nsendwho); sendwho[0] = 1;
            recvwho = malloc(sizeof(int) * nrecvwho); recvwho[0] = 0; recvwho[1] = 1;
            ssize = calloc(nsendwho,sizeof(int)); rsize = calloc(nrecvwho,sizeof(int));
            sendind = malloc(sizeof(int*) * nsendwho); recvind = malloc(sizeof(int*) * nrecvwho);  
            for(i=0; i < nsendwho; ++i) {sendind[i] = malloc(sizeof(int) * 3); ssize[i] = 3*embDataUnitSize;}
            for(i=0; i < nrecvwho; ++i) {recvind[i] = malloc(sizeof(int) * 3); rsize[i] = 3*embDataUnitSize;}
            sendind[0][0] = 21; sendind[0][1] = 22; sendind[0][2] = 23;  
            recvind[0][0] = 1; recvind[0][1] = 2; recvind[0][2] = 3; recvind[1][0] = 14; recvind[1][1] = 15; recvind[1][2] = 16; 
            break;
        case 3:
            nsendwho = 1; nrecvwho = 1;
            sendwho = malloc(sizeof(int) * nsendwho); sendwho[0] = 0;
            recvwho = malloc(sizeof(int) * nrecvwho); recvwho[0] = 1;
            ssize = calloc(nsendwho,sizeof(int)); rsize = calloc(nrecvwho,sizeof(int));
            sendind = malloc(sizeof(int*) * nsendwho); recvind = malloc(sizeof(int*) * nrecvwho);  
            for(i=0; i < nsendwho; ++i) {sendind[i] = malloc(sizeof(int) * 3); ssize[i] = 3*embDataUnitSize;}
            for(i=0; i < nrecvwho; ++i) {recvind[i] = malloc(sizeof(int) * 3); rsize[i] = 3*embDataUnitSize;}
            sendind[0][0] = 31; sendind[0][1] = 32; sendind[0][2] = 33;  
            recvind[0][0] = 17; recvind[0][1] = 18; recvind[0][2] = 19; 
            break;
        case 4:
            nsendwho = 1; nrecvwho = 2;
            sendwho = malloc(sizeof(int) * nsendwho); sendwho[0] = 6;
            recvwho = malloc(sizeof(int) * nrecvwho); recvwho[0] = 5; recvwho[1] = 7;
            ssize = calloc(nsendwho,sizeof(int)); rsize = calloc(nrecvwho,sizeof(int));
            sendind = malloc(sizeof(int*) * nsendwho); recvind = malloc(sizeof(int*) * nrecvwho);  
            for(i=0; i < nsendwho; ++i) {sendind[i] = malloc(sizeof(int) * 3); ssize[i] = 3*embDataUnitSize;}
            for(i=0; i < nrecvwho; ++i) {recvind[i] = malloc(sizeof(int) * 3); rsize[i] = 3*embDataUnitSize;}
            sendind[0][0] = 41; sendind[0][1] = 42; sendind[0][2] = 43;  
            recvind[0][0] = 51; recvind[0][1] = 52; recvind[0][2] = 53; recvind[1][0] = 71; recvind[1][1] = 72; recvind[1][2] = 73; 
            break;

        case 5:
            nrecvwho = 1; nsendwho = 2;
            recvwho = malloc(sizeof(int) * nrecvwho); recvwho[0] = 6;
            sendwho = malloc(sizeof(int) * nsendwho); sendwho[0] = 4; sendwho[1] = 6;
            ssize = calloc(nsendwho,sizeof(int)); rsize = calloc(nrecvwho,sizeof(int));
            sendind = malloc(sizeof(int*) * nsendwho); recvind = malloc(sizeof(int*) * nrecvwho);  
            for(i=0; i < nsendwho; ++i) {sendind[i] = malloc(sizeof(int) * 3); ssize[i] = 3*embDataUnitSize;}
            for(i=0; i < nrecvwho; ++i) {recvind[i] = malloc(sizeof(int) * 3); rsize[i] = 3*embDataUnitSize;}
            sendind[0][0] = 51; sendind[0][1] = 52; sendind[0][2] = 53;  sendind[1][0] = 51; sendind[1][1] = 52; sendind[1][2] = 53;
            recvind[0][0] = 61; recvind[0][1] = 62; recvind[0][2] = 63; 
            break;
        case 6: 
            nsendwho = 1; nrecvwho = 2;
            sendwho = malloc(sizeof(int) * nsendwho); sendwho[0] = 5;
            recvwho = malloc(sizeof(int) * nrecvwho); recvwho[0] = 4; recvwho[1] = 5;
            ssize = calloc(nsendwho,sizeof(int)); rsize = calloc(nrecvwho,sizeof(int));
            sendind = malloc(sizeof(int*) * nsendwho); recvind = malloc(sizeof(int*) * nrecvwho);  
            for(i=0; i < nsendwho; ++i) {sendind[i] = malloc(sizeof(int) * 3); ssize[i] = 3*embDataUnitSize;}
            for(i=0; i < nrecvwho; ++i) {recvind[i] = malloc(sizeof(int) * 3); rsize[i] = 3*embDataUnitSize;}
            sendind[0][0] = 61; sendind[0][1] = 62; sendind[0][2] = 63;  
            recvind[0][0] = 41; recvind[0][1] = 42; recvind[0][2] = 43; recvind[1][0] = 51; recvind[1][1] = 52; recvind[1][2] = 53; 
            break;
        case 7:
            nsendwho = 1; nrecvwho = 1;
            sendwho = malloc(sizeof(int) * nsendwho); sendwho[0] = 4;
            recvwho = malloc(sizeof(int) * nrecvwho); recvwho[0] = 0;
            ssize = calloc(nsendwho,sizeof(int)); rsize = calloc(nrecvwho,sizeof(int));
            sendind = malloc(sizeof(int*) * nsendwho); recvind = malloc(sizeof(int*) * nrecvwho);  
            for(i=0; i < nsendwho; ++i) {sendind[i] = malloc(sizeof(int) * 3); ssize[i] = 3*embDataUnitSize;}
            for(i=0; i < nrecvwho; ++i) {recvind[i] = malloc(sizeof(int) * 3); rsize[i] = 3*embDataUnitSize;}
            sendind[0][0] = 71; sendind[0][1] = 72; sendind[0][2] = 73;  
            recvind[0][0] = 1; recvind[0][1] = 2; recvind[0][2] = 3; 
            break;
        default:
            break;

    }
    int *GtoLmap = malloc(sizeof(int) *  (100));
    for (i = 0; i < 100; ++i) GtoLmap[i] = -1;
    int tcnt = 0, nlrows = 0, j;
    for (i = 0; i < nsendwho; ++i) for (j = 0; j < 3; ++j) {if(GtoLmap[sendind[i][j]] == -1) GtoLmap[sendind[i][j]] = tcnt++;}
    for (i = 0; i < nrecvwho; ++i) for (j = 0; j < 3; ++j) {if(GtoLmap[recvind[i][j]] == -1) GtoLmap[recvind[i][j]] = tcnt++;}
    for (i = 0; i < 100; ++i) if(GtoLmap[i] > -1) nlrows++;
    assert(nlrows == tcnt);
    vptcomm_real_t *data = malloc(sizeof(*data) * nlrows * embDataUnitSize );
    for (i = 0; i < nlrows*embDataUnitSize; ++i) {
        //data[i] = ((int)((i/40)) == mypid ? ((vptcomm_real_t)mypid/10): -1) ;
        data[i] = (vptcomm_real_t)mypid/10;
    }
        
    printf( "%s\n", "done init");
    indsmap = malloc(sizeof(int) * tcnt );
    for (i = 0; i < 100; ++i) if(GtoLmap[i] > -1) indsmap[GtoLmap[i]] = i;
    for (i = 0; i < nsendwho; ++i) for (j = 0; j < 3; ++j) {sendind[i][j] = GtoLmap[sendind[i][j]];}
    for (i = 0; i < nrecvwho; ++i) for (j = 0; j < 3; ++j) {recvind[i][j] = GtoLmap[recvind[i][j]]; }
    int * HI = malloc(sizeof(*HI) * npes);
    for (i = 0; i < npes; ++i) {
       //HI[i] = npes-i; 
       HI[i] = i;
    }
    
    vptcomm_real_t *sbuff, *recvbuff, **sendp, **recvp;
    sbuff = malloc(sizeof(vptcomm_real_t) * 3 * nsendwho * embDataUnitSize);
    recvbuff = malloc(sizeof(vptcomm_real_t) * 3 * nrecvwho * embDataUnitSize);
    vptcomm_real_t *ptr = sbuff;
    for (i = 0; i < nsendwho; ++i){
        for (j = 0; j < ssize[i]/embDataUnitSize; ++j) {
            memcpy(ptr, &data[sendind[i][j] * embDataUnitSize], sizeof(*data)*embDataUnitSize);
            ptr+= embDataUnitSize;
        }
    }
    sendp = malloc(sizeof(*sendp) * nsendwho); 
    recvp = malloc(sizeof(*recvp) * nrecvwho); 
    int *sendLocs, *recvLocs;
    sendLocs = malloc(sizeof(*sendLocs) * (nsendwho+1)); recvLocs = malloc(sizeof(*recvLocs) * (nrecvwho+1));
    memcpy(sendLocs+1, ssize, nsendwho); memcpy(recvLocs+1, rsize, nrecvwho); sendLocs[0]=0; recvLocs[0]=0;
    for (i = 1; i <= nsendwho; ++i) sendLocs[i]+= sendLocs[i-1]; for (i = 1; i <= nrecvwho; ++i) recvLocs[i]+= recvLocs[i-1]; 
    for (i = 0; i < nsendwho; ++i) sendp[i] = &sbuff[sendLocs[i]]; for (i = 0; i < nrecvwho; ++i)  recvp[i] = &recvbuff[recvLocs[i]]; 
    printf( "%s\n", "done setup");


    //HI[0] = 2; HI[1] = 7; HI[2] = 4; HI[3] = 5; HI[4] = 0; HI[5] = 6; HI[6] = 1; HI[7] = 3;
    int dsizes[3]={2,2,2};
    //STFW_init_ExpRed_instance(0, 2, dsizes, nsendwho, sendwho, nrecvwho, recvwho, ssend, sendp, srecv,recvp, ;
    STFW_init_instance(0, 3, dsizes, nsendwho, sendwho, nrecvwho, recvwho, ssize, sendp, rsize, recvp, HI);
    
    printf( "%s\n", "done stfw init");

    MPI_Barrier(MPI_COMM_WORLD);

    STFW_Comm(0);
    MPI_Barrier(MPI_COMM_WORLD);
    printf( "STFW_comm done\n");


    ptr = recvbuff;
    for (i = 0; i < nrecvwho; ++i){
        for (j = 0; j < rsize[i]/embDataUnitSize; ++j) {
            memcpy(&data[recvind[i][j]*embDataUnitSize], ptr, sizeof(*data)*embDataUnitSize);
            ptr+= embDataUnitSize;
        }
    }

    for (i = 0; i < nlrows; ++i) {
        printf( "[%d] ", i);
        for (j = 0; j < embDataUnitSize; ++j) {
            printf( "%0.2f ", data[(i*embDataUnitSize)+j]);
        }
        printf( "\n");
    }
    free(data);
    MPI_Barrier(MPI_COMM_WORLD);
    printf( "%s\n", "done comm");
    STFW_finalize();
    printf( "%s\n", "done STFW_finalize");
    MPI_Finalize();
    return 0;
}

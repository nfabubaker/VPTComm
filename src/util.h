#ifndef TP_UTIL_H
#define TP_UTIL_H

#include <time.h>
#include <stdio.h>
#include "vptcomm.h"

#ifdef NA_DBG

extern FILE *dbgfp;
extern FILE *sdbgfp;
extern char dbg_fname[1024];
extern char dbg_fn[1024];
extern void na_log(FILE *fp, const char* format, ...);

#endif

typedef struct _tmr_t{                                                                                                                                                                    
    struct timespec ts_beg;
    struct timespec ts_end;
    double elapsed;
} tmr_t; 

    void setintzero(vptcomm_idx_t *arr, size_t size);
    void setreal_tzero(vptcomm_real_t *arr, size_t size);
    void setrealzero(vptcomm_real_t *arr, int size);
    long get_wc_time ( void );
    vptcomm_idx_t convert(vptcomm_idx_t value);
    void substring(char *text, char out[1024]);
    void substring_b(char *dst, char *src);
    void stop_timer(tmr_t *t);
    void start_timer(tmr_t *t);
    void init_timer (tmr_t *t);

#endif

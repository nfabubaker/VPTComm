#include "util.h"
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdarg.h>
#include <mpi.h>
#include <assert.h>
/*void *myMalloc(size_t size)
  {
  void *p = malloc(size);
  if(p == NULL)
  {
  printf("malloc couldn't allocate %d byte memory\n", size);
  exit(1);
  }
  else
  return p;

  }*/

void na_log(FILE *fp, const char *format, ...) {
    va_list args;
    va_start(args, format);
    vfprintf(fp, format, args);
    fflush(fp);
    va_end(args);
}
void *myMalloc(size_t size){
    return aligned_alloc(8*64, size);
}

void setintzero(vptcomm_idx_t *arr, size_t size)
{
    size_t i;
    for(i = 0; i < size; i++)
        arr[i] = 0;
}

/* NA: Just testing things for setreal_tzero, source: https://bytes.com/topic/c/answers/222353-safe-zero-float-array-memset */

void setreal_tzero(vptcomm_real_t *arr, size_t size){
    vptcomm_real_t *first = arr, *last = arr+size;
    for(; first!=last; ++first)
        *first = 0.0;
}
void setrealzero(vptcomm_real_t *arr, int size){
    vptcomm_real_t *first = arr, *last = arr+size;
    for(; first!=last; ++first)
        *first = 0.0;
}

void start_timer(tmr_t *t){                                                                                                                                                               
    clock_gettime(CLOCK_MONOTONIC, &t->ts_beg);                                                                                                                                           
    return;                                                                                                                                                                               
}                                                                                                                                                                                         
void stop_timer(tmr_t *t){                                                                                                                                                                
    clock_gettime(CLOCK_MONOTONIC, &t->ts_end);                                                                                                                                           
    t->elapsed += (1000000000.0 * (double) (t->ts_end.tv_sec - t->ts_beg.tv_sec) + (double) (t->ts_end.tv_nsec - t->ts_beg.tv_nsec));                                                       
    return;                                                                                                                                                                               
}                


/* long get_wc_time ( void )
 * {
 *     static struct timeval twclk ;
 *     gettimeofday(&twclk, NULL) ;
 *     return(twclk.tv_sec*1000000 + twclk.tv_usec) ;
 * }
 */

vptcomm_idx_t convert(vptcomm_idx_t value)
{
    char arr[4];
    arr[0] = (char)(value >> 24);
    arr[1] = (char)(value >> 16);
    arr[2] = (char)(value >> 8);
    arr[3] = (char)(value);

    return arr[3] << 24 | (arr[2] & 0xFF) << 16 | (arr[1] & 0xFF) << 8 | (arr[0] & 0xFF);

}

void substring(char *text, char out[1024])
{
    char *ptr = text;
    char *prevptr = NULL;

    while( (ptr = strstr(ptr,"/")))
    {
        prevptr = ptr++;
    }
    prevptr++;
    vptcomm_idx_t sl = strlen(prevptr);
    strncpy(out, prevptr, sl);
    out[sl] = '\0';

}


void substring_b(char *dst, char *src){
    char *ptr = src;
    char *prevptr = NULL;

    while( (ptr = strstr(ptr, "_"))){
        prevptr = ptr++;
    }

    //prevptr;
    vptcomm_idx_t sl = strlen(src) - strlen(prevptr);
    strncpy(dst, src, sl);
    dst[sl] = '\0';
}




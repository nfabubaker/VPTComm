//include "dbg_print.h"
#include "vptcomm.h"
#include "vpt.h"
#include "stfw.h"
#include <assert.h>
#include <mpi.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
/* #define DEBUG */

void init_cm_p(stfw_comm *cm){
    cm->st_np = 0; cm->fw_np = 0;
    cm->fw_buf = NULL; cm->fw_procs = NULL; cm->fw_locs = NULL;
    cm->fw_offset = NULL; cm->fw_offset_d = NULL; cm->fw_imap = NULL;
    cm->fw_imap_d = NULL; cm->st_buf = NULL; cm->st_procs = NULL;
    cm->st_locs = NULL; cm->st_reqs = NULL; cm->st_stts = NULL;
    cm->st_reqs_d = NULL; cm->st_stts_d = NULL;

    cm->pck_cnt = NULL; cm->pck_sizes = NULL; cm->pck_cp = NULL;

#ifdef D_SPEC_STFW
    cm->pck_fw_idx = NULL; cm->pck_src = NULL; cm->pck_dest = NULL;
#endif
    cm->unit_RepCnt = NULL;
    cm->unit_ptr = NULL;
}



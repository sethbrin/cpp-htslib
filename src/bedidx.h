#ifndef BEDIDX_H_
#define BEDIDX_H_

#ifdef _WIN32
#define drand48() ((double)rand() / RAND_MAX)
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int n, m;
    uint64_t *a;
    int *idx;
} bed_reglist_t;

#include "htslib/khash.h"
KHASH_MAP_INIT_STR(reg, bed_reglist_t)

#define LIDX_SHIFT 13

typedef kh_reg_t reghash_t;

void bed_destroy(void *_h);

int *bed_index_core(int n, uint64_t *a, int *n_idx);
void bed_index(void *_h);
int bed_overlap_core(const bed_reglist_t *p, int beg, int end);
int bed_overlap(const void *_h, const char *chr, int beg, int end);
void *bed_read(const char *fn);



#ifdef __cplusplus
}
#endif

#endif

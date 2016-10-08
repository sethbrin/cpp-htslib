#ifndef BAM2BCF_INDEL_H
#define BAM2BCF_INDEL_H

#ifdef __cplusplus
extern "C" {
#endif

#include "htslib/hts.h"
#include "htslib/sam.h"
#include "bam2bcf.h"

void *bcf_call_add_rg(void *_hash, const char *hdtext, const char *list);
void bcf_call_del_rghash(void *_hash);

int bcf_call_gap_prep(int n, int *n_plp, bam_pileup1_t **plp, int pos, bcf_callaux_t *bca, const char *ref, const void *rghash);

#ifdef __cplusplus
}
#endif

#endif

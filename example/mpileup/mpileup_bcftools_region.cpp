/*  mpileup_bcftools_region.cpp -- combine mpileup and bcftools mcall commands.

    Copyright (C) 2008-2015 Genome Research Ltd.
    Portions copyright (C) 2009-2012 Broad Institute.

    Author: Heng Li <lh3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <config.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>
#include <errno.h>
#include <sys/stat.h>
#include <getopt.h>
#include <htslib/faidx.h>
#include <htslib/kstring.h>
#include <htslib/khash_str2int.h>
#include <htslib/sam.h>
#include "bedidx.h"
#include "bam2bcf_indel.h"
#include "sam_utils.h"
#include "sam_opts.h"
#include "kthread.h"
#include "call.h"
#include "ploidy.h"

#include <easehts/utils.h>
#include <easehts/vcf.h>
#include <vector>
#include <thread>
#include <mutex>
#include <list>

static inline int printw(int c, FILE *fp)
{
    char buf[16];
    int l, x;
    if (c == 0) return fputc('0', fp);
    for (l = 0, x = c < 0? -c : c; x > 0; x /= 10) buf[l++] = x%10 + '0';
    if (c < 0) buf[l++] = '-';
    buf[l] = 0;
    for (x = 0; x < l/2; ++x) {
        int y = buf[x]; buf[x] = buf[l-1-x]; buf[l-1-x] = y;
    }
    fputs(buf, fp);
    return 0;
}

static inline void pileup_seq(FILE *fp, const bam_pileup1_t *p, int pos, int ref_len, const char *ref)
{
    int j;
    if (p->is_head) {
        putc('^', fp);
        putc(p->b->core.qual > 93? 126 : p->b->core.qual + 33, fp);
    }
    if (!p->is_del) {
        int c = p->qpos < p->b->core.l_qseq
            ? seq_nt16_str[bam_seqi(bam_get_seq(p->b), p->qpos)]
            : 'N';
        if (ref) {
            int rb = pos < ref_len? ref[pos] : 'N';
            if (c == '=' || seq_nt16_table[c] == seq_nt16_table[rb]) c = bam_is_rev(p->b)? ',' : '.';
            else c = bam_is_rev(p->b)? tolower(c) : toupper(c);
        } else {
            if (c == '=') c = bam_is_rev(p->b)? ',' : '.';
            else c = bam_is_rev(p->b)? tolower(c) : toupper(c);
        }
        putc(c, fp);
    } else putc(p->is_refskip? (bam_is_rev(p->b)? '<' : '>') : '*', fp);
    if (p->indel > 0) {
        putc('+', fp); printw(p->indel, fp);
        for (j = 1; j <= p->indel; ++j) {
            int c = seq_nt16_str[bam_seqi(bam_get_seq(p->b), p->qpos + j)];
            putc(bam_is_rev(p->b)? tolower(c) : toupper(c), fp);
        }
    } else if (p->indel < 0) {
        printw(p->indel, fp);
        for (j = 1; j <= -p->indel; ++j) {
            int c = (ref && (int)pos+j < ref_len)? ref[pos+j] : 'N';
            putc(bam_is_rev(p->b)? tolower(c) : toupper(c), fp);
        }
    }
    if (p->is_tail) putc('$', fp);
}

#include <assert.h>
#include "bam2bcf.h"
#include "sample.h"

#define MPLP_BCF        1
#define MPLP_VCF        (1<<1)
#define MPLP_NO_COMP    (1<<2)
#define MPLP_NO_ORPHAN  (1<<3)
#define MPLP_REALN      (1<<4)
#define MPLP_NO_INDEL   (1<<5)
#define MPLP_REDO_BAQ   (1<<6)
#define MPLP_ILLUMINA13 (1<<7)
#define MPLP_IGNORE_RG  (1<<8)
#define MPLP_PRINT_POS  (1<<9)
#define MPLP_PRINT_MAPQ (1<<10)
#define MPLP_PER_SAMPLE (1<<11)
#define MPLP_SMART_OVERLAPS (1<<12)

void *bed_read(const char *fn);
void bed_destroy(void *_h);
int bed_overlap(const void *_h, const char *chr, int beg, int end);

typedef struct {
  char *ref;
  int ref_id;
  int ref_len;
} ref_t;

typedef struct {
  ref_t* refs;
} mplp_ref_t;

void ref_init(ref_t& ref) {
  ref.ref = NULL;
  ref.ref_id = -1;
  ref.ref_len = 0;
}

void ref_destory(ref_t& ref) {
  free(ref.ref);
  ref.ref = NULL;
}

void mplp_ref_init(mplp_ref_t* ref, int n) {
  assert(n>0);
  ref->refs = (ref_t*)malloc(n * sizeof(ref_t));
  for (int i = 0; i < n; i++) {
    ref_init(ref->refs[i]);
  }
}

void mplp_ref_destory(mplp_ref_t* ref, int n) {
  for (int i = 0; i < n; i++) {
    ref_destory(ref->refs[i]);
  }
  free(ref->refs);
}

typedef struct {
    int min_mq, flag, min_baseQ, capQ_thres, max_depth, max_indel_depth, fmt_flag, all;
    int rflag_require, rflag_filter;
    int openQ, extQ, tandemQ, min_support; // for indels
    double min_frac; // for indels
    char *reg, *pl_list, *fai_fname, *output_fname;
    void *bed, *rghash;
    int argc;
    char **argv;
    sam_global_args ga;
    int nthreads;
    int max_lreads;
    int chunk_size;
    faidx_t *fai;
    mplp_ref_t ref;
    pthread_mutex_t mutex;

} mplp_conf_t;

typedef struct {
    samFile *fp;
    bam_hdr_t *h;
    const mplp_conf_t *conf;
} mplp_aux_t;

typedef struct {
    int n;
    int *n_plp, *m_plp;
    bam_pileup1_t **plp;
} mplp_pileup_t;

typedef struct {
  mplp_conf_t* conf;
  mplp_aux_t** data;
  bam_hdr_t* h;
  int n;
  char** fn;
  bcf_hdr_t* bcf_hdr;
  bam_sample_t* sm;
  void* rghash;
  int tid_idx;
  int tid;
  int beg;
  int end;

  call_t call;
  ncic::easehts::VCFWriter* writer;

} ktp_aux_t;

// buffer record
typedef struct {
  mplp_aux_t* data;
  void* worker_data;
  hts_itr_t *iter;
  hts_idx_t *idx;
  samFile *fp;
  bam1_t* b;
  bool b_set;
} mplp_record_t;


// store the data
typedef struct {
  int tid;
  int beg;
  int end;
  mplp_record_t** record_buf_vec;

  ktp_aux_t* aux;
#ifdef MPILEUP_STATS
  int total_reads;
#endif
} ktp_worker_t;

typedef struct {
  // Each element of worker_data_vec is use by a thread
  std::vector<ktp_worker_t*> worker_data_vec;
} ktp_workers_t;

// double buffer avoid lock
// NOTE no three reference process in a time, otherwise thread-safety
static int mplp_get_ref(ktp_aux_t *shared, int tid,  char **ref, int *ref_len) {
    ref_t* r = shared->conf->ref.refs + tid;

    //printf("get ref %d {%d/%p, %d/%p}\n", tid, r->ref_id[0], r->ref[0], r->ref_id[1], r->ref[1]);

    if (!shared->conf->fai) {
        *ref = NULL;
        return 0;
    }

    // Do we need to reference count this so multiple mplp_aux_t can
    // track which references are in use?
    // For now we just cache the last two. Sufficient?
    if (r->ref) {
        *ref = r->ref;
        *ref_len = r->ref_len;
        return 1;
    }

    pthread_mutex_lock(&shared->conf->mutex);
    if (r->ref) {
      *ref = r->ref;
      *ref_len = r->ref_len;
      pthread_mutex_unlock(&shared->conf->mutex);
      return 1;
    }
    r->ref_id = tid;
    r->ref = faidx_fetch_seq(shared->conf->fai,
                             shared->h->target_name[r->ref_id],
                             0,
                             INT_MAX,
                             &r->ref_len);

    if (!r->ref) {
        r->ref = NULL;
        r->ref_id = -1;
        r->ref_len = 0;
        *ref = NULL;
        pthread_mutex_unlock(&shared->conf->mutex);
        return 0;
    }

    *ref = r->ref;
    *ref_len = r->ref_len;
    pthread_mutex_unlock(&shared->conf->mutex);
    return 1;
}

static void
print_empty_pileup(FILE *fp, const mplp_conf_t *conf, const char *tname,
                   int pos, int n, const char *ref, int ref_len)
{
    int i;
    fprintf(fp, "%s\t%d\t%c", tname, pos+1, (ref && pos < ref_len)? ref[pos] : 'N');
    for (i = 0; i < n; ++i) {
        fputs("\t0\t*\t*", fp);
        if (conf->flag & MPLP_PRINT_MAPQ) fputs("\t*", fp);
        if (conf->flag & MPLP_PRINT_POS) fputs("\t*", fp);
    }
    putc('\n', fp);
}

static int mplp_func(void *data, bam1_t *b)
{
  mplp_record_t* record_buf = (mplp_record_t*)data;
  mplp_aux_t *ma = (mplp_aux_t*)record_buf->data;
  ktp_worker_t* worker_data = (ktp_worker_t*)record_buf->worker_data;
  int ret = -1, skip = 0, ref_len;
  char *ref;
  do {
    int has_ref;
    if (record_buf->b_set) {
      bam_copy1(b, record_buf->b);
      record_buf->b_set = false;
      ret = 1;
    } else {
      ret = sam_itr_next(record_buf->fp, record_buf->iter, b);
      if (ret < 0) break;
    }
#ifdef MPILEUP_STATS
    worker_data->total_reads ++;
#endif
    // The 'B' cigar operation is not part of the specification, considering as obsolete.
    //  bam_remove_B(b);
    if (b->core.tid < 0 || (b->core.flag&BAM_FUNMAP)) { // exclude unmapped reads
      skip = 1;
      continue;
    }
    if (ma->conf->rflag_require && !(ma->conf->rflag_require&b->core.flag)) { skip = 1; continue; }
    if (ma->conf->rflag_filter && ma->conf->rflag_filter&b->core.flag) { skip = 1; continue; }
    if (ma->conf->bed && ma->conf->all == 0) { // test overlap
      skip = !bed_overlap(ma->conf->bed, ma->h->target_name[b->core.tid], b->core.pos, bam_endpos(b));
      if (skip) continue;
    }
    if (ma->conf->rghash) { // exclude read groups
      uint8_t *rg = bam_aux_get(b, "RG");
      skip = (rg && khash_str2int_get(ma->conf->rghash, (const char*)(rg+1), NULL)==0);
      if (skip) continue;
    }
    if (ma->conf->flag & MPLP_ILLUMINA13) {
      int i;
      uint8_t *qual = bam_get_qual(b);
      for (i = 0; i < b->core.l_qseq; ++i)
        qual[i] = qual[i] > 31? qual[i] - 31 : 0;
    }

    if (ma->conf->fai && b->core.tid >= 0) {
      has_ref = mplp_get_ref(worker_data->aux, b->core.tid, &ref, &ref_len);
      if (has_ref && ref_len <= b->core.pos) { // exclude reads outside of the reference sequence
        fprintf(stderr,"[%s] Skipping because %d is outside of %d [ref:%d]\n",
                __func__, b->core.pos, ref_len, b->core.tid);
        skip = 1;
        continue;
      }
    } else {
      has_ref = 0;
    }

    skip = 0;
    if (has_ref && (ma->conf->flag&MPLP_REALN)) sam_prob_realn(b, ref, ref_len, (ma->conf->flag & MPLP_REDO_BAQ)? 7 : 3);
    if (has_ref && ma->conf->capQ_thres > 10) {
      int q = sam_cap_mapq(b, ref, ref_len, ma->conf->capQ_thres);
      if (q < 0) skip = 1;
      else if (b->core.qual > q) b->core.qual = q;
    }
    if (b->core.qual < ma->conf->min_mq) skip = 1;
    else if ((ma->conf->flag&MPLP_NO_ORPHAN) && (b->core.flag&BAM_FPAIRED) && !(b->core.flag&BAM_FPROPER_PAIR)) skip = 1;
  } while (skip);
  return ret;
}

static void group_smpl(mplp_pileup_t *m, bam_sample_t *sm, kstring_t *buf,
                       int n, char *const*fn, int *n_plp, const bam_pileup1_t **plp, int ignore_rg)
{
    int i, j;
    memset(m->n_plp, 0, m->n * sizeof(int));
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n_plp[i]; ++j) {
            const bam_pileup1_t *p = plp[i] + j;
            uint8_t *q;
            int id = -1;
            q = ignore_rg? NULL : bam_aux_get(p->b, "RG");
            if (q) id = bam_smpl_rg2smid(sm, fn[i], (char*)q+1, buf);
            //if (id < 0) id = bam_smpl_rg2smid(sm, fn[i], 0, buf);
            if (id < 0) {
              if (n == 1) {
                id = 0;
              } else {
                id = bam_smpl_rg2smid(sm, fn[i], 0, buf);
              }
            }
            if (id < 0 || id >= m->n) {
                assert(q); // otherwise a bug
                fprintf(stderr, "[%s] Read group %s used in file %s but absent from the header or an alignment missing read group.\n", __func__, (char*)q+1, fn[i]);
                exit(EXIT_FAILURE);
            }
            if (m->n_plp[id] == m->m_plp[id]) {
                m->m_plp[id] = m->m_plp[id]? m->m_plp[id]<<1 : 8;
                m->plp[id] = (bam_pileup1_t*)realloc(m->plp[id], sizeof(bam_pileup1_t) * m->m_plp[id]);
            }
            m->plp[id][m->n_plp[id]++] = *p;
        }
    }
}

/**
 * Write VCF header
 */
static void write_header(mplp_conf_t*& conf, int n, char **fn,
                         bcf_hdr_t* bcf_hdr,
                         bam_sample_t*& sm,
                         void*& rghash, bam_hdr_t*& h,
                         mplp_aux_t** data) {
  // Add sample
  // read the header of each file in the list and initialize data
  for (int i = 0; i < n; ++i) {
    bam_hdr_t *h_tmp;
    data[i] = (mplp_aux_t*)calloc(1, sizeof(mplp_aux_t));
    data[i]->fp = sam_open_format(fn[i], "rb", &conf->ga.in);
    if ( !data[i]->fp )
    {
      fprintf(stderr, "[%s] failed to open %s: %s\n", __func__, fn[i], strerror(errno));
      exit(EXIT_FAILURE);
    }
    if (hts_set_opt(data[i]->fp, CRAM_OPT_DECODE_MD, 0)) {
      fprintf(stderr, "Failed to set CRAM_OPT_DECODE_MD value\n");
      exit(EXIT_FAILURE);
    }
    if (conf->fai_fname && hts_set_fai_filename(data[i]->fp, conf->fai_fname) != 0) {
      fprintf(stderr, "[%s] failed to process %s: %s\n",
              __func__, conf->fai_fname, strerror(errno));
      exit(EXIT_FAILURE);
    }
    data[i]->conf = conf;

    h_tmp = sam_hdr_read(data[i]->fp);
    if ( !h_tmp ) {
      fprintf(stderr,"[%s] fail to read the header of %s\n", __func__, fn[i]);
      exit(EXIT_FAILURE);
    }
    bam_smpl_add(sm, fn[i], (conf->flag&MPLP_IGNORE_RG)? 0 : h_tmp->text);
    // Collect read group IDs with PL (platform) listed in pl_list (note: fragile, strstr search)
    rghash = bcf_call_add_rg(rghash, h_tmp->text, conf->pl_list);

    if (i == 0) data[i]->h = h = h_tmp; // save the header of the first file
    else {
      // FIXME: check consistency between h and h_tmp
      bam_hdr_destroy(h_tmp);

      // we store only the first file's header; it's (alleged to be)
      // compatible with the i-th file's target_name lookup needs
      data[i]->h = h;
    }
  }

  // write the VCF header
  if (conf->flag & MPLP_BCF)
  {
    kstring_t str = {0,0,NULL};

    ksprintf(&str, "##samtoolsVersion=%s+htslib-%s\n",samtools_version(),hts_version());
    bcf_hdr_append(bcf_hdr, str.s);

    str.l = 0;
    ksprintf(&str, "##samtoolsCommand=samtools mpileup");
    for (int i=1; i<conf->argc; i++) ksprintf(&str, " %s", conf->argv[i]);
    kputc('\n', &str);
    bcf_hdr_append(bcf_hdr, str.s);

    if (conf->fai_fname)
    {
      str.l = 0;
      ksprintf(&str, "##reference=file://%s\n", conf->fai_fname);
      bcf_hdr_append(bcf_hdr, str.s);
    }

    // Translate BAM @SQ tags to BCF ##contig tags
    // todo: use/write new BAM header manipulation routines, fill also UR, M5
    for (int i=0; i<h->n_targets; i++)
    {
      str.l = 0;
      ksprintf(&str, "##contig=<ID=%s,length=%d>", h->target_name[i], h->target_len[i]);
      bcf_hdr_append(bcf_hdr, str.s);
    }
    free(str.s);
    bcf_hdr_append(bcf_hdr,"##ALT=<ID=*,Description=\"Represents allele(s) other than observed.\">");
    bcf_hdr_append(bcf_hdr,"##INFO=<ID=INDEL,Number=0,Type=Flag,Description=\"Indicates that the variant is an INDEL.\">");
    bcf_hdr_append(bcf_hdr,"##INFO=<ID=IDV,Number=1,Type=Integer,Description=\"Maximum number of reads supporting an indel\">");
    bcf_hdr_append(bcf_hdr,"##INFO=<ID=IMF,Number=1,Type=Float,Description=\"Maximum fraction of reads supporting an indel\">");
    bcf_hdr_append(bcf_hdr,"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw read depth\">");
    bcf_hdr_append(bcf_hdr,"##INFO=<ID=VDB,Number=1,Type=Float,Description=\"Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)\",Version=\"3\">");
    bcf_hdr_append(bcf_hdr,"##INFO=<ID=RPB,Number=1,Type=Float,Description=\"Mann-Whitney U test of Read Position Bias (bigger is better)\">");
    bcf_hdr_append(bcf_hdr,"##INFO=<ID=MQB,Number=1,Type=Float,Description=\"Mann-Whitney U test of Mapping Quality Bias (bigger is better)\">");
    bcf_hdr_append(bcf_hdr,"##INFO=<ID=BQB,Number=1,Type=Float,Description=\"Mann-Whitney U test of Base Quality Bias (bigger is better)\">");
    bcf_hdr_append(bcf_hdr,"##INFO=<ID=MQSB,Number=1,Type=Float,Description=\"Mann-Whitney U test of Mapping Quality vs Strand Bias (bigger is better)\">");
#if CDF_MWU_TESTS
    bcf_hdr_append(bcf_hdr,"##INFO=<ID=RPB2,Number=1,Type=Float,Description=\"Mann-Whitney U test of Read Position Bias [CDF] (bigger is better)\">");
    bcf_hdr_append(bcf_hdr,"##INFO=<ID=MQB2,Number=1,Type=Float,Description=\"Mann-Whitney U test of Mapping Quality Bias [CDF] (bigger is better)\">");
    bcf_hdr_append(bcf_hdr,"##INFO=<ID=BQB2,Number=1,Type=Float,Description=\"Mann-Whitney U test of Base Quality Bias [CDF] (bigger is better)\">");
    bcf_hdr_append(bcf_hdr,"##INFO=<ID=MQSB2,Number=1,Type=Float,Description=\"Mann-Whitney U test of Mapping Quality vs Strand Bias [CDF] (bigger is better)\">");
#endif
    bcf_hdr_append(bcf_hdr,"##INFO=<ID=SGB,Number=1,Type=Float,Description=\"Segregation based metric.\">");
    bcf_hdr_append(bcf_hdr,"##INFO=<ID=MQ0F,Number=1,Type=Float,Description=\"Fraction of MQ0 reads (smaller is better)\">");
    bcf_hdr_append(bcf_hdr,"##INFO=<ID=I16,Number=16,Type=Float,Description=\"Auxiliary tag used for calling, see description of bcf_callret1_t in bam2bcf.h\">");
    bcf_hdr_append(bcf_hdr,"##INFO=<ID=QS,Number=R,Type=Float,Description=\"Auxiliary tag used for calling\">");
    bcf_hdr_append(bcf_hdr,"##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"List of Phred-scaled genotype likelihoods\">");
    if ( conf->fmt_flag&B2B_FMT_DP )
      bcf_hdr_append(bcf_hdr,"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Number of high-quality bases\">");
    if ( conf->fmt_flag&B2B_FMT_DV )
      bcf_hdr_append(bcf_hdr,"##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"Number of high-quality non-reference bases\">");
    if ( conf->fmt_flag&B2B_FMT_DPR )
      bcf_hdr_append(bcf_hdr,"##FORMAT=<ID=DPR,Number=R,Type=Integer,Description=\"Number of high-quality bases observed for each allele\">");
    if ( conf->fmt_flag&B2B_INFO_DPR )
      bcf_hdr_append(bcf_hdr,"##INFO=<ID=DPR,Number=R,Type=Integer,Description=\"Number of high-quality bases observed for each allele\">");
    if ( conf->fmt_flag&B2B_FMT_DP4 )
      bcf_hdr_append(bcf_hdr,"##FORMAT=<ID=DP4,Number=4,Type=Integer,Description=\"Number of high-quality ref-fwd, ref-reverse, alt-fwd and alt-reverse bases\">");
    if ( conf->fmt_flag&B2B_FMT_SP )
      bcf_hdr_append(bcf_hdr,"##FORMAT=<ID=SP,Number=1,Type=Integer,Description=\"Phred-scaled strand bias P-value\">");
    if ( conf->fmt_flag&B2B_FMT_AD )
      bcf_hdr_append(bcf_hdr,"##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">");
    if ( conf->fmt_flag&B2B_FMT_ADF )
      bcf_hdr_append(bcf_hdr,"##FORMAT=<ID=ADF,Number=R,Type=Integer,Description=\"Allelic depths on the forward strand\">");
    if ( conf->fmt_flag&B2B_FMT_ADR )
      bcf_hdr_append(bcf_hdr,"##FORMAT=<ID=ADR,Number=R,Type=Integer,Description=\"Allelic depths on the reverse strand\">");
    if ( conf->fmt_flag&B2B_INFO_AD )
      bcf_hdr_append(bcf_hdr,"##INFO=<ID=AD,Number=R,Type=Integer,Description=\"Total allelic depths\">");
    if ( conf->fmt_flag&B2B_INFO_ADF )
      bcf_hdr_append(bcf_hdr,"##INFO=<ID=ADF,Number=R,Type=Integer,Description=\"Total allelic depths on the forward strand\">");
    if ( conf->fmt_flag&B2B_INFO_ADR )
      bcf_hdr_append(bcf_hdr,"##INFO=<ID=ADR,Number=R,Type=Integer,Description=\"Total allelic depths on the reverse strand\">");

    for (int i=0; i<sm->n; i++)
      bcf_hdr_add_sample(bcf_hdr, sm->smpl[i]);
    bcf_hdr_add_sample(bcf_hdr, NULL);
    // End of BCF header creation
    if ( bcf_hdr->dirty ) bcf_hdr_sync(bcf_hdr);

  }
  else {
    ERROR("Unsupported, current only support MPLP_BCF or MPLP_VCF");
  }

}

static void ktp_worker_destory(ktp_worker_t* worker_data, int n) {
  for (int i = 0; i < n; i++) {
    mplp_record_t* record_buf = worker_data->record_buf_vec[i];
    if (record_buf->iter) hts_itr_destroy(record_buf->iter);
    hts_idx_destroy(record_buf->idx);
    bam_destroy1(record_buf->b);
    sam_close(record_buf->fp);
    delete record_buf;
  }
  delete []worker_data->record_buf_vec;
  delete worker_data;
}

/**
 * Decide if the current tid has reads
 */
static bool decide_next_tid(ktp_aux_t* shared, ktp_worker_t* worker_data,
                            const std::vector<int>& sorted_tids) {
  int min_val = INT_MAX;
  bam1_t* record = bam_init1();
  bool res = false;
  do {
    shared->tid_idx ++;
    if (shared->tid_idx >= sorted_tids.size()) {
      break;
    }
    shared->tid = sorted_tids[shared->tid_idx];
    for (int i = 0; i < shared->n; i++) {
      mplp_record_t* record_buf = worker_data->record_buf_vec[i];
      // NOTE the end should +1
      if ( (record_buf->iter=sam_itr_queryi(
                  record_buf->idx, shared->tid,
                  0,
                  shared->h->target_len[shared->tid])) == 0) {
        fprintf(stderr, "[E::%s] fail to parse reference '%s' with %s\n",
                __func__, shared->h->target_name[shared->tid], shared->fn[i]);
        exit(EXIT_FAILURE);
      }

      int ret = sam_itr_next(record_buf->fp, record_buf->iter, record);
      if (ret >= 0) {
        int cur = record->core.pos;
        if (cur < min_val) {
          min_val = cur;
        }
      }
      hts_itr_destroy(record_buf->iter);
      record_buf->iter = nullptr;
    }

    if (min_val == INT_MAX) {
      continue;
    }

    shared->beg = min_val;
    res = true;
    break;
  } while (true);
  bam_destroy1(record);
  return res;
}

/**
 * Read record from bam files
 * Store the read record in the region tid:beg-end
 *
 */
static bool read_records(ktp_aux_t* shared, ktp_worker_t* worker_data,
                         const std::vector<int>& sorted_tids,
                         std::mutex& mtx) {
  bam_hdr_t* h = shared->h;

  do {
    int beg, end, tid;
    {
      std::lock_guard<std::mutex> lock(mtx);
      if (shared->tid_idx >= static_cast<int>(sorted_tids.size())) {
        return false;
      }
      beg = shared->end + 1;
      if (shared->tid == -1 ||
          beg >= h->target_len[shared->tid]) {
        do {
          // decide if the current tid has reads
          // add set the beg
          if (decide_next_tid(shared, worker_data, sorted_tids)) {
            beg = shared->beg;
            break;
          } else {
            return false;
          }
        } while (true);
      }
      tid = shared->tid;
      end = beg + shared->conf->chunk_size - 1;
      shared->beg = beg;
      shared->end = end;
      worker_data->beg = beg;
      worker_data->end = end;
      worker_data->tid = tid;
#ifdef MPILEUP_STATS
      worker_data->total_reads = 0;
#endif
    }

    bool flag = false;
    for (int i = 0; i < shared->n; i++) {
      mplp_record_t* record_buf = worker_data->record_buf_vec[i];
      // NOTE the end should +1
      if (record_buf->iter) hts_itr_destroy(record_buf->iter);
      if ( (record_buf->iter=sam_itr_queryi(
                  record_buf->idx, tid,
                  beg - shared->conf->max_lreads,
                  end + shared->conf->max_lreads+1)) == 0) {
        fprintf(stderr, "[E::%s] fail to parse region '%d:%d-%d' with %s\n",
                __func__, tid, beg, end, shared->fn[i]);
        exit(EXIT_FAILURE);
      }
      int ret = sam_itr_next(record_buf->fp, record_buf->iter, record_buf->b);

      if (ret >= 0) {
        record_buf->b_set = true;
        flag = true;
      } else {
        hts_itr_destroy(record_buf->iter);
        record_buf->iter = nullptr;
      }

      record_buf->worker_data = worker_data;
      record_buf->data = worker_data->aux->data[i];
      worker_data->record_buf_vec[i] = record_buf;
    }

    if (flag) {
      return true;
    }
  } while (true);
  return false;
}

static void call_deep_copy(call_t* dst, call_t* src) {
  memcpy(dst, src, sizeof(call_t));
  dst->nqsum = 5;
  dst->qsum  = (float*) malloc(sizeof(float)*dst->nqsum); // will be expanded later if ncessary
  dst->nals_map = 5;
  dst->als_map  = (int*) malloc(sizeof(int)*dst->nals_map);
  dst->npl_map  = 5*(5+1)/2;     // will be expanded later if necessary
  dst->pl_map   = (int*) malloc(sizeof(int)*dst->npl_map);
  dst->gts  = (int32_t*) calloc(bcf_hdr_nsamples(src->hdr)*2,sizeof(int32_t));   // assuming at most diploid everywhere
}

static bool mcall_process(call_t* call, bcf1_t* bcf_rec) {
  bcf_unpack(bcf_rec, BCF_UN_STR);

  // Which allele is symbolic? All SNPs should have it, but not indels
  call->unseen = 0;
  for (int i=1; i<bcf_rec->n_allele; i++)
  {
    if ( bcf_rec->d.allele[i][0]=='X' ) { call->unseen = i; break; }  // old X
    if ( bcf_rec->d.allele[i][0]=='<' )
    {
      if ( bcf_rec->d.allele[i][1]=='X' && bcf_rec->d.allele[i][2]=='>' ) { call->unseen = i; break; } // old <X>
      if ( bcf_rec->d.allele[i][1]=='*' && bcf_rec->d.allele[i][2]=='>' ) { call->unseen = i; break; } // new <*>
    }
  }
  int is_ref = (bcf_rec->n_allele==1 || (bcf_rec->n_allele==2 && call->unseen>0)) ? 1 : 0;

  if ( is_ref && call->flag&CALL_VARONLY )
    return false;

  bcf_unpack(bcf_rec, BCF_UN_ALL);

  // Calling modes which output VCFs
  int ret;
  ret = mcall(call, bcf_rec);
  if ( ret==-1 ) error("Something is wrong\n");
  if ( (call->flag & CALL_VARONLY) && ret==0) return false;     // not a variant
  if (!bcf_rec) return false;
  return true;
}

/**
 * The core process function, process a region
 *
 * [beg, end]
 *
 * @return -1 0 1
 * -1 means some error occured
 *  0 means end
 *  1 means the iter do not cover any reads
 */
static void worker(ktp_aux_t* shared,
                   ktp_worker_t* worker_data,
                   errmod_t* em,
                   int* ret) {
#ifdef MPILEUP_STATS
  int region_total_depth = 0;
#endif

  fprintf(stderr, "Process %s:%d-%d\n",
          shared->h->target_name[((ktp_worker_t*)worker_data)->tid],
          ((ktp_worker_t*)worker_data)->beg,
          ((ktp_worker_t*)worker_data)->end);
  mplp_conf_t* conf = shared->conf;
  int n = shared->n;
  char** fn = shared->fn;
  bcf_hdr_t* bcf_hdr = shared->bcf_hdr;
  //bam_sample_t* sm = bam_smpl_init();
  bam_sample_t* sm = shared->sm;
  void* rghash = shared->rghash;
  bam_hdr_t* h = shared->h;

  // init sm
  //for (int i = 0; i < n; ++i) {
  //  bam_smpl_add(sm, fn[i],
  //               (shared->conf->flag&MPLP_IGNORE_RG)? 0 : h->text);
  //}

  bam_mplp_t iter;
  const bam_pileup1_t **plp;
  int* n_plp;
  char *ref;
  int i, pos, ref_len, max_depth, max_indel_depth;
  bcf_callaux_t *bca = NULL;
  bcf_callret1_t *bcr = NULL;
  bcf_call_t bc;
  call_t call;
  call_deep_copy(&call, &shared->call);

  int tid = worker_data->tid;
  int beg = worker_data->beg;
  int end = worker_data->end;


  mplp_pileup_t gplp;
  kstring_t buf;
  memset(&buf, 0, sizeof(kstring_t));
  memset(&bc, 0, sizeof(bcf_call_t));
  memset(&gplp, 0, sizeof(mplp_pileup_t));
  plp = (const bam_pileup1_t**)calloc(n, sizeof(bam_pileup1_t*));
  n_plp = (int*)calloc(n, sizeof(int));

  // Initialise the calling algorithm
  bca = bcf_call_init2(-1., conf->min_baseQ, em);
  //bca = bcf_call_init(-1., conf->min_baseQ);
  bcr = (bcf_callret1_t*)calloc(sm->n, sizeof(bcf_callret1_t));
  bca->rghash = rghash;
  bca->openQ = conf->openQ, bca->extQ = conf->extQ, bca->tandemQ = conf->tandemQ;
  bca->min_frac = conf->min_frac;
  bca->min_support = conf->min_support;
  bca->per_sample_flt = conf->flag & MPLP_PER_SAMPLE;

  bc.bcf_hdr = bcf_hdr;
  bc.n = sm->n;
  bc.PL = (int32_t*)malloc(15 * sm->n * sizeof(*bc.PL));
  if (conf->fmt_flag)
  {
    assert( sizeof(float)==sizeof(int32_t) );
    bc.DP4 = (int32_t*)malloc(sm->n * sizeof(int32_t) * 4);
    bc.fmt_arr = (uint8_t*)malloc(sm->n * sizeof(float)); // all fmt_flag fields
    if ( conf->fmt_flag&(B2B_INFO_DPR|B2B_FMT_DPR|B2B_INFO_AD|B2B_INFO_ADF|B2B_INFO_ADR|B2B_FMT_AD|B2B_FMT_ADF|B2B_FMT_ADR) )
    {
      // first B2B_MAX_ALLELES fields for total numbers, the rest per-sample
      bc.ADR = (int32_t*) malloc((sm->n+1)*B2B_MAX_ALLELES*sizeof(int32_t));
      bc.ADF = (int32_t*) malloc((sm->n+1)*B2B_MAX_ALLELES*sizeof(int32_t));
      for (int i=0; i<sm->n; i++)
      {
        bcr[i].ADR = bc.ADR + (i+1)*B2B_MAX_ALLELES;
        bcr[i].ADF = bc.ADF + (i+1)*B2B_MAX_ALLELES;
      }
    }
  }


  // allocate data storage proportionate to number of samples being studied sm->n
  gplp.n = sm->n;
  gplp.n_plp = (int*)calloc(sm->n, sizeof(int));
  gplp.m_plp = (int*)calloc(sm->n, sizeof(int));
  gplp.plp = (bam_pileup1_t**)calloc(sm->n, sizeof(bam_pileup1_t*));

  //fprintf(stderr, "[%s] %d samples in %d input files\n", __func__, sm->n, n);
  // init pileup
  iter = bam_mplp_init(n, mplp_func, (void**)worker_data->record_buf_vec);
  if ( conf->flag & MPLP_SMART_OVERLAPS ) bam_mplp_init_overlaps(iter);
  max_depth = conf->max_depth;
  if (max_depth * sm->n > 1<<20)
    fprintf(stderr, "(%s) Max depth is above 1M. Potential memory hog!\n", __func__);
  if (max_depth * sm->n < 8000) {
    max_depth = 8000 / sm->n;
    //fprintf(stderr, "<%s> Set max per-file depth to %d\n", __func__, max_depth);
  }
  max_indel_depth = conf->max_indel_depth * sm->n;
  bam_mplp_set_maxcnt(iter, max_depth);
  *ret = 0;
  // if the iter cover any reads
  bool cover = false;

  bcf1_t *bcf_rec = bcf_init1();
  // begin pileup
  while ( (*ret=bam_mplp_auto(iter, &tid, &pos, n_plp, plp)) > 0) {
    if (pos < beg || pos > end) continue; // out of the region requested

    cover = true;
    mplp_get_ref(shared, tid, &ref, &ref_len);
    //printf("tid=%d len=%d ref=%p/%s\n", tid, ref_len, ref, ref);
    int total_depth, _ref0, ref16;
    if (conf->bed && tid >= 0 && !bed_overlap(conf->bed, h->target_name[tid], pos, pos+1)) continue;
    for (i = total_depth = 0; i < n; ++i) total_depth += n_plp[i];
#ifdef MPILEUP_STATS
    region_total_depth += total_depth;
#endif
    group_smpl(&gplp, sm, &buf, n, fn, n_plp, plp, conf->flag & MPLP_IGNORE_RG);
    _ref0 = (ref && pos < ref_len)? ref[pos] : 'N';
    ref16 = seq_nt16_table[_ref0];
    bcf_callaux_clean(bca, &bc);
    for (i = 0; i < gplp.n; ++i)
      bcf_call_glfgen(gplp.n_plp[i], gplp.plp[i], ref16, bca, bcr + i);
    bc.tid = tid; bc.pos = pos;
    bcf_call_combine(gplp.n, bcr, bca, ref16, &bc);
    bcf_clear1(bcf_rec);
    bcf_call2bcf(&bc, bcf_rec, bcr, conf->fmt_flag, 0, 0);

    //bcf_write1(bcf_fp, bcf_hdr, bcf_rec);
    if (mcall_process(&call, bcf_rec)) {
      shared->writer->Add(std::unique_ptr<ncic::easehts::VariantContext>(
              new ncic::easehts::VariantContext(bcf_rec)));
    }
    // call indels; todo: subsampling with total_depth>max_indel_depth instead of ignoring?
    if (!(conf->flag&MPLP_NO_INDEL) && total_depth < max_indel_depth && bcf_call_gap_prep(gplp.n, gplp.n_plp, gplp.plp, pos, bca, ref, rghash) >= 0)
    {
      bcf_callaux_clean(bca, &bc);
      for (i = 0; i < gplp.n; ++i)
        bcf_call_glfgen(gplp.n_plp[i], gplp.plp[i], -1, bca, bcr + i);
      if (bcf_call_combine(gplp.n, bcr, bca, -1, &bc) >= 0) {
        bcf_clear1(bcf_rec);
        bcf_call2bcf(&bc, bcf_rec, bcr, conf->fmt_flag, bca, ref);
        //bcf_write1(bcf_fp, bcf_hdr, bcf_rec);
        if (mcall_process(&call, bcf_rec)) {
          shared->writer->Add(std::unique_ptr<ncic::easehts::VariantContext>(
                  new ncic::easehts::VariantContext(bcf_rec)));
        }
      }
    }
  }
  // the past will set ret to 0/-1, 0 means cover the reads
  // -1 means some error occurred
  if (!cover) {
    *ret = 1;
  }
  //bam_smpl_destroy(sm);

  for (i = 0; i < gplp.n; ++i) free(gplp.plp[i]);
  free(gplp.plp); free(gplp.n_plp); free(gplp.m_plp);
  free(buf.s);
  // clean up
  free(bc.tmp.s);
  bcf_call_destroy2(bca);
  free(bc.PL);
  free(bc.DP4);
  free(bc.ADR);
  free(bc.ADF);
  free(bc.fmt_arr);
  free(bcr);
  bcf_destroy1(bcf_rec);

  bam_mplp_destroy(iter);

  free(plp); free(n_plp);
  mcall_destroy(&call);
#ifdef MPILEUP_STATS
  fprintf(stderr, "Process done %s:%d-%d\treads:%d\tdepth:%d\n",
          shared->h->target_name[((ktp_worker_t*)worker_data)->tid],
          ((ktp_worker_t*)worker_data)->beg,
          ((ktp_worker_t*)worker_data)->end,
          ((ktp_worker_t*)worker_data)->total_reads,
          region_total_depth);
#else
  fprintf(stderr, "Process done %s:%d-%d\n",
          shared->h->target_name[((ktp_worker_t*)worker_data)->tid],
          ((ktp_worker_t*)worker_data)->beg,
          ((ktp_worker_t*)worker_data)->end);
#endif
}

// The num_reads and length of the reference
typedef struct {
  uint64_t num_reads;
  uint32_t length;
  int tid;
} reference_depth_t;

/**
 * Calc the sorted_tids which is sorted by reads of reference / length of
 * reference
 * here we just use the fomular to substitute  depth
 */
static std::vector<int> sorted_reference_by_depth(const ktp_aux_t* shared) {
  std::vector<reference_depth_t> sorted_tids;
  bam_hdr_t* h = shared->h;

  for (int tid = 0; tid < h->n_targets; tid++) {
    sorted_tids.push_back({0, h->target_len[tid], tid});

  }


  for (int i = 0; i < shared->n; i++) {
    hts_idx_t *idx = hts_idx_load(shared->fn[i], HTS_FMT_BAI);
    if (idx == NULL) {
      fprintf(stderr, "[%s] fail to load index for %s\n", __func__, shared->fn[i]);
      exit(EXIT_FAILURE);
    }

    for (int tid = 0; tid < h->n_targets; tid++) {
      uint64_t mapped, unmapped;
      hts_idx_get_stat(idx, tid, &mapped, &unmapped);

      sorted_tids[tid].num_reads += (mapped + unmapped);
    }
  }

  std::sort(sorted_tids.begin(), sorted_tids.end(),
            [](const reference_depth_t& lhs,
               const reference_depth_t& rhs)->bool {
            return lhs.num_reads/static_cast<double>(lhs.length) >
              rhs.num_reads/static_cast<double>(rhs.length);
            });
  std::vector<int> res;
  for (const auto& item : sorted_tids) {
    if (item.num_reads == 0) {
      break;
    } else {
      res.push_back(item.tid);
    }
  }
  return res;
};

static void mem_process(ktp_aux_t* shared) {
  std::vector<int> sorted_tids = sorted_reference_by_depth(shared);
  std::mutex mtx;
  std::vector<std::thread> threads;
  // 0.83 == CALL_DEFTHETA
  errmod_t* em = errmod_init(1.0 - 0.83);
  for (int i = 0; i < shared->conf->nthreads; i++) {
    threads.emplace_back([&shared, &mtx, &sorted_tids, em]() {
      ktp_worker_t* worker_data = new ktp_worker_t;
      worker_data->aux = shared;
      worker_data->record_buf_vec = new mplp_record_t*[shared->n];
      hts_idx_t* idxs[shared->n];
      for (int i = 0; i < shared->n; i++) {
        worker_data->record_buf_vec[i] = new mplp_record_t;
        mplp_record_t* record_buf = worker_data->record_buf_vec[i];

        record_buf->iter = nullptr;
        record_buf->b = bam_init1();
        record_buf->b_set = false;
        record_buf->fp = sam_open_format(shared->fn[i], "rb", &shared->conf->ga.in);
        if ( !record_buf->fp )
        {
          fprintf(stderr, "[%s] failed to open %s: %s\n", __func__, shared->fn[i], strerror(errno));
          exit(EXIT_FAILURE);
        }
        if (hts_set_opt(record_buf->fp, CRAM_OPT_DECODE_MD, 0)) {
          fprintf(stderr, "Failed to set CRAM_OPT_DECODE_MD value\n");
          exit(EXIT_FAILURE);
        }
        if (shared->conf->fai_fname && hts_set_fai_filename(record_buf->fp, shared->conf->fai_fname) != 0) {
          fprintf(stderr, "[%s] failed to process %s: %s\n",
                  __func__, shared->conf->fai_fname, strerror(errno));
          exit(EXIT_FAILURE);
        }
        hts_idx_t *idx = sam_index_load(record_buf->fp, shared->fn[i]);
        if (idx == NULL) {
          fprintf(stderr, "[%s] fail to load index for %s\n", __func__, shared->fn[i]);
          exit(EXIT_FAILURE);
        }
        record_buf->idx = idx;
      }

      while (true) {
        if (read_records(shared, worker_data, sorted_tids, mtx)) {
          int ret;
          worker(shared, worker_data, em, &ret);
        } else {
          break;
        }
      }

      ktp_worker_destory(worker_data, shared->n);
    });
  }

  for (int i = 0; i < shared->conf->nthreads; i++) {
    threads[i].join();
  }
  errmod_destroy(em);
}

static void init_aux_call(ktp_aux_t* aux) {
  memset(&aux->call, 0, sizeof(call_t));
  aux->call.prior_type = -1;
  aux->call.indel_frac = -1;
  aux->call.theta      = 1.1e-3;
  aux->call.pref       = 0.5;
  aux->call.min_perm_p = 0.01;
  aux->call.min_lrt    = 1;
  aux->call.trio_Pm_SNPs = 1 - 1e-8;
  aux->call.trio_Pm_ins  = aux->call.trio_Pm_del  = 1 - 1e-9;
  // the current used call flag is vm
  aux->call.flag |= CALL_VARONLY;
  aux->call.hdr = aux->bcf_hdr;
  //ploidy_t *ploidy = ploidy_init_string("* * * 0 0\n* * * 1 1\n* * * 2 2\n",2);
  //int nsex = ploidy_nsex(ploidy);
  //ploidy_destroy(ploidy);
  int nsamples = bcf_hdr_nsamples(aux->call.hdr);
  aux->call.ploidy = (uint8_t*)malloc(nsamples);
  for (int i = 0; i < nsamples; i++) {
    aux->call.ploidy[i] = 2;
  }

  mcall_init(&aux->call);

  //bcf_hdr_remove(aux->call.hdr, BCF_HL_INFO, "QS");
  //bcf_hdr_remove(aux->call.hdr, BCF_HL_INFO, "I16");
}

/*
 * Performs pileup
 * @param conf configuration for this pileup
 * @param n number of files specified in fn
 * @param fn filenames
 *
 * NOTE current version is only used with following conditions:
 * 1. (conf.flag & MPLP_BCF) & 1 = 1
 * 2. only perform the bam files, but no sam files
 */
static int mpileup(mplp_conf_t *conf, int n, char **fn)
{
  fprintf(stderr, "nthreads:%d -- max_lreads:%d\n", conf->nthreads, conf->max_lreads);
  extern void *bcf_call_add_rg(void *rghash, const char *hdtext, const char *list);
  extern void bcf_call_del_rghash(void *rghash);
  bam_hdr_t *h = NULL; /* header of first file in input list */
  void *rghash = NULL;
  mplp_aux_t **data;

  data = (mplp_aux_t**)calloc(n, sizeof(mplp_aux_t*));

  bam_sample_t *sm = NULL;
  sm = bam_smpl_init();

  if (n == 0) {
    fprintf(stderr,"[%s] no input file/data given\n", __func__);
    exit(EXIT_FAILURE);
  }

  const char *mode;
  if ( conf->flag & MPLP_VCF )
    mode = (conf->flag&MPLP_NO_COMP)? "wu" : "wz";   // uncompressed VCF or compressed VCF
  else
    mode = (conf->flag&MPLP_NO_COMP)? "wub" : "wb";  // uncompressed BCF or compressed BCF
  ncic::easehts::VCFWriter writer(conf->output_fname?
                                      conf->output_fname : "-",
                                      mode);

  bcf_hdr_t* bcf_hdr = writer.GetHeader()->GetRawHeader();
  write_header(conf, n, fn, bcf_hdr, sm, rghash, h, data);
  mplp_ref_init(&conf->ref, h->n_targets);

  ktp_aux_t aux;
  aux.conf = conf;
  aux.data = data;
  aux.h = h;
  aux.n = n;
  aux.fn = fn;
  aux.bcf_hdr = bcf_hdr;
  aux.sm = sm;
  aux.rghash = rghash;
  aux.tid_idx = -1;
  aux.tid = -1;
  aux.beg = 0;
  aux.end = -1;
  aux.writer = &writer;
  // here will change the header
  init_aux_call(&aux);
  writer.WriteHeader();

  mem_process(&aux);

  bam_smpl_destroy(sm);
  bcf_call_del_rghash(rghash);
  mplp_ref_destory(&conf->ref, h->n_targets);
  bam_hdr_destroy(h);
  for (int i = 0; i < n; ++i) {
    sam_close(data[i]->fp);
    free(data[i]);
  }
  free(data);

  mcall_destroy(&aux.call);
  free(aux.call.ploidy);
  return 0;
}

static int is_url(const char *s)
{
  static const char uri_scheme_chars[] =
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+.-";
  return s[strspn(s, uri_scheme_chars)] == ':';
}

#define MAX_PATH_LEN 1024
int read_file_list(const char *file_list,int *n,char **argv[])
{
  char buf[MAX_PATH_LEN];
  int len, nfiles = 0;
  char **files = NULL;
  struct stat sb;

  *n = 0;
  *argv = NULL;

  FILE *fh = fopen(file_list,"r");
  if ( !fh )
  {
    fprintf(stderr,"%s: %s\n", file_list,strerror(errno));
    return 1;
  }

  files = (char**)calloc(nfiles,sizeof(char*));
  nfiles = 0;
  while ( fgets(buf,MAX_PATH_LEN,fh) )
  {
    // allow empty lines and trailing spaces
    len = strlen(buf);
    while ( len>0 && isspace(buf[len-1]) ) len--;
    if ( !len ) continue;

    // check sanity of the file list
    buf[len] = 0;
    if (! (is_url(buf) || stat(buf, &sb) == 0))
    {
      // no such file, check if it is safe to print its name
      int i, safe_to_print = 1;
      for (i=0; i<len; i++)
        if (!isprint(buf[i])) { safe_to_print = 0; break; }
      if ( safe_to_print )
        fprintf(stderr,"The file list \"%s\" appears broken, could not locate: %s\n", file_list,buf);
      else
        fprintf(stderr,"Does the file \"%s\" really contain a list of files and do all exist?\n", file_list);
      return 1;
    }

    nfiles++;
    files = (char**)realloc(files,nfiles*sizeof(char*));
    files[nfiles-1] = strdup(buf);
  }
  fclose(fh);
  if ( !nfiles )
  {
    fprintf(stderr,"No files read from %s\n", file_list);
    return 1;
  }
  *argv = files;
  *n    = nfiles;
  return 0;
}
#undef MAX_PATH_LEN

int parse_format_flag(const char *str)
{
  int i, flag = 0, n_tags;
  char **tags = hts_readlist(str, 0, &n_tags);
  for(i=0; i<n_tags; i++)
  {
    if ( !strcasecmp(tags[i],"DP") ) flag |= B2B_FMT_DP;
    else if ( !strcasecmp(tags[i],"DV") ) { flag |= B2B_FMT_DV; fprintf(stderr, "[warning] tag DV functional, but deprecated. Please switch to `AD` in future.\n"); }
    else if ( !strcasecmp(tags[i],"SP") ) flag |= B2B_FMT_SP;
    else if ( !strcasecmp(tags[i],"DP4") ) { flag |= B2B_FMT_DP4; fprintf(stderr, "[warning] tag DP4 functional, but deprecated. Please switch to `ADF` and `ADR` in future.\n"); }
    else if ( !strcasecmp(tags[i],"DPR") ) { flag |= B2B_FMT_DPR; fprintf(stderr, "[warning] tag DPR functional, but deprecated. Please switch to `AD` in future.\n"); }
    else if ( !strcasecmp(tags[i],"INFO/DPR") ) { flag |= B2B_INFO_DPR; fprintf(stderr, "[warning] tag INFO/DPR functional, but deprecated. Please switch to `INFO/AD` in future.\n"); }
    else if ( !strcasecmp(tags[i],"AD") ) flag |= B2B_FMT_AD;
    else if ( !strcasecmp(tags[i],"ADF") ) flag |= B2B_FMT_ADF;
    else if ( !strcasecmp(tags[i],"ADR") ) flag |= B2B_FMT_ADR;
    else if ( !strcasecmp(tags[i],"INFO/AD") ) flag |= B2B_INFO_AD;
    else if ( !strcasecmp(tags[i],"INFO/ADF") ) flag |= B2B_INFO_ADF;
    else if ( !strcasecmp(tags[i],"INFO/ADR") ) flag |= B2B_INFO_ADR;
    else
    {
      fprintf(stderr,"Could not parse tag \"%s\" in \"%s\"\n", tags[i], str);
      exit(EXIT_FAILURE);
    }
    free(tags[i]);
  }
  if (n_tags) free(tags);
  return flag;
}

static void print_usage(FILE *fp, const mplp_conf_t *mplp)
{
  char *tmp_require = bam_flag2str(mplp->rflag_require);
  char *tmp_filter  = bam_flag2str(mplp->rflag_filter);

  // Display usage information, formatted for the standard 80 columns.
  // (The unusual string formatting here aids the readability of this
  // source code in 80 columns, to the extent that's possible.)

  fprintf(fp,
          "\n"
          "Usage: samtools mpileup [options] in1.bam [in2.bam [...]]\n"
          "\n"
          "Input options:\n"
          "  -6, --illumina1.3+      quality is in the Illumina-1.3+ encoding\n"
          "  -A, --count-orphans     do not discard anomalous read pairs\n"
          "  -b, --bam-list FILE     list of input BAM filenames, one per line\n"
          "  -B, --no-BAQ            disable BAQ (per-Base Alignment Quality)\n"
          "  -C, --adjust-MQ INT     adjust mapping quality; recommended:50, disable:0 [0]\n"
          "  -d, --max-depth INT     max per-file depth; avoids excessive memory usage [%d]\n", mplp->max_depth);
  fprintf(fp,
          "  -E, --redo-BAQ          recalculate BAQ on the fly, ignore existing BQs\n"
          "  -f, --fasta-ref FILE    faidx indexed reference sequence file\n"
          "  -G, --exclude-RG FILE   exclude read groups listed in FILE\n"
          "  -l, --positions FILE    skip unlisted positions (chr pos) or regions (BED)\n"
          "  -q, --min-MQ INT        skip alignments with mapQ smaller than INT [%d]\n", mplp->min_mq);
  fprintf(fp,
          "  -Q, --min-BQ INT        skip bases with baseQ/BAQ smaller than INT [%d]\n", mplp->min_baseQ);
  fprintf(fp,
          "  -r, --region REG        region in which pileup is generated\n"
          "  -R, --ignore-RG         ignore RG tags (one BAM = one sample)\n"
          "  --rf, --incl-flags STR|INT  required flags: skip reads with mask bits unset [%s]\n", tmp_require);
  fprintf(fp,
          "  --ff, --excl-flags STR|INT  filter flags: skip reads with mask bits set\n"
          "                                            [%s]\n", tmp_filter);
  fprintf(fp,
          "  -x, --ignore-overlaps   disable read-pair overlap detection\n"
          "\n"
          "Output options:\n"
          "  -o, --output FILE       write output to FILE [standard output]\n"
          "  -g, --BCF               generate genotype likelihoods in BCF format\n"
          "  -v, --VCF               generate genotype likelihoods in VCF format\n"
          "\n"
          "Output options for mpileup format (without -g/-v):\n"
          "  -O, --output-BP         output base positions on reads\n"
          "  -s, --output-MQ         output mapping quality\n"
          "  -a                      output all positions (including zero depth)\n"
          "  -a -a (or -aa)          output absolutely all positions, including unused ref. sequences\n"
          "\n"
          "Output options for genotype likelihoods (when -g/-v is used):\n"
          "  -t, --output-tags LIST  optional tags to output:\n"
          "               DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR []\n"
          "  -u, --uncompressed      generate uncompressed VCF/BCF output\n"
          "\n"
          "SNP/INDEL genotype likelihoods options (effective with -g/-v):\n"
          "  -e, --ext-prob INT      Phred-scaled gap extension seq error probability [%d]\n", mplp->extQ);
  fprintf(fp,
          "  -F, --gap-frac FLOAT    minimum fraction of gapped reads [%g]\n", mplp->min_frac);
  fprintf(fp,
          "  -h, --tandem-qual INT   coefficient for homopolymer errors [%d]\n", mplp->tandemQ);
  fprintf(fp,
          "  -I, --skip-indels       do not perform indel calling\n"
          "  -L, --max-idepth INT    maximum per-file depth for INDEL calling [%d]\n", mplp->max_indel_depth);
  fprintf(fp,
          "  -m, --min-ireads INT    minimum number gapped reads for indel candidates [%d]\n", mplp->min_support);
  fprintf(fp,
          "  -o, --open-prob INT     Phred-scaled gap open seq error probability [%d]\n", mplp->openQ);
  fprintf(fp,
          "  -p, --per-sample-mF     apply -m and -F per-sample for increased sensitivity\n"
          "  -P, --platforms STR     comma separated list of platforms for indels [all]\n");
  fprintf(fp,
          "  -@, --threads INT       Number of BAM/CRAM compression threads [%d]\n"
          "  -j --max_lreads INT     maximum alignment length of reads, alignmentEnd - alignmentStart[%d]\n"
          "  -c, --chunk_size INT    Number of reads should be process in single worker[%d]\n"
          , mplp->nthreads, mplp->max_lreads, mplp->chunk_size);
  sam_global_opt_help(fp, "-.--.");
  fprintf(fp,
          "\n"
          "Notes: Assuming diploid individuals.\n");

  free(tmp_require);
  free(tmp_filter);
}

int main(int argc, char *argv[])
{
  int c;
  const char *file_list = NULL;
  char **fn = NULL;
  int nfiles = 0, use_orphan = 0;
  mplp_conf_t mplp;
  memset(&mplp, 0, sizeof(mplp_conf_t));
  mplp.min_baseQ = 13;
  mplp.capQ_thres = 0;
  mplp.max_depth = 250; mplp.max_indel_depth = 250;
  mplp.openQ = 40; mplp.extQ = 20; mplp.tandemQ = 100;
  mplp.min_frac = 0.002; mplp.min_support = 1;
  mplp.flag = MPLP_NO_ORPHAN | MPLP_REALN | MPLP_SMART_OVERLAPS;
  mplp.argc = argc; mplp.argv = argv;
  mplp.rflag_filter = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP;
  mplp.output_fname = NULL;
  mplp.all = 0;
  mplp.nthreads = 1;
  // The default max read length is 300
  mplp.max_lreads = 300;
  mplp.chunk_size = 100000;
  sam_global_args_init(&mplp.ga);

  static const struct option lopts[] =
  {
    SAM_OPT_GLOBAL_OPTIONS('-', 0, '-', '-', 0),
    {"rf", required_argument, NULL, 1},   // require flag
    {"ff", required_argument, NULL, 2},   // filter flag
    {"incl-flags", required_argument, NULL, 1},
    {"excl-flags", required_argument, NULL, 2},
    {"output", required_argument, NULL, 3},
    {"open-prob", required_argument, NULL, 4},
    {"illumina1.3+", no_argument, NULL, '6'},
    {"count-orphans", no_argument, NULL, 'A'},
    {"bam-list", required_argument, NULL, 'b'},
    {"no-BAQ", no_argument, NULL, 'B'},
    {"no-baq", no_argument, NULL, 'B'},
    {"adjust-MQ", required_argument, NULL, 'C'},
    {"adjust-mq", required_argument, NULL, 'C'},
    {"max-depth", required_argument, NULL, 'd'},
    {"redo-BAQ", no_argument, NULL, 'E'},
    {"redo-baq", no_argument, NULL, 'E'},
    {"fasta-ref", required_argument, NULL, 'f'},
    {"exclude-RG", required_argument, NULL, 'G'},
    {"exclude-rg", required_argument, NULL, 'G'},
    {"positions", required_argument, NULL, 'l'},
    {"region", required_argument, NULL, 'r'},
    {"ignore-RG", no_argument, NULL, 'R'},
    {"ignore-rg", no_argument, NULL, 'R'},
    {"min-MQ", required_argument, NULL, 'q'},
    {"min-mq", required_argument, NULL, 'q'},
    {"min-BQ", required_argument, NULL, 'Q'},
    {"min-bq", required_argument, NULL, 'Q'},
    {"ignore-overlaps", no_argument, NULL, 'x'},
    {"BCF", no_argument, NULL, 'g'},
    {"bcf", no_argument, NULL, 'g'},
    {"VCF", no_argument, NULL, 'v'},
    {"vcf", no_argument, NULL, 'v'},
    {"output-BP", no_argument, NULL, 'O'},
    {"output-bp", no_argument, NULL, 'O'},
    {"output-MQ", no_argument, NULL, 's'},
    {"output-mq", no_argument, NULL, 's'},
    {"output-tags", required_argument, NULL, 't'},
    {"uncompressed", no_argument, NULL, 'u'},
    {"ext-prob", required_argument, NULL, 'e'},
    {"gap-frac", required_argument, NULL, 'F'},
    {"tandem-qual", required_argument, NULL, 'h'},
    {"skip-indels", no_argument, NULL, 'I'},
    {"max-idepth", required_argument, NULL, 'L'},
    {"min-ireads ", required_argument, NULL, 'm'},
    {"per-sample-mF", no_argument, NULL, 'p'},
    {"per-sample-mf", no_argument, NULL, 'p'},
    {"platforms", required_argument, NULL, 'P'},
    {"threads", required_argument, NULL, '@'},
    {"max_lreads", required_argument, NULL, 'j'},
    {"chunk_size", required_argument, NULL, 'c'},
    {NULL, 0, NULL, 0}
  };
  while ((c = getopt_long(argc, argv, "Agf:r:l:q:Q:uRC:BDSd:L:b:P:po:e:h:Im:F:EG:6OsVvxt:a@:j:c:",lopts,NULL)) >= 0) {
    switch (c) {
      case 'x': mplp.flag &= ~MPLP_SMART_OVERLAPS; break;
      case  1 :
                mplp.rflag_require = bam_str2flag(optarg);
                if ( mplp.rflag_require<0 ) { fprintf(stderr,"Could not parse --rf %s\n", optarg); return 1; }
                break;
      case  2 :
                mplp.rflag_filter = bam_str2flag(optarg);
                if ( mplp.rflag_filter<0 ) { fprintf(stderr,"Could not parse --ff %s\n", optarg); return 1; }
                break;
      case  3 : mplp.output_fname = optarg; break;
      case  4 : mplp.openQ = atoi(optarg); break;
      case 'f':
                mplp.fai = fai_load(optarg);
                if (mplp.fai == NULL) return 1;
                mplp.fai_fname = optarg;
                break;
      case 'd': mplp.max_depth = atoi(optarg); break;
      case 'r': mplp.reg = strdup(optarg); break;
      case 'l':
                // In the original version the whole BAM was streamed which is inefficient
                //  with few BED intervals and big BAMs. Todo: devise a heuristic to determine
                //  best strategy, that is streaming or jumping.
                mplp.bed = bed_read(optarg);
                if (!mplp.bed) { print_error_errno("mpileup", "Could not read file \"%s\"", optarg); return 1; }
                break;
      case 'P': mplp.pl_list = strdup(optarg); break;
      case 'p': mplp.flag |= MPLP_PER_SAMPLE; break;
      case 'g': mplp.flag |= MPLP_BCF; break;
      case 'v': mplp.flag |= MPLP_BCF | MPLP_VCF; break;
      case 'u': mplp.flag |= MPLP_NO_COMP | MPLP_BCF; break;
      case 'B': mplp.flag &= ~MPLP_REALN; break;
      case 'D': mplp.fmt_flag |= B2B_FMT_DP; fprintf(stderr, "[warning] samtools mpileup option `-D` is functional, but deprecated. Please switch to `-t DP` in future.\n"); break;
      case 'S': mplp.fmt_flag |= B2B_FMT_SP; fprintf(stderr, "[warning] samtools mpileup option `-S` is functional, but deprecated. Please switch to `-t SP` in future.\n"); break;
      case 'V': mplp.fmt_flag |= B2B_FMT_DV; fprintf(stderr, "[warning] samtools mpileup option `-V` is functional, but deprecated. Please switch to `-t DV` in future.\n"); break;
      case 'I': mplp.flag |= MPLP_NO_INDEL; break;
      case 'E': mplp.flag |= MPLP_REDO_BAQ; break;
      case '6': mplp.flag |= MPLP_ILLUMINA13; break;
      case 'R': mplp.flag |= MPLP_IGNORE_RG; break;
      case 's': mplp.flag |= MPLP_PRINT_MAPQ; break;
      case 'O': mplp.flag |= MPLP_PRINT_POS; break;
      case 'C': mplp.capQ_thres = atoi(optarg); break;
      case 'q': mplp.min_mq = atoi(optarg); break;
      case 'Q': mplp.min_baseQ = atoi(optarg); break;
      case 'b': file_list = optarg; break;
      case 'o': {
        char *end;
        long value = strtol(optarg, &end, 10);
        // Distinguish between -o INT and -o FILE (a bit of a hack!)
        if (*end == '\0') mplp.openQ = value;
        else mplp.output_fname = optarg;
      }
                break;
      case 'e': mplp.extQ = atoi(optarg); break;
      case 'h': mplp.tandemQ = atoi(optarg); break;
      case 'A': use_orphan = 1; break;
      case 'F': mplp.min_frac = atof(optarg); break;
      case 'm': mplp.min_support = atoi(optarg); break;
      case 'L': mplp.max_indel_depth = atoi(optarg); break;
      case 'G': {
        FILE *fp_rg;
        char buf[1024];
        mplp.rghash = khash_str2int_init();
        if ((fp_rg = fopen(optarg, "r")) == NULL)
          fprintf(stderr, "(%s) Fail to open file %s. Continue anyway.\n", __func__, optarg);
        while (!feof(fp_rg) && fscanf(fp_rg, "%s", buf) > 0) // this is not a good style, but forgive me...
          khash_str2int_inc(mplp.rghash, strdup(buf));
        fclose(fp_rg);
      }
                break;
      case 't': mplp.fmt_flag |= parse_format_flag(optarg); break;
      case 'a': mplp.all++; break;
      case '@': mplp.nthreads = atoi(optarg); break;
      case 'j': mplp.max_lreads = atoi(optarg); break;
      case 'c': mplp.chunk_size = atoi(optarg); break;
      default:
                if (parse_sam_global_opt(c, optarg, lopts, &mplp.ga) == 0) break;
                /* else fall-through */
      case '?':
                print_usage(stderr, &mplp);
                return 1;
    }
  }
  if (!mplp.fai && mplp.ga.reference) {
    mplp.fai_fname = mplp.ga.reference;
    mplp.fai = fai_load(mplp.fai_fname);
    if (mplp.fai == NULL) return 1;
  }
  pthread_mutex_init(&mplp.mutex, 0);

  if ( !(mplp.flag&MPLP_REALN) && mplp.flag&MPLP_REDO_BAQ )
  {
    fprintf(stderr,"Error: The -B option cannot be combined with -E\n");
    return 1;
  }
  if (use_orphan) mplp.flag &= ~MPLP_NO_ORPHAN;
  if (argc == 1)
  {
    print_usage(stderr, &mplp);
    return 1;
  }
  if (mplp.nthreads < 1) mplp.nthreads = 1;
  int ret;
  if (file_list) {
    if ( read_file_list(file_list,&nfiles,&fn) ) return 1;
    ret = mpileup(&mplp,nfiles,fn);
    for (c=0; c<nfiles; c++) free(fn[c]);
    free(fn);
  }
  else
    ret = mpileup(&mplp, argc - optind, argv + optind);
  if (mplp.rghash) khash_str2int_destroy_free(mplp.rghash);
  free(mplp.reg); free(mplp.pl_list);
  if (mplp.fai) fai_destroy(mplp.fai);
  pthread_mutex_destroy(&mplp.mutex);

  if (mplp.bed) bed_destroy(mplp.bed);
  return ret;
}

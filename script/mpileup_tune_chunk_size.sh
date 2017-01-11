#!/bin/bash
#
# Tune the best chunk_size of mpileup
# Author: PingZeng(zengping@ncic.ac.cn)
#

# XXX the dir of the current script
SCRIPT_DIR=$(cd $(dirname $(readlink -f $0)); pwd )


# Set the following path
BASE_DIR=/home/zp/work/rnaseq/cpp-htslib/test/
DATA_DIR=${BASE_DIR}/localtestdata/mpileup/
FASTA=${BASE_DIR}/localtestdata/references/references/ucsc.hg19.fa
BAMFILES="${DATA_DIR}/accepted_hits_sorted_chr1.bam"

BIN=${SCRIPT_DIR}/../build/debug/bin/mpileup_bcftools_region

chunk_sizes=(10000 25000 50000 100000 200000 400000 800000)
for chunk_size in ${chunk_sizes[@]}
do
  echo "==========chunk_size: $chunk_size start============"
  date +'%Y-%m-%d %H:%M:%S'
  ${BIN} -uvf $FASTA $BAMFILES \
  -@ 24 --chunk_size $chunk_size \
  -o ${DATA_DIR}/tune_chunk_size.vcf &> /dev/null

  date +'%Y-%m-%d %H:%M:%S'
  rm ${DATA_DIR}/tune_chunk_size.vcf
  echo "==========chunk_size: $chunk_size end============"
done


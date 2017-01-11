#!/bin/bash
#
# Check the correstness of mutect
# Author: PingZeng(zengping@ncic.ac.cn)
#

# XXX the dir of the current script
SCRIPT_DIR=$(cd $(dirname $(readlink -f $0)); pwd )

# Set the multithread you want to check correstness
nthreads=(3)

#########################################################################

BIN=${SCRIPT_DIR}/../build/debug/bin/mpileup_bcftools_region

DATA_DIR=${SCRIPT_DIR}/../test/pileup/
FASTA=${DATA_DIR}/mpileup.ref.fa
BAMFILES="${DATA_DIR}/mpileup.1.bam ${DATA_DIR}/mpileup.2.bam ${DATA_DIR}/mpileup.3.bam"

for thread in ${nthreads[@]}
do
  echo "=================Check mpileup correstness on thread:$thread==============="
  ${BIN} -uvf $FASTA $BAMFILES \
  -@ 24 --chunk_size 1000 \
  -o ${DATA_DIR}/out.vcf &> /dev/null

  sort -k2n ${DATA_DIR}/out.vcf > ${DATA_DIR}/out-sorted.vcf

  diff ${DATA_DIR}/out-sorted.vcf ${DATA_DIR}/result.vcf

  if [ $? -eq 0 ]; then
    echo "=================mpileup correst on thread:$thread============="
    rm ${DATA_DIR}/out*.vcf
  else
    echo "=================mpileup not correst on thread:$thread============="
    exit
  fi
done


#!/bin/bash
#
# Check the correstness of mutect
# Author: PingZeng(zengping@ncic.ac.cn)
#

# XXX the dir of the current script
SCRIPT_DIR=$(cd $(dirname $(readlink -f $0)); pwd )

# Set the following path
REFERENCE=/home/zp/work/rnaseq/cpp-htslib/test/localtestdata/references/references/human_g1k_v37.fasta
DBSNP=/home/zp/work/rnaseq/cpp-htslib/test/localtestdata/references/vcf/dbsnp_138.b37.compressed.vcf
COSMIC=/home/zp/work/rnaseq/cpp-htslib/test/localtestdata/references/vcf/b37_cosmic_v54_120711.compressed.vcf

# Set the multithread you want to check correstness
nthreads=(24 12)

#########################################################################

BIN=${SCRIPT_DIR}/../build/debug/bin/mutect
DATA_DIR=${SCRIPT_DIR}/../test/mutect/

for thread in ${nthreads[@]}
do
  echo "=================Check mutect correstness on thread:$thread==============="
  ${BIN} \
    --intervals ${DATA_DIR}/mutect_panle_372.interval_list \
    --reference_sequence ${REFERENCE} \
    --nthreads $thread \
    --fraction_contamination 0.00 \
    --dbsnp ${DBSNP} \
    --cosmic ${COSMIC} \
    --normal_file ${DATA_DIR}/MG225_normal_sorted_X.bam \
    --tumor_file ${DATA_DIR}/MG225_tumor_sorted_X.bam \
    --out ${DATA_DIR}/out-cpp.txt \
    --vcf ${DATA_DIR}/out-cpp.vcf &> /dev/null

  sort -k2n ${DATA_DIR}/out-cpp.txt > ${DATA_DIR}/out-sorted.txt
  sort -k2n ${DATA_DIR}/out-cpp.vcf > ${DATA_DIR}/out-sorted.vcf

  diff ${DATA_DIR}/out-sorted.txt ${DATA_DIR}/result.txt &&
    diff ${DATA_DIR}/out-sorted.vcf ${DATA_DIR}/result.vcf

  if [ $? -eq 0 ]; then
    echo "=================mutect correst on thread:$thread============="
    rm ${DATA_DIR}/out-*.txt ${DATA_DIR}/out-*.vcf
  else
    echo "=================mutect not correst on thread:$thread============="
    exit
  fi
done


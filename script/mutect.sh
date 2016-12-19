USE_LARGE=1

BASE_DIR=/home/zp/work/rnaseq/cpp-htslib/test/localtestdata/
if [ $USE_LARGE -eq 1 ]; then
  DATA_DIR=${BASE_DIR}/large/
  normal=MG225_normal_sorted.bam
  tumor=MG225_tumor_sorted.bam
else
  DATA_DIR=${BASE_DIR}/middle/
  normal=MG225_normal_sorted_X.bam
  tumor=MG225_tumor_sorted_X.bam
fi

date +'%Y-%m-%d %H:%M:%S'
./bin/mutect \
 --intervals ${DATA_DIR}/mutect_panle_372.interval_list \
 --reference_sequence /mnt/sdb_1T/GATK_Reference/reference/human_g1k_v37.fasta \
 --nthreads 48 \
 --fraction_contamination 0.00 \
 --dbsnp ${BASE_DIR}/references/vcf/dbsnp_138.b37.compressed.vcf \
 --cosmic ${BASE_DIR}/references/vcf/b37_cosmic_v54_120711.compressed.vcf \
 --normal_file ${DATA_DIR}/${normal} \
 --tumor_file ${DATA_DIR}/${tumor} \
 --out ${DATA_DIR}/out.txt \
 --vcf ${DATA_DIR}/out.vcf
date +'%Y-%m-%d %H:%M:%S'

USE_LARGE=1

if [ $USE_LARGE -eq 1 ]; then
  BASE_DIR=/home/zp/work/rnaseq/cpp-htslib/test/large/
  normal=MG225_normal_sorted.bam
  tumor=MG225_tumor_sorted.bam
else
  BASE_DIR=/home/zp/work/rnaseq/cpp-htslib/test/
  normal=MG225_normal_sorted_X.bam
  tumor=MG225_tumor_sorted_X.bam
fi

date +'%Y-%m-%d %H:%M:%S'
./bin/mutect \
 --intervals ${BASE_DIR}/mutect_panle_372.interval_list \
 --reference_sequence /mnt/sdb_1T/GATK_Reference/reference/human_g1k_v37.fasta \
 --nthreads 48 \
 --downsampling \
 --normal_file ${BASE_DIR}/${normal} \
 --tumor_file ${BASE_DIR}/${tumor} \
 --out ${BASE_DIR}/out.txt \
 --vcf ${BASE_DIR}/out.vcf
date +'%Y-%m-%d %H:%M:%S'

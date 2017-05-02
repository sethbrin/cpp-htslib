USE_LARGE=0

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

echo "==================mutect origin start====="
date +'%Y-%m-%d %H:%M:%S'

#java -Xmx2g -jar ${DATA_DIR}/mutect-1.1.7.jar \
#  --analysis_type MuTect \
#  --downsampling_type NONE \
#  --intervals ${DATA_DIR}/mutect_panle_372.interval_list \
#  --reference_sequence /home/zp/work/rnaseq/cpp-htslib/test/localtestdata/references/references/human_g1k_v37.fasta \
#  --fraction_contamination 0.00 \
#  --dbsnp ${BASE_DIR}/references/vcf/dbsnp_138.b37.vcf \
#  --cosmic ${BASE_DIR}/references/vcf/b37_cosmic_v54_120711.vcf \
#  --input_file:normal ${DATA_DIR}/${normal} \
#  --input_file:tumor ${DATA_DIR}/${tumor} \
#  --out ${DATA_DIR}/out-java.txt \
#  --vcf ${DATA_DIR}/out-java.vcf

echo "==================mutect origin end====="
date +'%Y-%m-%d %H:%M:%S'

echo "==================mutect start====="
date +'%Y-%m-%d %H:%M:%S'

./bin/mutect \
  --intervals ${DATA_DIR}/mutect_panle_372.interval_list \
  --reference_sequence /home/zp/work/rnaseq/cpp-htslib/test/localtestdata/references/references/human_g1k_v37.fasta \
  --nthreads 12 \
  --fraction_contamination 0.00 \
  --dbsnp ${BASE_DIR}/references/vcf/dbsnp_138.b37.compressed.vcf \
  --cosmic ${BASE_DIR}/references/vcf/b37_cosmic_v54_120711.compressed.vcf \
  --normal_file ${DATA_DIR}/${normal} \
  --tumor_file ${DATA_DIR}/${tumor} \
  --out ${DATA_DIR}/out-cpp.txt \
  --vcf ${DATA_DIR}/out-cpp.vcf

date +'%Y-%m-%d %H:%M:%S'
echo "==================mutect end====="

USE_LARGE=1

BASE_DIR=/home/zp/work/rnaseq/cpp-htslib/test/

if [ $USE_LARGE -eq 1 ]; then
  DATA_DIR=${BASE_DIR}/localtestdata/mpileup/
  FASTA=/home/zp/work/rnaseq/cpp-htslib/test/localtestdata/references/references/ucsc.hg19.fa
  BAMFILES="${DATA_DIR}/accepted_hits_sorted_chr1.small.bam"
else
  DATA_DIR=${BASE_DIR}/pileup/
  FASTA=${DATA_DIR}/mpileup.ref.fa
  BAMFILES="${DATA_DIR}/mpileup.1.bam ${DATA_DIR}/mpileup.2.bam ${DATA_DIR}/mpileup.3.bam"
fi

date +'%Y-%m-%d %H:%M:%S'
./bin/mpileup -uvf $FASTA $BAMFILES -@ 2 -o ${DATA_DIR}/out.vcf
#${DATA_DIR}/samtools mpileup -uvf $FASTA $BAMFILES -o ${DATA_DIR}/mpileup.vcf

#${BASE_DIR}/samtools mpileup -ugf $FASTA $BAMFILES \
#  | ${BASE_DIR}/bcftools call -O v -o ${DATA_DIR}/out.vcf -vm
date +'%Y-%m-%d %H:%M:%S'

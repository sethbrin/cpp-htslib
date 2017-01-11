USE_LARGE=1

BASE_DIR=/home/zp/work/rnaseq/cpp-htslib/test/

if [ $USE_LARGE -eq 1 ]; then
  DATA_DIR=${BASE_DIR}/localtestdata/mpileup/
  FASTA=/home/zp/work/rnaseq/cpp-htslib/test/localtestdata/references/references/ucsc.hg19.fa
  BAMFILES="${DATA_DIR}/accepted_hits_sorted_chr1.bam"
  #BAMFILES="${DATA_DIR}/accepted_hits_sorted.bam"
  #BAMFILES="${DATA_DIR}/accepted_hits_sorted_chr1.small.bam"
else
  DATA_DIR=${BASE_DIR}/pileup/
  FASTA=${DATA_DIR}/mpileup.ref.fa
  BAMFILES="${DATA_DIR}/mpileup.1.bam ${DATA_DIR}/mpileup.2.bam ${DATA_DIR}/mpileup.3.bam"
fi

#chunk_sizes=(10000 25000 50000 100000 200000 400000 800000)
chunk_sizes=(400000 500000)
for chunk_size in ${chunk_sizes[@]}
do
  echo "==========chunk_size: $chunk_size start============"
  date +'%Y-%m-%d %H:%M:%S'
  ./bin/mpileup_bcftools_region -uvf $FASTA $BAMFILES -@ 24 --chunk_size $chunk_size -o ${DATA_DIR}/mpileup_bcftools_region-${chunk_size}.vcf
  date +'%Y-%m-%d %H:%M:%S'
  echo "==========chunk_size: $chunk_size end============"
done
exit
#${DATA_DIR}/samtools mpileup -ugf $FASTA $BAMFILES \
#  | ${DATA_DIR}/bcftools call -O v -o ${DATA_DIR}/mpileup.vcf -vm
#${DATA_DIR}/samtools mpileup -uvf $FASTA $BAMFILES -o ${DATA_DIR}/mpileup-all.vcf
./bin/mpileup_bcftools_region -uvf $FASTA $BAMFILES -@ 24 --chunk_size 50000 -o ${DATA_DIR}/mpileup_bcftools_region.vcf
#cat ${DATA_DIR}/mpileup2.vcf \
#  | ../../../bcftools/bcftools call -O v -o ${DATA_DIR}/mpileup.vcf -vm
date +'%Y-%m-%d %H:%M:%S'
exit

echo "mpileup2 thread1 end"
date +'%Y-%m-%d %H:%M:%S'
./bin/mpileup2 -uvf $FASTA $BAMFILES -@ 1 --chunk_size 500000 \
  | ${DATA_DIR}/bcftools call -O v -o ${DATA_DIR}/out.vcf -vm
date +'%Y-%m-%d %H:%M:%S'
echo "mpileup2 thread end"


echo "samtools mpileup thread4 start"
date +'%Y-%m-%d %H:%M:%S'
${DATA_DIR}/samtools mpileup -ugf $FASTA $BAMFILES \
  | ${DATA_DIR}/bcftools call -O v -o ${DATA_DIR}/mpileup.vcf -vm
date +'%Y-%m-%d %H:%M:%S'
echo "samtools mpileup thread4 end"

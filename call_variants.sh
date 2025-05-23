REF=$1
recal_dir=$2
vcf_dir=$3
mkdir -p $vcf_dir
for bam in $recal_dir/*.recal.bam
do
  sample=$(basename $bam .recal.bam)
  gatk HaplotypeCaller -R $REF -I $bam -O $vcf_dir/${sample}.vcf.gz
  tabix -p vcf $vcf_dir/${sample}.vcf.gz
done

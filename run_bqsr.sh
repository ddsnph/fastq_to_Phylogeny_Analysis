REF=$1
bam_dir=$2
recal_dir=$3
sites_dir=$4
mkdir -p $recal_dir
for bam in $bam_dir/*.dedup.bam
do
  sample=$(basename $bam .dedup.bam)
  gatk BaseRecalibrator -R $REF -I $bam \
    --known-sites $sites_dir/hapmap_3.3.hg38.vcf.gz \
    --known-sites $sites_dir/omni2.5.hg38.vcf.gz \
    --known-sites $sites_dir/phase1.snps.high_confidence.hg38.vcf.gz \
    --known-sites $sites_dir/Mills_and_1000G_indels.hg38.vcf.gz \
    -O $recal_dir/${sample}.recal.table
  gatk ApplyBQSR -R $REF -I $bam --bqsr-recal-file $recal_dir/${sample}.recal.table \
    -O $recal_dir/${sample}.recal.bam
done

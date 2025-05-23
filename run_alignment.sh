REF=$1
fastq_dir=$2
bam_dir=$3
mkdir -p $bam_dir
for S in $(ls $fastq_dir/*_1.fastq.gz | sed 's/_1.fastq.gz//')
do
  sample=$(basename $S)
  bwa mem $REF ${S}_1.fastq.gz ${S}_2.fastq.gz \
    | samtools view -bS - \
    | samtools sort -o $bam_dir/${sample}.sorted.bam
  samtools index $bam_dir/${sample}.sorted.bam
  picard MarkDuplicates I=$bam_dir/${sample}.sorted.bam \
    O=$bam_dir/${sample}.dedup.bam M=$bam_dir/${sample}.metrics.txt REMOVE_DUPLICATES=true
  samtools index $bam_dir/${sample}.dedup.bam
done

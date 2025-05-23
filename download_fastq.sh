cancer_list=$1
normal_list=$2
outdir=$3
mkdir -p $outdir
for S in $(cat $cancer_list $normal_list)
do
  fasterq-dump $S --split-files -O $outdir
  gzip $outdir/${S}_1.fastq
  gzip $outdir/${S}_2.fastq
done

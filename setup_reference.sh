REF=$1
mkdir -p resources
cp $REF resources/hs38.fa
bwa index resources/hs38.fa
samtools faidx resources/hs38.fa

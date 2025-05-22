# fastq_to_Phylogeny_Analysis
This is a repository showing the process I took to take fastq files to vcf files to perform phylogeny analysis with them.

Steps:
I used these steps from this website to guide my steps: https://www.hadriengourle.com/wrangling-genomics/aio/
Then I cloned this repository: https://github.com/MU-Data-Science/GAF/tree/main/post-vcf-pipeline

From here: I created a conda environemnt with bwa, samtools, gatk4, picard, sra-tools, python version 3.8, bcftools, tabix, iqtree2, and vcf-kit 
I actitvated and worked in this environment

I grabbed cancer data fastq from here, I just needed the SRA number: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA820526
I grabbed normal data fastq from here, once again I just needed the SRA number: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SAMN00000538&o=acc_s%3Aa

Within the data directory I created a list 

**Major Issues I ran Into**

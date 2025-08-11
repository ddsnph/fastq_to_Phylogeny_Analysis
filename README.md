# fastq_to_Phylogeny_Analysis

This repo documents how I went from FASTQ reads to per-sample VCFs, and then attempted a phylogeny from those VCFs. I ran everything on CloudLab and used conda for tools.

## References I followed

* Wrangling Genomics notes I used as a guide  
  https://www.hadriengourle.com/wrangling-genomics/aio/
* Post-VCF pipeline ideas and structure  
  https://github.com/MU-Data-Science/GAF/tree/main/post-vcf-pipeline

## Data sources

* Cancer FASTQ accessions: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA820526  
* Normal FASTQ accessions: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SAMN00000538&o=acc_s%3Aa

I saved IDs in two files under `data/`:

* `data/cancer_acc.txt`  
* `data/normal_acc.txt`  

One accession per line.

## Environment

I created a conda env with: BWA, Samtools, GATK4, Picard, SRA-tools, Python 3.8, bcftools, tabix, IQ-TREE 2, vcf-kit.

Compute ran on a CloudLab node.

## Folder layout I used

```text
data/          # FASTQs and accession lists
resources/     # hs38 reference and known-sites VCFs
results/
  ├── bam/     # sorted, deduped, recalibrated BAMs
  └── vcf/     # per-sample VCFs and .tbi
scripts/       # shell scripts I ran (listed below)
```
---

## Scripts I Ran and What They Do
All scripts live in `scripts/`.

### `download_fastq.sh`
- Reads IDs from `data/cancer_acc.txt` and `data/normal_acc.txt`
- Downloads with SRA Toolkit (`fasterq-dump --split-files`)
- Gzips the FASTQs

### `setup_reference.sh`
- Prepares the hs38 reference in `resources/hs38.fa`
- Builds indexes:
  - `bwa index`
  - `samtools faidx`
  - Picard `CreateSequenceDictionary`
- Known sites in `resources/` (with `.tbi` index files):
  - `hapmap_3.3.hg38.vcf.gz`
  - `omni2.5.hg38.vcf.gz`
  - `phase1.snps.high_confidence.hg38.vcf.gz`
  - `Mills_and_1000G_indels.hg38.vcf.gz`

### `run_alignment.sh`
- Aligns each sample with BWA-MEM to hs38
- Converts to BAM, sorts, and indexes with Samtools
- Marks duplicates with Picard
- Outputs to `results/bam/<sample>/`

### `run_bqsr.sh`
- Runs GATK BaseRecalibrator with the known sites
- Applies BQSR to create `<sample>.recal.bam`

### `call_variants.sh`
- Runs GATK HaplotypeCaller per sample
- Produces compressed and indexed VCFs in `results/vcf/`

---

## What Worked
- Read downloads completed
- hs38 reference and indexes created
- Alignments finished, duplicates marked, BQSR applied
- Per-sample VCFs written to `results/vcf/` with `.tbi` indexes

---

## What I Attempted Next for Phylogeny
**Planned steps:**
1. Merge per-sample VCFs with `bcftools merge` to get a multi-sample VCF
2. Convert the merged VCF to an alignment of variant sites
3. Use vcf-kit to export a sites-only alignment, or
4. Build per-sample consensus FASTA from the merged VCF and align, then run IQ-TREE
5. Build a tree with IQ-TREE 2 using model GTR+G and ultrafast bootstrap

---

## What Blocked Completion
- Project storage ran out during downstream steps, causing:
  - Truncated or incomplete files during `bcftools merge`
  - Tabix index failures
  - Segmentation faults in `bcftools` due to truncated inputs
- Whole-genome consensus FASTAs were very large, which increased space pressure
- IQ-TREE failed when fed complete genome consensus files produced this way
- Since a clean merged VCF and alignment were not created, the tree step could not run

---

## Notes for Reproduction
- Watch free space closely; remove `fasterq.tmp.*` and intermediate BAMs once VCFs exist
- Make sure `hs38.fa` has `.fai` and `.dict`, and every known-sites VCF has a `.tbi`
- If `bcftools merge` warns about a missing BGZF EOF marker or segfaults:
  - One or more inputs are truncated
  - Recreate that VCF after freeing space
  - Re-index with `tabix` before merging
- For phylogeny, prefer a sites-only alignment derived from the merged VCF rather than whole-genome consensus FASTAs

---

## Acknowledgements
- CloudLab for compute resources  
- Wrangling Genomics notes: [https://www.hadriengourle.com/wrangling-genomics/aio/](https://www.hadriengourle.com/wrangling-genomics/aio/)  
- MU Data Science post-VCF pipeline: [https://github.com/MU-Data-Science/GAF/tree/main/post-vcf-pipeline](https://github.com/MU-Data-Science/GAF/tree/main/post-vcf-pipeline)

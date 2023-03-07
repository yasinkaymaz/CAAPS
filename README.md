# Cancer Associated Alternative Polyadenylations (CAAPS) Project

The main goal of this project is to determine alternative polyadenylation sites of transcripts specific to certain cell types in the tumor micro-environment using single-cell RNA sequencing data.

The project utilizes the current analysis pipeline written as snakemake rule to perform several scripts designed to detect cell type specific alternative polyadenylation switches between coditions of samples.

Current version can:

- Download from SRA database of NCBI.
- Align reads to reference transcriptome using `STAR solo` (needs >30GB of RAM).
- Run standard `Seurat` pipeline to normalize data and generate de novo cell clusters.
- (Optional) Identify cell types with a reference dataset using `HieRFIT`.
- Detect and eliminate internal priming events (IPEs) due to arbitrary oligo-dT hybridizations (utilizes `cleanUpdTseq`).
- Calls PAS peaks after Poisson denoising.
- Provides test statistics for switching PAS using NxN Chi-square.
- Determines genes with 3'-UTR shorthening and lenghtening. 


Needs conda and snakemake to run:

<code> conda install -c bioconda snakemake </code>

To run the analysis:

<code>snakemake --cores 2 -s SnakeFile.smk --use-conda 
snakemake --cores 2 -s SplitBams.smk --use-conda 
snakemake --cores 2 -s ClustersPAS.smk --use-conda 
bash scripts/make_PAS_table.v3.sh
Rscript scripts/PASseq_chisquire_beta_v4.R <code>

You can adjust the config file to change the accession entries

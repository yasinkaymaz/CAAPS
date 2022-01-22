
#Download example toy dataset:

wget http://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/2.1.0/pbmc8k/pbmc8k_possorted_genome_bam.bam


cellranger mkfastq --id=tiny-bcl \
                     --run=~/Documents/CAAPS/exampleData/cellranger-tiny-bcl-1.2.0/ \
                     --csv=cellranger-tiny-bcl-simple-1.2.0.csv

#https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/bam

#https://davetang.org/muse/2018/06/06/10x-single-cell-bam-files/

conda activate CAAPS
fasterq-dump SRR8325947 --split-files \
--threads 6 \
--mem 2048MB \
--progress \
--details

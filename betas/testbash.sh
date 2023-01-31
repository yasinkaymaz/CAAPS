
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


STAR --genomeDir starsolo \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist 10x_V3_whitelist.txt \
    --soloUMIlen 12 \
    --readFilesIn ${wd}/2270183_P7_2_S2_L001_R2_001.fastq.gz,${wd}/2270183_P7_2_S2_L002_R2_001.fastq.gz ${wd}/2270183_P7_2_S2_L001_R1_001.fastq.gz,${wd}/2270183_P7_2_S2_L002_R1_001.fastq.gz \
    --runThreadN 20 \
    --outFileNamePrefix s1 \
    --outSAMtype BAM SortedByCoordinate \
    --outReadsUnmapped elp1_s2_Unmapped \
    --twopassMode Basic \
    --chimSegmentMin 20 \
    --readFilesCommand zcat \
    --clipAdapterType CellRanger4 \
    --outFilterScoreMin 30 \
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
    --soloUMIfiltering MultiGeneUMI_CR \
    --soloUMIdedup 1MM_CR



python /Users/yasinkaymaz/Dropbox/codes/endSeqTools/endSeq_tools.py makeSamplesTable \
-sf /Users/yasinkaymaz/Documents/data/PASseqData/file_contains_sample_fastqs_ex.txt \
--gtf /Users/yasinkaymaz/Documents/data/gencode.v42.annotation.gtf

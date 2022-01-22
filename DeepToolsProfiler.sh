
conda activate DAB



bamCoverage -b output.sorted.bam \
-o output.sorted.bw \
--binSize 20 \
--normalizeUsing CPM \
--smoothLength 60 \
--extendReads 150 \
--centerReads \
-p 6 2> bamCoverage.log



computeMatrix scale-regions \
-R ~/Documents/data/GRCh38_Genes_GENCODEv38.bed.gz \
-S ~/Documents/CAAPS/output.sorted.bw \
--skipZeros \
-p 6 \
--regionBodyLength 2000 \
-a 500 -b 500 \
-o ~/Documents/CAAPS/output_scRNAseq_sites.gz

plotProfile -m ~/Documents/CAAPS/output_scRNAseq_sites.gz \
-out scRNAseq_profile.png \
--perGroup  --plotTitle "" \
--samplesLabel "scRNAseq" \
-T "scRNAseq read alignment sites"  -z "" \
--startLabel "" \
--endLabel "" \
--colors red

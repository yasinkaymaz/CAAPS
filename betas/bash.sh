python ~/Dropbox/codes/endSeqTools/endSeq_tools.py trimPolyAstrecth -f TDP120h1.fastq


python ~/Dropbox/codes/endSeqTools/endSeq_tools.py genomeAlign -f TDP120h1_Atrimmed.fastq

python ~/Dropbox/codes/endSeqTools/endSeq_tools.py samProcess --sam /Users/yasinkaymaz/Documents/data/PASseqData/TDP120h1_Atrimmed.sam

bamtools split -in TDP120h1_Atrimmed.bam -reference
samtools sort TDP120h1_Atrimmed.REF_chr22.bam -o TDP120h1_Atrimmed.REF_chr22.sorted.bam
samtools index TDP120h1_Atrimmed.REF_chr22.sorted.bam

python ~/Dropbox/codes/endSeqTools/endSeq_tools.py genomeCov --bam /Users/yasinkaymaz/Documents/data/PASseqData/TDP120h1_Atrimmed.REF_chr22.bam --re 3

#python ~/Dropbox/codes/endSeqTools/endSeq_tools.py NBclassifier --bed /Users/yasinkaymaz/Documents/data/PASseqData/TDP120h1_Atrimmed.REF_chr22_stranded_read_count.bed
Rscript --vanilla ~/Dropbox/codes/endSeqTools/src/cleanUpdTSeq.R /Users/yasinkaymaz/Documents/data/PASseqData/TDP120h1_Atrimmed.REF_chr22_stranded_read_count.bed nbpout.nbpred human

#python ~/Dropbox/codes/endSeqTools/endSeq_tools.py clusterPeaks --bed /Users/yasinkaymaz/Documents/data/PASseqData/TDP120h1_Atrimmed.REF_chr22_stranded_read_count.bed --NBprobs /Users/yasinkaymaz/Documents/data/PASseqData/nbpout.nbpred
python ~/Dropbox/codes/endSeqTools/endSeq_tools.py clusterPeaks --bed /Users/yasinkaymaz/Documents/data/PASseqData/TDP120h1_chr22_Atrimmed_sorted_stranded_read_count.bed --NBprobs /Users/yasinkaymaz/Documents/data/PASseqData/TDP120h1_chr22_Atrimmed_sorted_stranded_read_count.nbpred
python ~/Dropbox/codes/endSeqTools/endSeq_tools.py annotatePeaks --bed TDP120h1_chr22_Atrimmed_sorted_stranded_CPM_TruePeaks_GB.bed

#####

conda install -c bioconda STAR=2.7.1a



https://data.humancellatlas.org/explore/projects/c4077b3c-5c98-4d26-a614-246d12c2e5d7/get-curl-command

curl --location --fail 'https://service.azul.data.humancellatlas.org/manifest/files?catalog=dcp23&format=curl&filters=%7B%22fileFormat%22%3A+%7B%22is%22%3A+%5B%22fastq.gz%22%5D%7D%2C+%22projectId%22%3A+%7B%22is%22%3A+%5B%22c4077b3c-5c98-4d26-a614-246d12c2e5d7%22%5D%7D%2C+%22genusSpecies%22%3A+%7B%22is%22%3A+%5B%22Homo+sapiens%22%5D%7D%7D&objectKey=manifests%2F4355dc68-c197-5f3f-96b6-868d146f3a55.911e27d8-7978-5ffa-942b-c6dda54e8ae5.curlrc' | curl --config -

mkdir humanSTARindex\

STAR --runMode genomeGenerate \
--runThreadN 10 \
--genomeDir ./ \
--genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--sjdbGTFfile Homo_sapiens.GRCh38.99.gtf



STAR --genomeDir humanSTARindex \
--readFilesIn SRR8325947_reads/SRR8325947_2.fastq SRR8325947_reads/SRR8325947_1.fastq \
--runThreadN 16 \
--outSAMattributes CR CY UR UY \
--soloType Droplet \
--outSAMtype BAM Unsorted \
--soloCBwhitelist humanINDEX/737K-august-2016.txt \
--outFileNamePrefix SRR8325947_STAR/SRR8325947 > SRR8325947_STAR/STARsolo_SRR8325947.log 2>&1 | tee SRR8325947_STAR/STARsolo_SRR8325947.log



python ~/Dropbox/codes/endSeqTools/endSeq_tools.py samProcess \
--sam /Users/yasinkaymaz/Documents/data/PASseqData/SRR8325947_STAR/SRR8325947Aligned.out.sam

bamtools split -in SRR8325947Aligned.out_sorted.bam -reference

samtools index SRR8325947Aligned.out_sorted.REF_22.bam


python ~/Dropbox/codes/endSeqTools/endSeq_tools.py genomeCov --bam /Users/yasinkaymaz/Documents/data/PASseqData/SRR8325947_STAR/chrs/SRR8325947Aligned.out_sorted.REF_22.bam --re 3

#Added chr before 22 to make chr22
awk '{print "chr"$0}' /Users/yasinkaymaz/Documents/data/PASseqData/SRR8325947_STAR/chrs/SRR8325947Aligned.out_sorted.REF_22_stranded_read_count.bed > /Users/yasinkaymaz/Documents/data/PASseqData/SRR8325947_STAR/chrs/SRR8325947Aligned.out_sorted.REF_22_stranded_read_count.chr.bed
#python ~/Dropbox/codes/endSeqTools/endSeq_tools.py NBclassifier --bed /Users/yasinkaymaz/Documents/data/PASseqData/SRR8325947_STAR/chrs/SRR8325947Aligned.out_sorted.REF_22_stranded_read_count.bed
Rscript --vanilla ~/Dropbox/codes/endSeqTools/src/cleanUpdTSeq.R /Users/yasinkaymaz/Documents/data/PASseqData/SRR8325947_STAR/chrs/SRR8325947Aligned.out_sorted.REF_22_stranded_read_count.chr.bed nbpout.nbpred human

#python ~/Dropbox/codes/endSeqTools/endSeq_tools.py clusterPeaks --bed /Users/yasinkaymaz/Documents/data/PASseqData/SRR8325947_STAR/chrs/SRR8325947Aligned.out_sorted.REF_22_stranded_read_count.chr.bed --NBprobs /Users/yasinkaymaz/Documents/data/PASseqData/SRR8325947_STAR/chrs/nbpout.nbpred
~/Dropbox/codes/endSeqTools/src/cluster.sh 40 /Users/yasinkaymaz/Documents/data/PASseqData/SRR8325947_STAR/chrs/SRR8325947Aligned.out_sorted.REF_22_Atrimmed_sorted_stranded_read_count.bed /Users/yasinkaymaz/Documents/data/PASseqData/SRR8325947_STAR/chrs/nbpout.nbpred



tail -n+2 TDP120h1_chr22_Atrimmed_sorted_stranded_read_count.bed|\
clusterBed -s -d 40 -i stdin|\
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$1"_"$7}'|\
groupBy -g 7 -c 1,2,2,3,3,4,5,5,5,5,5,5,5,6 -o collapse,collapse,min,collapse,max,collapse,collapse,max,min,sum,median,mean,stdev,collapse|\
expandCols -c 2,3,5,7,8,15 |\
awk '{if($8 == $9)printf ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.02f\t%.02f\t%s\n", $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15)}'|\
awk '{print $2"\t"$4"\t"$6"\t"$7"\t"$8"\t"$15"\t"$3"\t"$5"\t"$1"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$6-$4}'|\
sort -k1,1 -k2,2g |\
groupBy -g 9 -c 1,2,3,4,5,6,7,8,10,11,12,13,14,15,16 -o first,first,first,first,first,first,first,last,first,first,first,first,first,first,first|\
awk '{print $2"\t"$3"\t"$4"\t"$1"\t"$6"\t"$7"\t"$8"\t"$9"\t"$5"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16}' > tmp1.bed

cut -f1,7 data8325947_metaClusters.txt > my_barcodes.txt


samtools sort SRR8325947Aligned.out.bam -o SRR8325947Aligned.sorted.bam -n


samtools view -h SRR8325947Aligned.sorted.bam |\
sed -e '/^@SQ/s/SN\:/SN\:chr/' -e '/^[^@]/s/\t/\tchr/2'|\
awk -F ' ' '$7=($7=="=" || $7=="*"?$7:sprintf("chr%s",$7))' |\
tr " " "\t"|less

samtools view -h SRR8325947Aligned.sorted.bam |\
sed -e '/^@SQ/s/SN\:/SN\:chr/' -e '/^[^@]/s/\t/\tchr/2'|\
awk -F '\t' 'BEGIN{OFS="\t";} { if ($1 !~ /^@/ && $7 != "=" && $7 != "*") {$7 = "chr"$7; print$0;} else print$0}'|less


samtools view -h SRR8325947Aligned.sorted.bam |sed -e '/^@SQ/s/SN\:/SN\:chr/' -e '/^[^@]/s/\t/\tchr/2'|awk -F '\t' 'BEGIN{OFS="\t";} { if ($1 !~ /^@/ && $7 != "=" && $7 != "*") {$7 = "chr"$7; print$0;} else print$0}'|samtools view - -b -o SRR8325947Aligned.sorted.chr.bam

samtools sort SRR8325947Aligned.sorted.chr.bam -o SRR8325947Aligned.sorted.bam
samtools index SRR8325947Aligned.sorted.bam

sinto filterbarcodes -b SRR8325947Aligned.sorted.bam -c ../my_barcodes.txt --barcodetag "CR"

samtools index 0.bam



######
bamtools merge -list bams.list -out merged.bam

samtools sort merged.bam -o merged_sorted.bam

bamtools split -in merged_sorted.bam -reference

samtools index merged_sorted.REF_1.bam
python ~/Dropbox/codes/endSeqTools/endSeq_tools.py genomeCov --bam /Users/yasinkaymaz/Documents/data/PASseqData/SRR832_backup/merged_sorted.REF_1.bam --re 3

awk '{print "chr"$0}' merged_sorted.REF_1_stranded_read_count.bed > merged_sorted.REF_1_Atrimmed_sorted_stranded_read_count.bed

#Rscript --vanilla ~/Dropbox/codes/endSeqTools/src/cleanUpdTSeq.R /Users/yasinkaymaz/Documents/data/PASseqData/SRR832_backup/merged_sorted.REF_1_Atrimmed_sorted_stranded_read_count.bed nbpout.nbpred human

Rscript --vanilla ~/Dropbox/codes/endSeqTools/src/SetNBpred.1.0.R \
/Users/yasinkaymaz/Documents/data/PASseqData/SRR832_backup/merged_sorted.REF_1_Atrimmed_sorted_stranded_read_count.bed \
/Users/yasinkaymaz/Documents/data/PASseqData/SRR832_backup/merged_sorted.REF_1_Atrimmed_sorted_stranded_read_count.nbpred human


~/Dropbox/codes/endSeqTools/src/cluster.sh 40 \
/Users/yasinkaymaz/Documents/data/PASseqData/SRR832_backup/merged_sorted.REF_1_Atrimmed_sorted_stranded_read_count.bed \
/Users/yasinkaymaz/Documents/data/PASseqData/SRR832_backup/merged_sorted.REF_1_Atrimmed_sorted_stranded_read_count.nbpred


conda activate MedGen
snakemake --cores 2 -s SnakeFile.smk --reason
snakemake --cores 2 -s SplitBams.smk --reason
snakemake --cores 2 -s ClustersPAS.smk --reason

ls -1 *_STAR/*_*_Atrimmed_sorted_stranded_read_count.nbpred|sed 's/_Atrimmed_sorted_stranded_read_count.nbpred//g'|awk '{print $0".fastq"}' > master_samples.list.txt

conda activate CAAPS


key="Tumor1A_vs_Healthy_Endothelial"
key="Brain_met_IV_vs_TumorIIA_Cancer"


bash ~/Dropbox/codes/CAAPS/scripts/make_PAS_table.v3.sh \
"${key}"_samples.txt \
humanSTARindex/Homo_sapiens.GRCh38.99.gtf  \
~/Dropbox/codes/CAAPS/scripts/

Rscript ~/Dropbox/codes/CAAPS/scripts/PASseq_chisquire_beta_v4.R \
mergedPeaks_annotated_CPM-3.bed \
"${key}"_APA_switches.txt \
3 2 \
/Users/yasinkaymaz/Documents/data/PASseqData
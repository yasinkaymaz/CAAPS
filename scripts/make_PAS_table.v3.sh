#!/bin/bash


WORKINGDIR=`dirname $1`
echo "Working directory is: $WORKINGDIR . Please check the output here!"
SamplesFile=`basename $1`
echo "Major Samples File is: '$SamplesFile'"
#This file is the gene annotation file in GTF format from which gene boundaries will be extracted.
GeneAnnotationFile=$2
echo "Gene annotation file which will be used: $GeneAnnotationFile"

scriptDir=$3


GeneSpaceFile=${GeneAnnotationFile%.gtf}.gsf
GenesFile=${GeneAnnotationFile%.gtf}.ann

if [ -f $WORKINGDIR/Samples_stats.txt ];
	then
	rm $WORKINGDIR/Samples_stats.txt
else
	echo "Creating a new Samples_stats file!"
fi

if [ -f $GeneSpaceFile ];
	then
	echo "File $GeneSpaceFile exists. No need to generate gene space file (.gsf) again!"
#$scriptDir/GTF_parser_For_GeneAnnotation.sh $GeneAnnotationFile $GenesFile
else
	echo "File $GeneSpaceFile does not exist. Please make sure that GTF file includes gene entries! Processing GTF file to create gene space file (.gsf)..."
$scriptDir/GTF_parser_For_GeneSpace.sh $GeneAnnotationFile $GeneSpaceFile
$scriptDir/GTF_parser_For_GeneAnnotation.sh $GeneAnnotationFile $GenesFile
   echo "Creating .gsf file is complete!"
fi

FileNametoMerge='Atrimmed_sorted_stranded_read_count.bed'
ReadCountFile='Atrimmed_sorted_stranded_read_count.bed'
NBpredFile='Atrimmed_sorted_stranded_read_count.nbpred'

echo "Files that will be used are *_$ReadCountFile and *_$NBpredFile of each sample stored with unique names in $SamplesFile"

while read line
do
Fastq=$line
basename=${Fastq%.fastq}

if [ -f "$basename"_genespace_lambda.bed.tmp ];
then
   echo "File "$basename"_genespace_lambda.bed.tmp exists. No need to generate again!"
else

echo "Combining read count and nbpred files for sample: $basename ..."
#Note: You don't need to calculate lambda for individual samples. But 
#tail -n+2|sed 's/\"//g'|\

cut -f2 "$basename"_$NBpredFile |\
sed 's/\"//g'|\
paste -d "\t" "$basename"_$ReadCountFile - |\
bedtools intersect -s -wa -wb -a - -b $GeneSpaceFile |\
sort -k11,11 |\
groupBy -g 11 -c 1,2,3,4,5,6,7,8,9,10,11,5 -o collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse,mean |\
expandCols -c 2,3,4,5,6,7,8,9,10,11,12 |\
cut -f 2-8,12,13 > "$basename"_genespace_lambda.bed.tmp

echo "Pooling the candidate PA sites..."
cat "$basename"_genespace_lambda.bed.tmp >> $WORKINGDIR/PApool.bed.tmp
echo "Done for $basename! Next..."
#*############################################################################################################################
cut -f2 "$basename"_$NBpredFile |sed 's/\"//g'|paste -d "\t" "$basename"_$ReadCountFile - |\
bedtools intersect -s -wa -v -a - -b $GeneAnnotationFile |\
bedtools intersect -s -wa -v -a - -b $GeneSpaceFile |\
awk '($7 < 0.001)'|cut -f1-6 |sort -k1,1 -k2,2g |uniq|clusterBed -s -d 40 |awk '{OFS="\t"; $(NF)=$1"_"$(NF); print $0}' |\
groupBy -g 7 -c 1,2,3,5,6,5 -o first,min,max,collapse,first,sum -prec 100|cut -f2-7|\
awk '{ print $1"\t"$2"\t"$3"\t""Intergenic.PAs"NR"\t"$6"\t"$5 }'| awk '($5 > 4)' > "$basename"_intergenic.PAs.bed.tmp

cat "$basename"_intergenic.PAs.bed.tmp >> $WORKINGDIR/PApool_intergenic.bed.tmp
#***##########################################################################################################################


fi

done < $WORKINGDIR/$SamplesFile
echo "Combining files is done for all samples!"

if [ -f $WORKINGDIR/PApool.bed_sorted_lambda_ppois_bed.tmp ];
then
   echo "File PApool.bed_sorted_lambda_ppois_bed.tmp exists. No need to generate again!"
else

echo "Sorting the Pooled PA data..."
sort -k1,1 -k2,2g $WORKINGDIR/PApool.bed.tmp > $WORKINGDIR/PApool.bed_sorted.tmp
#*############################################################################################################################
sort -k1,1 -k2,2g $WORKINGDIR/PApool_intergenic.bed.tmp > $WORKINGDIR/PApool_intergenic_sorted.bed.tmp
#***##########################################################################################################################
#-S 50% --parallel=6
echo "Now calculating the lambda parameter of each gene for combined noise removal..."
awk '{print $1"_"$2"_"$3"_"$6"\t"$5"\t"$7 }' $WORKINGDIR/PApool.bed_sorted.tmp |\
groupBy -g 1 -c 2,3 -o sum,min |\
awk -F _ '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}'|\
awk '{OFS="\t"; print $1,$2,$3,"poolID",$5,$4,$6}'|\
bedtools intersect -s -wa -wb -a - -b $GeneSpaceFile |\
sort -k11,11 |\
groupBy -g 11 -c 1,2,3,4,5,6,7,8,9,10,11,5 -o collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse,mean |\
expandCols -c 2,3,4,5,6,7,8,9,10,11,12 |\
cut -f 2-8,12,13 > $WORKINGDIR/PApool.bed_sorted_lambda.bed.tmp
echo "Lambda values have been calculated!"
#-S 50% --parallel=6
echo "Calculating a pvalue for a each PA site using Poisson distribution... This pvalue will be used for noise filtration..."
echo "Genomic locations that are likely Internal priming events are also filtered!"
Rscript $scriptDir/ppois.R $WORKINGDIR/PApool.bed_sorted_lambda.bed.tmp $WORKINGDIR/PApool.bed_sorted_lambda_ppois.bed.tmp
echo "Filtration is complete on Pooled data!"

echo "Now clustering closer genomic locations within 100nt to merge expression values..."
echo "PA regions of which cumulative read counts do not exceed to 100 will be excluded..."
#Column 7 is the PA sites read count total of each gene after clustering closer locations. Exclude low count PA sites
cut -f1-6 $WORKINGDIR/PApool.bed_sorted_lambda_ppois.bed.tmp |\
sort -k2,2g |uniq|\
clusterBed -s -d 40 |\
awk '{OFS="\t"; $(NF)=$1"_"$(NF); print $0}' |\
groupBy -g 7 -c 1,2,3,5,6,5 -o collapse,collapse,collapse,collapse,collapse,sum -prec 100|\
awk '($7 > 100)' |\
expandCols -c 2,3,4,5,6 |\
cut -f2-6 |\
awk '{ print $1"_"$2"_"$3"_"$5"\t"$4 }' > $WORKINGDIR/PApool.bed_sorted_lambda_ppois_bed.tmp
#*############################################################################################################################
awk '{print $1"_"$2"_"$3"_"$6"\t"$5 }' $WORKINGDIR/PApool_intergenic_sorted.bed.tmp > $WORKINGDIR/PApool_intergenic_sorted.bed.tmp.tmp
#***##########################################################################################################################

echo "Low expressed PA regions excluded!"
#-S 50% --parallel=6 -k1,1


fi


#Choose the value type in matrix. For each value type: 
array=( Count Percent CPM )

for valueType in "${array[@]}"
do
echo "Processing files for: $valueType"

echo $WORKINGDIR/PApool.bed_sorted_lambda_ppois_bed.tmp >> $WORKINGDIR/samples.tmp
#*############################################################################################################################
#################echo $WORKINGDIR/PApool_intergenic_sorted.bed.tmp.tmp >> $WORKINGDIR/samples_intergenic.tmp
#***##########################################################################################################################

while read line
do
Fastq=$line
echo "Sample file: $line"

basename=${Fastq%.fastq}
awk '{print $1"_"$2"_"$3"_"$6"\t"$5 }' "$basename"_$FileNametoMerge > "$basename"_"$FileNametoMerge".tmp

#Calculate the library size
total_c=`awk 'FNR==NR{a[$1]=$0;next}{ print $0"\t", a[$1]}' "$basename"_"$FileNametoMerge".tmp $WORKINGDIR/PApool.bed_sorted_lambda_ppois_bed.tmp |\
awk '{if(NF >3) print $3"\t"$4}'|\
awk -F _ '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}'|\
awk '{x=x+$5}END{print x}'`

if [ $valueType = "Count" ]
then
factor='$5*1'
valueType='Count'
echo "Now calculating the raw read counts of PAs for sample $basename"

elif [ $valueType = "Percent" ]
then
valueType='Percent'
factor='$5*100/$7'
echo "Now calculating the percent PA usages for sample $basename"

else
valueType='CPM'
factor='$5'"*1000000/$total_c"
echo "Now calculating the normalized PA expressions for sample $basename"
echo "Library size for $basename  is: $total_c" >> $WORKINGDIR/Samples_stats.txt

fi

echo "Now checking the sites that pass the Noise filtration and Internal Priming event filtration for sample $basename ..."
FilterLowCPM=1.0
awk 'FNR==NR{a[$1]=$0;next}{ print $0"\t", a[$1]}' "$basename"_"$FileNametoMerge".tmp $WORKINGDIR/PApool.bed_sorted_lambda_ppois_bed.tmp |\
awk '{if(NF >3) print $3"\t"$4}'|\
uniq |\
awk -F _ '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}'|\
awk '{OFS="\t"; print $1,$2,$3,"poolID",$5,$4}'|\
sort -k1,1 -k2,2g |\
clusterBed -s -d 40 |\
awk '{OFS="\t"; $(NF)=$1"_"$(NF); print $0}' |\
groupBy -g 7 -c 1,2,3,4,5,6,5 -o collapse,collapse,collapse,collapse,collapse,collapse,sum -prec 100|\
awk '{if($8*1000000/'$total_c' > '$FilterLowCPM') print $0}' |\
expandCols -c 2,3,4,5,6,7 |\
cut -f2-7 |\
bedtools intersect -s -wa -wb -a - -b $GeneSpaceFile |\
sort -k10,10|\
groupBy -g 10 -c 1,2,3,4,5,6,5 -o collapse,collapse,collapse,collapse,collapse,collapse,sum -prec 100|\
expandCols -c 2,3,4,5,6,7|\
cut -f2-|\
awk '{print $0"\t"'$factor'}'|\
awk '{OFS="\t"; print $1,$2,$3,$4,$8,$6}'|\
awk '{print $1"_"$2"_"$3"_"$6"\t"$5 }' > "$basename"_"$FileNametoMerge".filtered.tmp
#-S 50% --parallel=6
#*############################################################################################################################
awk 'FNR==NR{a[$1]=$0;next}{ print $0"\t", a[$1]}' "$basename"_"$FileNametoMerge".tmp $WORKINGDIR/PApool_intergenic_sorted.bed.tmp.tmp |\
awk '{if(NF >3) print $3"\t"$4}'|\
uniq |\
awk -F _ '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}'|\
awk '{OFS="\t"; print $1,$2,$3,"poolID",$5,$4}'|\
sort -k1,1 -k2,2g |\
clusterBed -s -d 40 |\
awk '{OFS="\t"; $(NF)=$1"_"$(NF); print $0}' |\
groupBy -g 7 -c 1,2,3,4,5,6,5 -o collapse,collapse,collapse,collapse,collapse,collapse,sum -prec 100|\
awk '{if($8*1000000/'$total_c' > '$FilterLowCPM') print $0}' |\
expandCols -c 2,3,4,5,6,7 |\
cut -f2-7 |\
awk '{print $1"_"$2"_"$3"_"$6"\t"$5 }' > "$basename"_"$FileNametoMerge".intergenic.tmp
#***##########################################################################################################################

echo "Filtration and Normalization steps are done!"
#awk '{printf ("%s\t%s\t%s\t%s\t%.02f\t%s\n", $1,$2,$3,$4,$8,$6)}'|\
#awk '{if($8*1000000/'$total_c' > '$FilterLowCPM') printf ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.02f\n", $1,$2,$3,$4,$5,$6,$7,$8)}' |\
echo "$basename"_"$FileNametoMerge".filtered.tmp >> $WORKINGDIR/samples.tmp
#*############################################################################################################################
echo "$basename"_"$FileNametoMerge".intergenic.tmp >> $WORKINGDIR/samples_intergenic.tmp
#***##########################################################################################################################

echo "Processing individual sample $basename is done! Starting to next sample..."

done < $WORKINGDIR/$SamplesFile



x=`wc -l $WORKINGDIR/samples.tmp|awk '{print $1}'`
let NoS=$x-1
echo "Number of samples:" $NoS

c=$(for a in `seq $x`;do let b=$a+6; printf ",""%s"$b ; done)
o=$(for a in `seq $x`;do printf ",""sum" ; done)

let g=$x+7

$scriptDir/merge_bed_files.pl $WORKINGDIR/samples.tmp |\
tail -n+2 |\
awk -F _ '{print $1"\t"$2"\t"$3"\t""ID""\t""1""\t"$4}'|\
sort -k1,1 -k2,2g |\
clusterBed -s -d 40 |\
awk '{OFS="\t"; $(NF)=$1"_"$(NF); print $0}' |\
groupBy -g $g -c 1,2,3,4,5,6$c -o first,min,max,concat,sum,first$o -prec 100|\
awk '{OFS="\t"; $5=$1; print $0}' |\
cut -f2- |uniq > $WORKINGDIR/mergedPeaks_"$valueType".bed

#*############################################################################################################################
#Because we are not including PApool_intergenic.bed file in samples_intergenic.tmp file:
x=`wc -l $WORKINGDIR/samples_intergenic.tmp|awk '{print $1}'`
let NoS=$x-1
echo "Number of samples:" $NoS

c=$(for a in `seq $x`;do let b=$a+6; printf ",""%s"$b ; done)
o=$(for a in `seq $x`;do printf ",""sum" ; done)

let g=$x+7

$scriptDir/merge_bed_files.pl $WORKINGDIR/samples_intergenic.tmp |\
tail -n+2 |\
awk -F _ '{print $1"\t"$2"\t"$3"\t""ID""\t""1""\t"$4}'|\
sort -k1,1 -k2,2g |\
clusterBed -s -d 40 |\
awk '{OFS="\t"; $(NF)=$1"_"$(NF); print $0}' |\
groupBy -g $g -c 1,2,3,4,5,6$c -o first,min,max,concat,sum,first$o -prec 100|\
awk '{OFS="\t"; $5=$1; print $0}' |\
cut -f2- |sort -k1,1 -k2,2g |uniq > $WORKINGDIR/mergedPeaks_of_IntergenicPAs_"$valueType".bed
#***##########################################################################################################################
x=`wc -l $WORKINGDIR/samples.tmp|awk '{print $1}'`
let NoS=$x-1
echo "Number of samples:" $NoS

c=$(for a in `seq $x`;do let b=$a+6; printf ",""%s"$b ; done)
o=$(for a in `seq $x`;do printf ",""sum" ; done)

let g=$x+7

#-S 50% --parallel=6
let n=$x+6
let m=$x+6+4
let t=$x+6+7
let k=$x+6+11
#
#***
module unload bedtools/2.22.0
module load bedtools/2.17.0

awk '($3 == "polyA_site")' /project/umw_jeffrey_bailey/share/Homo_sapiens/Gencode/v19/gencode.v19.polyAs.gtf > $WORKINGDIR/PAsites.gtf.tmp
bedtools intersect -s -wa -wb -a $WORKINGDIR/mergedPeaks_"$valueType".bed -b $GenesFile |uniq |bedtools window -w 20 -sm -a - -b $WORKINGDIR/PAsites.gtf.tmp > $WORKINGDIR/mergedPeaks_annotated_"$valueType".KnownPAS.bed.tmp
bedtools intersect -s -wa -wb -a $WORKINGDIR/mergedPeaks_"$valueType".bed -b $GenesFile |uniq |bedtools window -v -w 20 -sm -a - -b $WORKINGDIR/PAsites.gtf.tmp > $WORKINGDIR/mergedPeaks_annotated_"$valueType".NovelPAS.bed.tmp
cat $WORKINGDIR/mergedPeaks_annotated_"$valueType".KnownPAS.bed.tmp $WORKINGDIR/mergedPeaks_annotated_"$valueType".NovelPAS.bed.tmp| sort -k1,1 -k2,2g|cut -f1-$n,$m,$t,$k|uniq > $WORKINGDIR/mergedPeaks_annotated_"$valueType".bed

#***

#bedtools intersect -s -wa -wb -a $WORKINGDIR/mergedPeaks_"$valueType".bed -b $GeneSpaceFile |cut -f1-$n,$m |uniq > $WORKINGDIR/mergedPeaks_annotated_"$valueType".bed

#awk '{$('$x'+7)=$1":"$2"-"$3"_"$('$x'+7); OFS="\t"; print $0}' mergedPeaks_annotated_"$valueType".bed |cut -f8-|awk '{print $('$x')"\t"$0}'|cut -f1-$x |uniq > heatmap_annotated_"$valueType".txt

module unload bedtools/2.17.0
module load bedtools/2.22.0

rm $WORKINGDIR/samples.tmp
rm $WORKINGDIR/samples_intergenic.tmp

c=$(for a in `seq $NoS`;do let b=$a+7; printf ",""%s"$b ; done)
o=$(for a in `seq $NoS`;do printf ",""sum" ; done)
sort -k$g,$g $WORKINGDIR/mergedPeaks_annotated_Count.bed |\
groupBy -g $g -c 1,2,3,$g,7,6$c -o first,min,max,first,sum,first$o -prec 100|\
cut -f2- |uniq > $WORKINGDIR/mergedPeaks_annotated_GeneCounts.bed


#Rscript /project/umw_jeffrey_bailey/OTHERS/endSeq_Tools/src/postMergingProcess.R heatmap_annotated_CPM.txt heatmap_annotated_CPM_Filtered.txt $x


done

echo "Done"


module unload bedtools/2.22.0
module load bedtools/2.17.0
bedtools intersect -s -wa -wb -a /project/umw_jeffrey_bailey/share/Homo_sapiens/UCSC/hg19/utr.bed -b  $WORKINGDIR/mergedPeaks_CPM.bed |sort -u -k10,10 > $WORKINGDIR/3utr_CPM.bed 

awk '(index($14,"-") == 0){ print }' $WORKINGDIR/mergedPeaks_annotated_CPM.bed > $WORKINGDIR/mergedPeaks_annotated_CPM-2.bed
cut -f4 $WORKINGDIR/mergedPeaks_annotated_CPM-2.bed|sort |uniq -c|sort -k1,1gr|awk '{if($1 >1)print $2}' > $WORKINGDIR/.tmp
awk 'NR==FNR{n[$0];next} !($4 in n){print}' $WORKINGDIR/.tmp $WORKINGDIR/mergedPeaks_annotated_CPM-2.bed > $WORKINGDIR/mergedPeaks_annotated_CPM-3.bed

#This is how to run this script
#/project/umw_jeffrey_bailey/OTHERS/endSeq_Tools/src/make_samples_matrix.sh /project/umw_jeffrey_bailey/yk42w/results/PASseq/Wald/RNaseH/sortedNbpreds-beds/Test/file_contains_sample_fastqs.txt /project/umw_jeffrey_bailey/share/Homo_sapiens/Gencode/v19/gencode.v19.annotation.gtf



rm $WORKINGDIR/*.tmp
while read line
do
Fastq=$line
basename=${Fastq%.fastq}
rm "$basename"*.tmp
done < $WORKINGDIR/$SamplesFile



#CHANGE LOG:
#Line 93 first changed to min. see below.
#groupBy -g 1 -c 2,3 -o sum,min |\


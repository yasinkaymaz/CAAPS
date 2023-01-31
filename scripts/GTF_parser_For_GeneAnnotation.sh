#!/bin/bash

GTFfile=$1
GeneSpaceOutFile=$2
#Plus strand
#	if($3=="gene" && $18=="\"protein_coding\";" && $7=="+"){
awk '{
	if($3=="gene" && $7=="+" && ($18=="\"protein_coding\";" || $18=="\"snRNA\";" || $18=="\"lincRNA\";" || $18=="\"antisense\";" || $18=="\"miRNA\";" || $18=="\"processed_transcript\";") ){
		print $0
	}
}' $GTFfile |\
awk '{
	OFS="\t";
	print $1,$4,$5,$7,$14,$18,$5
}'|\
sort -k1,1 -k2,2gr |\
sed '
	s/[";]//g
'|\
awk '{	i=1;
	chr[NR]=$1;
	str[NR]=$2;
	end[NR]=$3;
	snd[NR]=$4;
	gen[NR]=$5;
	type[NR]=$6;
	TranscriptEnd[NR]=$7;
	if (NR >1){
		chr[NR]=$1;
		str[NR]=$2;
		end[NR]=$3;
		snd[NR]=$4;
		gen[NR]=$5;
		type[NR]=$6;
		TranscriptEnd[NR]=$7;
		while (end[NR] >= str[NR-i]){
				i++
		}
		if(chr[NR] == chr[NR-i]){
			maxDist=str[NR-i]-end[NR]
			if(maxDist > 5000){
				print chr[NR]"\t"str[NR]"\t"end[NR]+5000"\t"gen[NR]"\t"end[NR]-str[NR]+5000"\t"snd[NR]"\t"type[NR]"\t"TranscriptEnd[NR]
				i=1;						
			}
			else
				print chr[NR]"\t"str[NR]"\t"str[NR-i]"\t"gen[NR]"\t"str[NR-i]-str[NR]"\t"snd[NR]"\t"type[NR]"\t"TranscriptEnd[NR]
				i=1;
		}
	}
	else
		print chr[NR]"\t"str[NR]"\t"end[NR]+2000"\t"gen[NR]"\t"end[NR]+2000-str[NR]"\t"snd[NR]"\t"type[NR]"\t"TranscriptEnd[NR]
}'| awk '{print "chr"$0}' > plus.tmp


#minus strand
#	if($3=="gene" && $18=="\"protein_coding\";" && $7=="-"){
awk '{
	if($3=="gene" && $7=="-" && ($18=="\"protein_coding\";" || $18=="\"snRNA\";" || $18=="\"lincRNA\";" || $18=="\"antisense\";" || $18=="\"miRNA\";" || $18=="\"processed_transcript\";") ){
		print $0
	}
}' $GTFfile |\
awk '{
	OFS="\t";
	print $1,$4,$5,$7,$14,$18,$4
}'|\
sort -k1,1 -k2,2g |\
sed '
	s/[";]//g
'|\
awk '{	i=1;
	chr[NR]=$1;
	str[NR]=$2;
	end[NR]=$3;
	snd[NR]=$4;
	gen[NR]=$5;
	type[NR]=$6;
	TranscriptEnd[NR]=$7;
	if (NR >1){
		chr[NR]=$1;
		str[NR]=$2;
		end[NR]=$3;
		snd[NR]=$4;
		gen[NR]=$5;
		type[NR]=$6;
		TranscriptEnd[NR]=$7;
		while (end[NR-i] >= str[NR]){
				i++
		}
		if(chr[NR] == chr[NR-i]){
                        maxDist=end[NR]-end[NR-i]
                        if(maxDist > 5000){
				print chr[NR]"\t"str[NR]-5000"\t"end[NR]"\t"gen[NR]"\t"end[NR]-str[NR]+5000"\t"snd[NR]"\t"type[NR]"\t"TranscriptEnd[NR]
				i=1;
			}
			else
				print chr[NR]"\t"end[NR-i]"\t"end[NR]"\t"gen[NR]"\t"end[NR]-end[NR-i]"\t"snd[NR]"\t"type[NR]"\t"TranscriptEnd[NR]
				i=1;
		}
	}
	else
		print chr[NR]"\t"str[NR]-2000"\t"end[NR]"\t"gen[NR]"\t"end[NR]-str[NR]-2000"\t"snd[NR]"\t"type[NR]"\t"TranscriptEnd[NR]
}'| awk '{print "chr"$0}' > minus.tmp

cat plus.tmp minus.tmp |\
sort -k1,1 -k2,2g > $GeneSpaceOutFile
rm *.tmp



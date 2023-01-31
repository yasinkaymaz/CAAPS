#!/bin/bash

GTFfile=$1
GeneSpaceOutFile=$2
#Plus strand
#	if($3=="gene" && $7=="+" && ($22 == "\"KNOWN\";") && ($20=="\"protein_coding\";" || $20=="\"snRNA\";" || $20=="\"lincRNA\";" || $20=="\"antisense\";" || $20=="\"miRNA\";" || $20=="\"processed_transcript\";") ){
awk '{
	if($3=="gene" && $18=="\"protein_coding\";" && $7=="+"){
		print $0
	}
}' $GTFfile |\
awk '{
	OFS="\t";
	print $1,$4,$5,$7,$14
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
	if (NR >1){
		chr[NR]=$1;
		str[NR]=$2;
		end[NR]=$3;
		snd[NR]=$4;
		gen[NR]=$5;
		while (end[NR] >= str[NR-i]){
				i++
		}
		if(chr[NR] == chr[NR-i]){
			maxDist=str[NR-i]-end[NR]
			if(maxDist > 5000){
	                        print chr[NR]"\t"str[NR]"\t"end[NR]+5000"\t"gen[NR]"\t"end[NR]-str[NR]+5000"\t"snd[NR]"\t"
				i=1;
			}
			else
				print chr[NR]"\t"str[NR]"\t"str[NR-i]"\t"gen[NR]"\t"str[NR-i]-str[NR]"\t"snd[NR]"\t"
				i=1;
		}
	}
	else
		print chr[NR]"\t"str[NR]"\t"end[NR]+2000"\t"gen[NR]"\t"end[NR]+2000-str[NR]"\t"snd[NR]
}' | awk '{print "chr"$0}' > plus.tmp


#minus strand
#	if($3=="gene" && $7=="-" && ($22 == "\"KNOWN\";") && ($20=="\"protein_coding\";" || $20=="\"snRNA\";" || $20=="\"lincRNA\";" || $20=="\"antisense\";" || $20=="\"miRNA\";" || $20=="\"processed_transcript\";") ){
awk '{
	if($3=="gene" && $18=="\"protein_coding\";" && $7=="-"){
		print $0
	}
}' $GTFfile |\
awk '{
	OFS="\t";
	print $1,$4,$5,$7,$14
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
	if (NR >1){
		chr[NR]=$1;
		str[NR]=$2;
		end[NR]=$3;
		snd[NR]=$4;
		gen[NR]=$5;
		while (end[NR-i] >= str[NR]){
				i++
		}
		if(chr[NR] == chr[NR-i]){
                        maxDist=end[NR]-end[NR-i]
                        if(maxDist > 5000){
				print chr[NR]"\t"str[NR]-5000"\t"end[NR]"\t"gen[NR]"\t"end[NR]-str[NR]+5000"\t"snd[NR]
				i=1;
			}
			else
				print chr[NR]"\t"end[NR-i]"\t"end[NR]"\t"gen[NR]"\t"end[NR]-end[NR-i]"\t"snd[NR]
				i=1;
		}
	}
	else
		print chr[NR]"\t"str[NR]-2000"\t"end[NR]"\t"gen[NR]"\t"end[NR]-str[NR]-2000"\t"snd[NR]
}'| awk '{print "chr"$0}' > minus.tmp

cat plus.tmp minus.tmp |\
sort -k1,1 -k2,2g |cut -f1-6 > $GeneSpaceOutFile
rm *.tmp



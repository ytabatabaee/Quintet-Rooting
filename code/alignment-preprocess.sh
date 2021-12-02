#!/bin/bash
while getopts f:n:r: flag
do
    case "${flag}" in
        f) base_path=${OPTARG};;
        n) taxa_count=${OPTARG};;
        r) replicates=${OPTARG};;
    esac
done


echo "Path: $base_path";
echo "Total species: $taxa_count";
echo "Replicates: $replicates";

taxon_name=0

for replicate in $(seq 1 $replicates)
do
	echo "Replicate $replicate"
	
	concatenated="$base_path/R$replicate/R$replicate.concat.fasta"
	final="$base_path/R$replicate/R$replicate.fasta"
	
	singlelist=""
	
	for taxon in $(seq 1 $taxa_count)
	do
		let "taxon_name++"
		echo "Taxon $taxon_name"
		
		compressed="$base_path/R$replicate/$taxon_name/$taxon_name.fasta.gz"
		decompressed="$base_path/R$replicate/$taxon_name/$taxon_name.fasta"
		sorted="$base_path/R$replicate/$taxon_name/$taxon_name.sorted.fasta"
		singlelined="$base_path/R$replicate/$taxon_name/$taxon_name.singleline.fasta"
		
		unzip="gzip -dk $compressed"
		eval $unzip
		sort="seqkit sort $decompressed > $sorted"
		eval $sort
		
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' "$sorted" > "$singlelined"
		singlelist="$singlelist $singlelined"
	done
	concat="seqkit concat $singlelist > $concatenated"
	eval $concat
	awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' "$concatenated" > "$final"
done

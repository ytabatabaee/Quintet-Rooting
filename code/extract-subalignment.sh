#!/bin/bash
while getopts i:n:c:o: flag
do
    case "${flag}" in
        i) input_fasta=${OPTARG};;
        n) subset_taxa=${OPTARG};;
        c) test_cases=${OPTARG};;
        o) output_folder=${OPTARG};;
    esac
done

for case in $(seq 1 $test_cases)
do
	case_folder="$output_folder/$case"
	mkdir $case_folder
	temp="$case_folder/$case.multiline.fasta"
	output_fasta="$case_folder/$case.fasta"
	subset_label="$case_folder/$case.txt"
	seqkit shuffle -s "$RANDOM" -2 "$input_fasta"|seqkit head -n "$subset_taxa" > "$temp"
	awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' "$temp" > "$output_fasta"
	awk 'sub(/^>/, "")' "$output_fasta" > "$subset_label"
done


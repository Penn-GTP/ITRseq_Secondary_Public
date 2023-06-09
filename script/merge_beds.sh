#!/bin/bash

pattern=$1
outname=$2

cat *${pattern}*ref_peak_track.bed > ${outname}_cat.bed #concatenate the files specified

cat_num_lines=$(cat ${outname}_cat.bed | wc -l)

bedtools sort -i ${outname}_cat.bed > ${outname}_cat_sorted.bed #sort the concatenated bed file
rm ${outname}_cat.bed

bedtools merge -i ${outname}_cat_sorted.bed -d 44 -c 4,5,6 -o collapse,sum,distinct > ${outname}_merged.bed #merge the regions in the bed files
rm ${outname}_cat_sorted.bed

merged_num_lines=$(cat ${outname}_merged.bed | wc -l)

sort -nrk 5 ${outname}_merged.bed > ${outname}_merged_sorted.bed #sort by score
rm ${outname}_merged.bed

bedtools window -w 500 -a ${outname}_merged_sorted.bed -b ARCUS_target_Mmul10.bed -v > ${outname}_merged_sorted_offtarget.bed #the inverse of the step below, only keeping off-target sites

bedtools window -w 500 -a ${outname}_merged_sorted.bed -b ARCUS_target_Mmul10.bed -u > ${outname}_merged_sorted_ontarget.bed #finds any site within 500bp of the on-target site and puts it in a separate file

#awk '{if($3-$2 >= 120 && $3-$2 <= 200) print}' merged_sorted.bed > merged_sorted_filtered.bed #filter by min and max length
rm ${outname}_merged_sorted.bed

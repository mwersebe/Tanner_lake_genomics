#!/bin/bash

# Make A file with locus name and Allele "1" counts for each population:

bcftools view -S TOP_samples.txt Tanners_CHRall_filter_BiSnps_LDprune_nomiss.vcf.gz|bcftools +fill-tags|bcftools view -H|while read line; do echo $line|awk -F "\ " '{print$1"_"$2}'|tr "\n" "\t";echo $line|awk -F ";" '{print$2}'|sed -e "s/AC=//g"; done > TOP_allele_counts.txt

bcftools view -S MID_samples.txt Tanners_CHRall_filter_BiSnps_LDprune_nomiss.vcf.gz|bcftools +fill-tags|bcftools view -H|while read line; do echo $line|awk -F "\ " '{print$1"_"$2}'|tr "\n" "\t";echo $line|awk -F ";" '{print$2}'|sed -e "s/AC=//g"; done > MID_allele_counts.txt

bcftools view -S DEEP_samples.txt Tanners_CHRall_filter_BiSnps_LDprune_nomiss.vcf.gz|bcftools +fill-tags|bcftools view -H|while read line; do echo $line|awk -F "\ " '{print$1"_"$2}'|tr "\n" "\t";echo $line|awk -F ";" '{print$2}'|sed -e "s/AC=//g"; done > DEEP_allele_counts.txt

# Arrange these as required in the input:

cut -f 1 TOP_allele_counts.txt|head|while read line
do 
echo "${line}_Nr_alleles"|tr "\n" "\t"
echo "2"
grep "$line" DEEP_allele_counts.txt|tr "\n" "\t"
variable=$(grep "$line" DEEP_allele_counts.txt|awk -F "\t" '{print$2}')
let "result = 36 - $variable"
echo $result
grep "$line" MID_allele_counts.txt|tr "\n" "\t"
variable=$(grep "$line" MID_allele_counts.txt|awk -F "\t" '{print$2}')
let "result = 36 - $variable"
echo $result
grep "$line" TOP_allele_counts.txt|tr "\n" "\t"
variable=$(grep "$line" TOP_allele_counts.txt|awk -F "\t" '{print$2}')
let "result = 36 - $variable"
echo $result
done

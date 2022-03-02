#!/bin/bash

#Procedure to make DNA Arlequin file from VCF for Approximate Bayesian Computation:

#Run in conda environment with BCftools, Bedtools installed, bbmap, PGSpider. 

#Before starting: split multisample vcf/bcf into single samples. Ex)
	#bcftools +split input.vcf.gz -Oz -o Out_Dir
	#recommend spliting a VCF that was filtered to sites just overlapping with random bed regions.
#Sample random chuncks of the genome: Ex)
#If using WGS this might be something like RAD tags:
	#bedtools random -g genome_file -l 150 -n 10000 > random.bed

#Generate Sequences with varinats applied based on their haplotypes:
#Be sure to generate use bcftools norm on the split vcf files and index.
#requires bcftools and bedtools

ls *.vcf.gz|awk -F "." '{print$1}' |while read line
do
cat SC_F1-1A_hap13B_scaffold_LOD26_complete.fa| bcftools consensus -H 1 $line.vcf.gz > temp1
bedtools getfasta -fi temp1 -bed random.bed > ${line}_hap1.fa
cat SC_F1-1A_hap13B_scaffold_LOD26_complete.fa| bcftools consensus -H 2 $line.vcf.gz > temp2 
bedtools getfasta -fi temp2 -bed random.bed > ${line}_hap2.fa
rm temp*.fai
done

#Generate separate files for each Chromosome for each sample:
#requires bbmap

ls *_hap1.fa|awk -F "_" '{print$1"_"$2"_"$3}'|while read line
do
for i in {01..12}
do
grep -A 1 ">CHR$i:" ${line}_hap1.fa > temp1.fa
fuse.sh in=temp1.fa out=${line}_CHR${i}_1.fa pad=0 name=${line}_1
grep -A 1 ">CHR$i:" ${line}_hap2.fa > temp2.fa
fuse.sh in=temp2.fa out=${line}_CHR${i}_2.fa pad=0 name=${line}_2
done
done

#Generate single multifasta file:

cat Tan_*_*_CHR*_*.fa > all_loci.fa


#Generate the Arp file with PGSpider

PGDspider2-cli -inputfile all_loci.fa -inputformat FASTA -outputfile Tanners_Observed_all_loci.arp -outputformat ARLEQUIN -spid fasta2arp.spid

#Arp Requires light editing for matching your simulations

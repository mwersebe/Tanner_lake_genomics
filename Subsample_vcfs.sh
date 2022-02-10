#!/bin/bash

#Script to sample lines from a VCF and convert to GENEPOP format for NeEstimator or PLINK for GONE estimation.

#spit the vcf into the variant lists and header portions.

#USEAGE: Run in a conda environment were bcftools and PGDSipder are installed. Requires sample program from: https://github.com/alexpreynolds/sample.git.

echo "Begining Process"

bcftools view -H -Ov VARIANTFILE > vcf_tail.txt 

bcftools view -h -Ov VARIANTFILE > vcf.header

echo "Finished Parsing Input VCF"

#Sample loop to sample lines from the vcft_tail file to construct a set of new sub-sampled vcf files

echo "Writing Output Files"

for i in {01..10} 
do
/home/giovanni/Programs/sample/sample -p -k 600000 vcf_tail.txt > vcf_sample$i.txt
cat vcf.header vcf_sample$i.txt > sample$i.vcf
plink --vcf sample$i.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --recode --out sample$i
done

#Sub in for GenePop conversion.
#PGDSpider2-cli -inputfile sample$i.vcf -inputformat VCF -outputfile sample$i.genepop.txt -outputformat GENEPOP -spid vcf_to_genepop.spid


echo "DONE working with VARIANTFILE" 






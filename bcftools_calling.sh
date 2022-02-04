
#Parallelize calling across 2Mb chuncks of the genome: ~20 minutes on an OSCER processor.

bcftools mpileup -f SC_F1-1A_hap13B_scaffold_LOD26_complete.fa -Ou -r CHR01:12000000-140000000 *.rg.bam |bcftools call -mv -Ob -o CHR01_12000000_14000000.bcf --thread 4


#Index each bcf file and concat each into a chromosome specific file:

bcftools index -c CHR10_6000000_7571162.bcf

bcftools concat -a -d all -Ob --threads 2 -o Tanners_CHR10.bcf -f CHR10.bcf.list

#Merge all bcf files to into a whole genome file:

bcftools index -c Tanners_CHR*.bcf

bcftools concat -a -d all -Ob --threads 2 -o Tanners_CHRall.bcf Tanners_CHR*.bcf

#Filter on quality and missingness to confident Bi-allelic SNPs:

bcftools filter -g3 -G10 -e'%QUAL<10 || (DP<10 && %QUAL<30 && MQ<10) || (AC==0)' -Ou Tanners_CHRall.bcf |bcftools view -m2 -M2 -v snps -Ou|bcftools +fill-tags |bcftools filter -e'F_MISSING > 0.1' -Ob -o Tanners_CHRall_filter_BiSnps.bcf

#Filter to SNPs in intron.gff3 regions:

bcftools view -Oz Tanners_CHRall_filter_BiSnps.bcf |bedtools intersect -a stdin -b intron.gff3 -header -wa |bgzip > Tanners_CHRall_filter_BiSnp_intron.vcf.gz

#LD Prune data set with PLINK

plink --vcf Tanners_CHRall_filter_BiSnps_intron.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# -indep-pairwise 50 10 0.1 --out Tanners_CHRall_filter_BiSnps_intron_LDprune

plink --vcf Tanners_CHRall_filter_BiSnps_intron.vcf.gz --allow-extra-chr --set-missing-var-ids @:# --extract Tanners_CHRall_filter_BiSnps_intron_LDprune.prune.in --out Tanners_CHRall_filter_BiSnps_intron_LDprune --make-bed --recode vcf bgz

#Move vcf file into R for PCA analysis:

#see script. 

#Filter VCF for alleles that change in frequency:

bcftools view -H Tanners_22-24_snps.vcf.gz -Ov|while read line; do echo $line|awk -F "\ " '{print$1"\t"$2"\t"$4"\t"$5}' |tr "\n" "\t"; echo $line|awk -F ";" '{print$3"\t"$5"\t"$6}'; done > Tan_22-24_alleleFreqs.tsv

#use sed to remove the AF=, MAF= and AN=

#Filter to polymorphic SNPs in 20-22 cm pop:

awk '{if ($6 >=0.142 && $6 <= 1.0) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' Tan_22-24_alleleFreqs.tsv > Tan_22-24_polymorphic.tsv

#make the first column Grep-able 

sed -e 's/^1/CHR01/g' Tan_22-24_polymorphic.tsv|sed -e 's/^2/CHR02/g' |sed -e 's/^3/CHR03/g'|sed -e 's/^4/CHR04/g'|sed -e 's/^5/CHR05/g'|sed -e 's/^6/CHR06/g'|sed -e 's/^7/CHR07/g'|sed -e 's/^8/CHR08/g'|sed -e 's/^9/CHR09/g'|sed -e 's/^10/CHR10/g'|sed -e 's/^11/CHR11/g'|sed -e 's/^12/CHR12/g' > temp


#Search for the same snps in the other files:

cat Tan_22-24_polymorphic.tsv |awk -F "\t" '{print$1" "$2}' |while read chr basepair; do grep "$chr" Tan_LC_alleleFreqs.tsv |grep -w "$basepair"|awk -F "\t" '{print$5"\t"$6"\t"$7}';done > Tan_LC_polymorphic.tsv &

paste -d "\t" Tan_22-24_polymorphic.tsv Tan_10-12_polymorphic.tsv Tan_LC_polymorphic.tsv > Tanners_polymorphic.tsv




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

# grab the allele freqs in each temporal pop: SUB in TOP/MID for each:

bcftools view -S DEEP_samples.txt Tanners_CHRall_filter_BiSnps_LDprune_nomiss.vcf.gz|bcftools +fill-tags|bcftools view -H|while read line; do echo $line|awk -F "\ " '{print$1"\t"$2"\t"$4"\t"$5}' |tr "\n" "\t"; echo $line|awk -F ";" '{print$3"\t"$5"\t"$6}'; done |sed -e "s/AN=//g"|sed -e "s/AF=//g"|sed -e "s/M//g" > DEEP_allele_freqs.tsv

cut -f 6,7 {POP}_allele_freqs.tsv > {POP}.temp

paste -d "\t" DEEP_allele_freqs.tsv {POP}.temp > All_allele_freqs.tsv



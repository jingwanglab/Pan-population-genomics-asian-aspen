#!/bin/bash
set -e
indir="/xxx/GATK_vcf"
outdir_indel="/xxx/filter/indel"
outdir_mark="/xxx/filter/remove_indel_mark"
outdir_snp="/xxx/filter/snp"
PIGZ="/xxx/software/pigz-2.4"
TABIX="/xxx/software/tabix-0.2.6"
BEDTOOLS="/data/apps/bedtools/bedtools-2.26.0/bin"
for vcf in $indir/*.vcf.gz
do
	temp=${vcf##*/}
	snp=${temp%v*}snp
	outfile_snp=${temp%v*}snp.recode.vcf
	indel=${temp%v*}indel
	outfile_indel=${temp%v*}indel.recode.vcf
	b=${temp%v*}indel.recode.vcf.gz
	outfile_mark=${temp%v*}rm_indel_mark.vcf
	echo "cd $outdir_snp;$VCFTOOLS/vcftools --gzvcf $vcf --remove-indels  --recode --recode-INFO-all --out $snp;$PIGZ/pigz $outfile_snp;cd $outdir_indel;$VCFTOOLS/vcftools --gzvcf $vcf --keep-only-indels  --recode --recode-INFO-all --out $indel;$PIGZ/pigz $outfile_indel;cd $outdir_mark;$BEDTOOLS/bedtools window -header -a $vcf -b $outdir_indel/$b -w 5 -u > $outfile_mark;$TABIX/bgzip $outfile_mark" >> filter.sh 
done

#!/bin/bash
set -e
TABIX="/xxx/software/tabix-0.2.6"
indir="/xxx/INDEL/filter1_gatk/filter"
outdir="/xxx/INDEL/filter2"
for vcf in $indir/*.gz
do
	temp=${vcf##*/}
	chr=${temp%indel.*}indel.bialleleic
	chr_gz=${temp%indel.*}indel.bialleleic.recode.vcf
	cd $outdir;vcftools --gzvcf $vcf --max-alleles 2 --min-alleles 2 --minDP 5 --minGQ 10 --max-missing 0.8 --maf 0.00001 --recode --recode-INFO-all --out $chr;$TABIX/bgzip $chr_gz ;done


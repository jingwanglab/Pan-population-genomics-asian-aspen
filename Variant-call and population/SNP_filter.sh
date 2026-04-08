#!/bin/bash
set -e
TABIX="/xxx/software/tabix-0.2.6"
VCFTOOLS="/data/apps/vcftools"
indir="/xxx/filter/snp"
outdir="/xxx/filter/snp/filter"

for vcf in $indir/*.gz
do
	temp=${vcf##*/}
	chr=${temp%recode.*}biallelic
	chr_gz=${temp%recode.*}biallelic.recode.vcf
	echo "cd $outdir;$VCFTOOLS/vcftools --gzvcf $vcf --max-alleles 2 --min-alleles 2 --minDP 5 --minGQ 10 --max-missing 0.8 --maf 0.00001 --recode --recode-INFO-all --out $chr;$TABIX/bgzip $chr_gz" >> filter.sh
done
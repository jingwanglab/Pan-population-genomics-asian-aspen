#!/bin/bash
set -e
indir="/xxx/INDEL/pd_indel"
outdir="/xxx/INDEL/filter1_gatk/filter"
sh_file="/xxx/INDEL/pd_indel_filter.sh"
GATK="/xxx/software/gatk-4.1.7.0/gatk"
BGZIP="/xxx/software/tabix-0.2.6/bgzip"
TABIX="/xxx/software/tabix-0.2.6/tabix"
reference="/xxx/reference/ZL2020-SXHZ.fasta"
for file in $indir/*.vcf.gz
do
	outfile=`basename $file .recode.vcf.gz`.QD2.FS200.SOR10.MQ-12.5.RPRS-8.vcf.gz	
	echo "cd $indir;$TABIX -p vcf $file;cd $outdir;$GATK VariantFiltration -R $reference -V $file -O $outfile --filter-expression \"QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\" --filter-name \"Filter\"" >> $sh_file
done


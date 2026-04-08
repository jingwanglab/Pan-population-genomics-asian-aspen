#!/bin/bash
#GATK-HaplotypeCaller
set -e
INDIR="/xx/picard"
OUTDIR="/xx/GATK_gvcf"
juqun_file="/w/00/u/user206/juqun.txt"
JAVA="/data/apps/jdk/jdk1.8.0_131/bin/java"
reference="/xx/reference/ZL2020-SXHZ.fasta"
GATK="/data/apps/gatk/gatk-4.0.5.1"

#create outdir
cd $OUTDIR
while read line
do
        if [ ! -d $line ];then
        mkdir $line
        fi
done < $juqun_file

#create g.vcf
for indir in $INDIR/*
do
	for file in `find  $indir -name "*sort.rmdup.bam"`
	do
		outdir=${indir##*/}
		infile=$file
		temp=${infile##*/}
		outfile_gvcf=${temp%sort*}g.vcf.gz
		echo "source ~/.bashrc;cd $OUTDIR/$outdir;/data/apps/jdk/jdk1.8.0_131/bin/java -Xmx16g -jar /data/apps/gatk/gatk-4.0.5.1/gatk-package-4.0.5.1-local.jar HaplotypeCaller --native-pair-hmm-threads 6 -R $reference -I $infile -O $outfile_gvcf -ERC GVCF --heterozygosity 0.01" >>  /w/00/u/user206/gatk.sh
	done
done

#GATK-Combinegvcf
INDIR="/xx/GATK_gvcf"
OUTDIR="/xx/combinegvcf"
beddir="/xx/reference/bed"

#create variant
inputfile=""

#create inputfile
for dir in $INDIR/*
do
        for file in $dir/*.g.vcf.gz
        do
                inputfile="$inputfile --variant $file"
        done
done
#choose -L
for bed in $beddir/*
do
	temp=${bed##*/}
	outfile=${temp%.*}.g.vcf.gz
	echo "source ~/.bashrc;cd $OUTDIR;/data/apps/jdk/jdk1.8.0_131/bin/java -Xmx36g -jar /data/apps/gatk/gatk-4.0.5.1/gatk-package-4.0.5.1-local.jar CombineGVCFs -R $reference $inputfile -L $bed -O $outfile" >> /w/00/u/user206/combinegvcf.sh
done

#GATK-GenotypeGVCFs
INDIR="/xx/combinegvcf"
OUTDIR="/xx/GATK_vcf"
bed_dir="/xx/reference/bed"
for bed in $bed_dir/*
do
	temp=${bed##*/}
	chr=${temp%.*}
	inputfile=${temp%.*}.222.g.recode.vcf.gz
	outputfile=${temp%.*}.222.vcf.gz
	echo "cd $OUTDIR;source ~/.bashrc;$GATK/gatk --java-options \"-Xmx36g\" GenotypeGVCFs -R $reference -L $bed -V $INDIR/$inputfile -O $outputfile --heterozygosity 0.01 -all-sites --indel-heterozygosity 0.0025 --allow-old-rms-mapping-quality-annotation-data" >> genotype.sh		
done

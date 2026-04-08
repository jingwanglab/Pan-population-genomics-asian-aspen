#!/bin/bash
samtools="/usr_storage/software/samtools-1.11/bin/samtools"
indir="/xxx/nextdenovo"
minimap2="/xxx/software/conda_evns/LRS_SV/bin/minimap2"
sniffles="/xxx/software/conda_evns/LRS_SV/bin/sniffles"
cuteSV="/xxx/software/conda_evns/cutesv/bin/cuteSV"
SURVIVOR="/xxx/software/SURVIVOR-master/Debug/SURVIVOR"
ref="/xxx/data/ZL2020-SXHZ.fasta"
for file in $indir/*
do
	temp=${file##*/}
	if [[ $temp =~ "fasta" ]];then
		species=${temp%.*}
		echo "$minimap2 -ax map-ont -R '@RG\tID:$species\tSM:$species' --MD -Y  -L $ref $file  | $samtools  sort -O BAM -T ./  -o $species.minimap2.sort.bam">$species.sh
		echo "$samtools index $species.minimap2.sort.bam">>$species.sh
		echo "  $cuteSV $species.minimap2.sort.bam $ref $species.minimap2.cutesv.vcf  ./  --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --threads 16 --min_size 50  --max_size 1000000 --min_support 10 "
		echo "$sniffles -m $species.minimap2.sort.bam -v $species.minimap2.sniffles.vcf -l 50 -n 10 ">>$species.sh
	    echo "$cuteSV --threads 1 --genotype -Ivcf 1_cutesv_merged_SURVIVOR_50dist_typesave_LG.vcf  $species.minimap2.sort.bam $ref  1_$species.minimap2.cutesv.gt.vcf ./  --min_size 50  --max_size 500000 --min_support 10 "
 echo "$SURVIVOR  merge sniffles_vcf_files_raw_calls.txt  1000 1 1 -1 -1 -1 sniffles_merged_SURVIVOR_1kbpdist_typesave.vcf "
	echo "nohup $sniffles -m $species.minimap2.sort.bam -v $species.minimap2.sniffles.gt.vcf --Ivcf sniffles_merged_SURVIVOR_1kbpdist_typesave.vcf  -l 50 -n 10 2>$species.sn.gt.log & "
	fi
done
ls | grep 'cutesv.gt.vcf' > cutesv_gt.txt
ls | grep 'cutesv.vcf'> cutesv_raw.txt
ls | grep 'sniffles.vcf' > sniffles_raw.txt
ls | grep 'sniffles.gt.vcf'>sniffles_gt.txt
$SURVIVOR  merge sniffles_raw.txt  1000 1 1 -1 -1 -1 sniffles_merged_raw.vcf
$SURVIVOR  merge sniffles_gt.txt 1000 -1 1 -1 -1 -1 sniffles_merged_gt.vcf
awk '{split($8,a,";");split(a[3],b,"=");if ($0 ~ "#") print $0 ;else if( b[2] <= 1000000 && b[2] >= -1000000)   print $0}' sniffles_merged_gt.vcf > sniflles_1M_merged_gt.vcf
$SURVIVOR  merge cutesv_raw.txt 1000 1 1 -1 -1 -1 cutesv_merged_raw.vcf
$SURVIVOR  merge cutesv_gt.txt 1000 -1 1 -1 -1 -1 cutesv_merged_gt.vcf
#merge sv 
$SURVIVOR  merge sniffles.cutesv.gt.vcf 1000 2 1 1 -1 50  cutesv_sniffles_sv.vcf

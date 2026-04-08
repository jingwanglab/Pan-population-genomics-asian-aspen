#!/bin/bash
# conda activate /xxx/conda/envs/py35/
vcf="/xxx/snp_indel_sv.maf0.05.mis0.8.vcf"
plink --vcf $vcf --allow-extra-chr --allow-no-sex  --make-bed  --out gemma_input
paste -d " " <(cut -d " " -f 1,2,3,4,5 gemma_input.fam) <(sed '1d' clim.points |cut -d " " -f 4- ) > gemma_input1.fam
paste -d " " <(cut -d " " -f 1,2,3,4,5 gemma_input.fam) <(less clim.points|grep 
gemma -bfile gemma_input -gk 2 -o gemma_matrix
for i in {1..10}
do	
	num=$(($i+3))
	bio=`head -n 1 clim.points|cut -d " " -f $num`
	gemma -bfile gemma_input -n $i -k output/gemma_matrix.sXX.txt -lmm -o  $bio.gemma
done

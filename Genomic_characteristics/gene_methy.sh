for individual in {ZL2020-SXHZ,ZSP191aL,ZSP192L,ZSP194L,ZSP196L};
do
cd /XX/gene_methy/${individual}/

bedtools flank -i ${individual}_gene.bed -g ${individual}_genome.len  -l 0 -r 2000 -s > ${individual}_gene_down2k.bed
bedtools flank -i ${individual}_gene.bed -g ${individual}_genome.len  -l 2000 -r 0 -s > ${individual}_gene_up2k.bed

for m in{CHG,CpG,CHH}
do
        bedtools intersect -a ${individual}_gene_down2k.bed -b /XX/methy_call/${individual}/${m}.bed -wa -wb > ${individual}_gene_down2k_${m}_intersect.bed
        bedtools intersect -a ${individual}_gene_up2k.bed -b /XX/methy_call/${individual}/${m}.bed -wa -wb > ${individual}_gene_up2k_${m}_intersect.bed
        bedtools intersect -a ${individual}_gene.bed -b /XX/methy_call/${individual}/${m}.bed -wa -wb > ${individual}_gene_${m}_intersect.bed
done
done

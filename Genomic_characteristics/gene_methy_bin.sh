for individual in {ZL2020-SXHZ,ZSP191aL,ZSP192L,ZSP194L,ZSP196L};
do
cd /XX/gene_methy/${individual}/

bedtools flank -i ${individual}_gene.bed -g ${individual}_genome.len  -l 0 -r 2000 -s > ${individual}_gene_down2k.bed
bedtools flank -i ${individual}_gene.bed -g ${individual}_genome.len  -l 2000 -r 0 -s > ${individual}_gene_up2k.bed
awk '($6=="-")' ${individual}_gene_down2k.bed | bedtools makewindows -b - -n 20 -i srcwinnum -reverse > ${individual}_gene_down2k_m.bed
awk '($6=="+")' ${individual}_gene_down2k.bed | bedtools makewindows -b - -n 20 -i srcwinnum > ${individual}_gene_down2k_p.bed
awk '($6=="-")' ${individual}_gene_up2k.bed | bedtools makewindows -b - -n 20 -i srcwinnum -reverse > ${individual}_gene_up2k_m.bed
awk '($6=="+")' ${individual}_gene_up2k.bed | bedtools makewindows -b - -n 20 -i srcwinnum > ${individual}_gene_up2k_p.bed
awk '($6=="-")' ${individual}_gene.bed | bedtools makewindows -b - -n 30 -i srcwinnum -reverse > ${individual}_gene_m.bed
awk '($6=="+")' ${individual}_gene.bed | bedtools makewindows -b - -n 30 -i srcwinnum > ${individual}_gene_p.bed
cat ${individual}_gene_down2k_m.bed ${individual}_gene_down2k_p.bed > ${individual}_gene_down2k_20bin.bed
cat ${individual}_gene_up2k_m.bed ${individual}_gene_up2k_p.bed > ${individual}_gene_up2k_20bin.bed
cat ${individual}_gene_m.bed ${individual}_gene_p.bed > ${individual}_gene_30bin.bed

for m in {CHG,CpG,CHH}
do
        bedtools intersect -a ${individual}_gene_down2k_20bin.bed -b /XX/methy_call/${individual}/${m}.bed -wa -wb > ${individual}_gene_down2k_20bin_${m}_intersect.bed
        bedtools intersect -a ${individual}_gene_up2k_20bin.bed -b /XX/methy_call/${individual}/${m}.bed -wa -wb > ${individual}_gene_up2k_20bin_${m}_intersect.bed
        bedtools intersect -a ${individual}_gene_30bin.bed -b /XX/methy_call/${individual}/${m}.bed -wa -wb > ${individual}_gene_30bin_${m}_intersect.bed
done
done
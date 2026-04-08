#Pairwise comparison of five genomes by mummer
for i in {ZL2020-SXHZ,ZSP191aL,ZSP192L,ZSP194L,ZSP196L};for j in {ZL2020-SXHZ,ZSP191aL,ZSP192L,ZSP194L,ZSP196L}; do ../bin/nucmer --maxmatch -l 40 -g 90 -b 100 -c 200 -p ./$i_$j ../Chr/$i.Chr ../Chr/$j.Chr ;done ;done

#Filter alignments with length less than 100 bp and similarity less than 90%
for i in {ZSP191aL,ZSP192L,ZSP194L,ZSP196L};do ../bin/delta-filter -m -i 90 -l 100 ../mummer/nucmer/${i}_ZL2020-SXHZ.delta > ../mummer/delta-filter/${i}_ZL2020-SXHZ._m_i90_l100.delta & done

#Convert the alignment results into coords format for subsequent analysis
do ../bin/show-coords  -THrd ../mummer/delta-filter/${i}_ZL2020-SXHZ._m_i90_l100.delta > ${i}_ZL2020-SXHZ_m_i90_l100.coords & done

#*.wga.block.txt prepare the bed formated file "xx.wga.block.txt" from the .coords file, e.g: ZL2020-SXHZ.wga.block.txt
for i in {ZSP191aL,ZSP192L,ZSP194L,ZSP196L};do cd /xxx/pangenome/ZL2020-SXHZ/${i};awk '{print $10"\t"$1"\t"$2"\t"$11"\t"$3"\t"$4}' ${i}_ZL2020-SXHZ_m_i90_l100.coords > ${i}.wga.block.txt ;done

# prepare the bed formated file "xx.del.bed" and "xx.ins.bed".
for i in {ZSP191aL,ZSP192L,ZSP194L,ZSP196L};do cd /xxx/pangenome/ZL2020-SXHZ/${i};cut -f 1-3 ${i}.wga.block.txt |sort -k1,1 -k2,2n -k3,3n > ZL2020-SXHZ.${i}.ref.aln.bed;done
for i in {ZSP191aL,ZSP192L,ZSP194L,ZSP196L};do cd /xxx/pangenome/ZL2020-SXHZ/${i};cut -f 4-6 ${i}.wga.block.txt |awk '{if ($2>$3) print $1"\t"$3"\t"$2; else print }' |sort -k1,1 -k2,2n -k3,3n >ZL2020-SXHZ.${i}.qry.aln.bed ; done
for i in {ZSP191aL,ZSP192L,ZSP194L,ZSP196L};do cd /xxx/pangenome/ZL2020-SXHZ/${i};/usr_storage/software/bedtools2/bin/complementBed -i ZL2020-SXHZ.${i}.ref.aln.bed -g /xxx/pangenome/chrlen/${i}.len > ${i}.del.bed ; done
for i in {ZSP191aL,ZSP192L,ZSP194L,ZSP196L};do cd /xxx/pangenome/ZL2020-SXHZ/${i};/usr_storage/software/bedtools2/bin/complementBed -i ZL2020-SXHZ.${i}.qry.aln.bed -g /xxx/pangenome/chrlen/${i}.len > ${i}.ins.bed ; done

#Pan-genome computing based on genome sequences
for i in {ZSP191aL,ZSP192L,ZSP194L,ZSP196L};do python -u wga.pangenome.py -w /xxx/ZL2020-SXHZ/${i}; -o ./ -g ./chrlength/ ;done

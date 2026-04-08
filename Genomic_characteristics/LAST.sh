#Building an index
for i in {ZSP191aL,ZSP192L,ZSP194L,ZSP196L};do lastdb -P0 -uNEAR -R01  $i-NEAR ../Chr/${i}.Chr & done
#Determine replacement and gap frequencies
for i in {ZSP191aL,ZSP192L,ZSP194L,ZSP196L};do last-train -P0 --revsym --matsym --gapsym -E0.05 -C2 ../LAST/lastdb/${i}-NEAR ../Chr/${i}.Chr > $i-ZL2020-SXHZ.mat & done
#Many-to-one comparison
for i in {ZSP191aL,ZSP192L,ZSP194L,ZSP196L};do lastal -m50 -E0.05 -C2 -p ../LAST/last-train/$i-ZL2020-SXHZ.mat ../LAST/lastdb/${i}-NEAR ../Chr/ZL2020-SXHZ.Chr | last-split -m1 > ${i}-ZL2020-SXHZ-1.maf;done
#Aligned to reference genome
for i in {ZSP191aL,ZSP192L,ZSP194L,ZSP196L};do maf-swap ../LAST/lastal/${i}-ZL2020-SXHZ-1.maf |awk '/^s/ {$2 = (++s % 2 ? "ZL2020-SXHZ." : "${i}.") $2} 1' |last-split -m1 |maf-swap > ${i}-ZL2020-SXHZ-2.maf & done
#Discard single sequence alignments, convert the results to tab format, and discard those with an error rate > 10^-5
for i in {ZSP191aL,ZSP192L,ZSP194L,ZSP196L};do last-postmask ../LAST/maf-swap/${i}-ZL2020-SXHZ-2.maf |maf-convert -n tab |awk -F'=' '$2 <= 1e-5' > ${i}-ZL2020-SXHZ.tab & done





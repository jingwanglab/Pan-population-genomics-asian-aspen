#Pairwise comparison of five genomes
for i in {ZL2020-SXHZ,ZSP191aL,ZSP192L,ZSP194L,ZSP196L};for j in {ZL2020-SXHZ,ZSP191aL,ZSP192L,ZSP194L,ZSP196L}; do ../bin/nucmer --maxmatch -l 40 -g 90 -b 100 -c 200 -p ./$i_$j ../Chr/$i.Chr ../Chr/$j.Chr ;done ;done

#Filter alignments with length less than 100 bp and similarity less than 90%
for i in {ZSP191aL,ZSP192L,ZSP194L,ZSP196L};do ../bin/delta-filter -m -i 90 -l 100 ../mummer/nucmer/${i}_ZL2020-SXHZ.delta > ../mummer/delta-filter/${i}_ZL2020-SXHZ._m_i90_l100.delta & done

#Convert the alignment results into coords format for subsequent analysis
do ../bin/show-coords  -THrd ../mummer/delta-filter/${i}_ZL2020-SXHZ._m_i90_l100.delta > ${i}_ZL2020-SXHZ_m_i90_l100.coords & done

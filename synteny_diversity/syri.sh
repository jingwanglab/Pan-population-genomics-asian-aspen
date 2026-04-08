#Use the coords file of the mummer alignment results to perform SYRI pairwise alignment
for i in {ZL2020-SXHZ,ZSP191aL,ZSP192L,ZSP194L,ZSP196L}
do for j in {ZL2020-SXHZ,ZSP191aL,ZSP192L,ZSP194L,ZSP196L}
do cd ../syri/$j/${i};nohup /usr_storage/software/syri/syri/bin/syri -c ../mummer/show-coords/${i}_$j_m_i90_l100.coords -d ../mummer/delta-filter/${i}_${j}._m_i90_l100.delta -r ../data-re/Chr/$i.Chr -q ../data-re/Chr/${j}.Chr --nc 5 --all -k & done ;done
#To visualize SYRI results by plotsr
for i in {ZSP191aL,ZSP192L,ZSP194L,ZSP196L};do cd ../syri/ZL2020-SXHZ/${i};/usr_storage/software/syri/syri/bin/plotsr syri.out ../data-re/Chr/ZL2020-SXHZ.Chr ../data-re/Chr/${i}.Chr  -H 8 -W 5 -o pdf &done


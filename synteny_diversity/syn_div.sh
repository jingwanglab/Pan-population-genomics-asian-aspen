#Prepare the block file for the collinear area of genome pairwise alignment
for i in {ZL2020-SXHZ,ZSP191aL,ZSP192L,ZSP194L,ZSP196L};do for j in {ZL2020-SXHZ,ZSP191aL,ZSP192L,ZSP194L,ZSP196L};do cd ../syri/$j/${i};grep SYNAL syri.out|cut -f 1-3,6-8,10 > ../SD_HR/$j/${i}/${i}.wga.syn.block.txt ; done ;done

## Calculate all collinear regions in each genome separately
for i in {ZL2020-SXHZ,ZSP191aL,ZSP192L,ZSP194L,ZSP196L};do cd ../SD_HR/$i;/usr_storage/software/bedtools2/bin/multiIntersectBed -i ./*/*.wga.syn.block.txt -names ZSP191aL ZSP192L ZSP194L ZSP196L >ZL2020-SXHZ.syn.txt & done

## Calculate the alignment results on each chromosome of each genome
for k in {ZSP191aL,ZSP192L,ZSP194L,ZSP196L}; do for j in {1..19}; do /usr_storage/software/mummer-4.0.0beta2/bin/show-aligns -r ../SD_HR/ZL2020-SXHZ/$k/out_m_i90_l100.delta Chr$j Chr$j > ../SD_HR/ZL2020-SXHZ/$k/out_m_i90_l100.Chr$j.aligns & done &done

#Prepare files in croods format by get.all.syn.coord.pl
for i in {ZL2020-SXHZ,ZSP191aL,ZSP192L,ZSP194L,ZSP196L};do cd ../SD_HR/${i}; perl get.all.syn.coord.pl ./${i}.syn.all.txt ../SD_HR/${i}/ ./${i}.syn.all.coords.txt & done

#Use run.cal.syn.div.pl to get the colinear diversity of each chromosome
for i in {ZL2020-SXHZ,ZSP191aL,ZSP192L,ZSP194L,ZSP196L};do cd ../SD_HR/ZSP191aL/${i}; perl ./run.cal.syn.div.pl ./${i}.syn.all.coords.txt ../SD_HR/${i}/ ../SD_HR/chrlen/ ./calculate.syn.diversity.pairwise.chr.pl ./splitChr/ > np.chr;done

for i in {ZL2020-SXHZ,ZSP191aL,ZSP192L,ZSP194L,ZSP196L};do for j in {1..19};do cd ../SD_HR/$i;nohup perl ./syn.div.merge.pl ./splitChr/ Chr${j} ../SD_HR/Chrbed_5/ZL2020-SXHZ.leng.txt ./Chr${j}.syn.div.pos.txt;done;done

#Calculate collinear diversity using calculate.syn.diversity.window.pl with window and step 
for k in {1..19}; do perl calculate.syn.diversity.window.pl ./Chr$k.syn.div.pos.txt 10000 2000 ./Chr$k.syn.div.win10kb.step2kb.txt & done &
for k in {1..19}; do perl calculate.syn.diversity.window.pl ./Chr$k.syn.div.pos.txt 200000 200000 ./Chr$k.syn.div.win200kb.step200kb.txt & done &

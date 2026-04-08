bowtie2_path='/usr_storage/software/bowtie/bowtie2-2.3.2/bin/'
samtools_path='/usr_storage/software/samtools-1.11/bin/samtools'
for individual in {ZL2020-SXHZ,ZSP191aL,ZSP192L,ZSP194L,ZSP196L};
do
java -jar /usr_storage/software/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 8 /XX/methy_reads/${individual}_1.fq.gz /XX/methy_reads/${individual}_2.fq.gz /XX/methy_reads/${individual}_1.clean.fq.gz /XX/methy_reads/${individual}_1.un.fq.gz /XX/methy_reads/${individual}_2.clean.fq.gz /XX/methy_reads/${individual}_2.un.fq.gz ILLUMINACLIP:/data/apps/trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:36
bismark_genome_preparation --bowtie2 --path_to_aligner $bowtie2_path --verbose /XX/Populus_pangenome/methy_call/${individual}/
bismark --bowtie2 --path_to_bowtie2 $bowtie2_path -p 10 --genome_folder /XX/Populus_pangenome/methy_call/${individual}/ -N 1 -L 20 --quiet --un --ambiguous --sam -o /XX/Populus_pangenome/methy_call/${individual}/ -1 /XX/methy_reads/${individual}_1.clean.fq.gz -2 /XX/methy_reads/${individual}_2.clean.fq.gz
deduplicate_bismark -p --samtools_path $samtools_path /XX/Populus_pangenome/methy_call/${individual}/${individual}.clean_bismark_bt2_pe.sam
bismark_methylation_extractor -p --gzip --bedGraph --comprehensive --no_overlap --counts --buffer_size 20G --report --cytosine_report --CX --samtools_path $samtools_path --genome_folder /XX/Populus_pangenome/methy_call/${individual}/ /XX/Populus_pangenome/methy_call/${individual}/${individual}.clean_bismark_bt2_pe.deduplicated.bam
cd /XX/Populus_pangenome/methy_call/${individual}/
bismark2bedGraph --cx CpG_context_${individual}.clean_bismark_bt2_pe.deduplicated.txt.gz -o CpG.bedgraph
bismark2bedGraph --cx CHG_context_${individual}.clean_bismark_bt2_pe.deduplicated.txt.gz -o CHG.bedgraph
bismark2bedGraph --cx CHH_context_${individual}.clean_bismark_bt2_pe.deduplicated.txt.gz -o CHH.bedgraph
gzip -d CpG.bedgraph.gz.bismark.cov.gz;awk '{a=$5+$6; if (a>3) print$1"\t",$2=$2-1"\t",$3"\t",$4"\t",$5"\t",$6=$5+$6}' CpG.bedgraph.gz.bismark.cov|sort -k1,1 -k2n,2 > CpG.bed;sed -i 's/ //g' CpG.bed
gzip -d CHG.bedgraph.gz.bismark.cov.gz;awk '{a=$5+$6; if (a>3) print$1"\t",$2=$2-1"\t",$3"\t",$4"\t",$5"\t",$6=$5+$6}' CHG.bedgraph.gz.bismark.cov|sort -k1,1 -k2n,2 > CHG.bed;sed -i 's/ //g' CHG.bed
gzip -d CHH.bedgraph.gz.bismark.cov.gz;awk '{a=$5+$6; if (a>3) print$1"\t",$2=$2-1"\t",$3"\t",$4"\t",$5"\t",$6=$5+$6}' CHH.bedgraph.gz.bismark.cov|sort -k1,1 -k2n,2 > CHH.bed;sed -i 's/ //g' CHH.bed
done

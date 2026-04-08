# Pan-population-genomics-Asian-aspen

Scripts for Chen et al. Genomic signatures of climate adaptation and maladaptation revealed through the integration of pangenomics and population genomics in a keysone forest tree species of East Asia. submitted



1.1 Genome assembly and alignment
Genome assembly:
Four Nanopore sequenced genomes were assembled to the chromosome level using the reference genome.
Ragtag.sh - Script to assemble four Nanopore sequenced individuals to the chromosome level.

Genome alignment:
LAST.sh - Script to align the four chromosome level genome with the reference genome.

1.2 Genome annotation
TE annotation:
First, we performed command ‘perl EDTA.pl --genome genome.fasta --sensitive 1 -anno 1’ implemented in EDTA software to annotate TE. Then, we employed TEsorter to reclassfy reclassify TEs that were annotated as “LTR/unknown” by EDTA, using the command 'TEsorter LTR_unknown.fasta -db rexdb-plant -p 6'.

Gene annotation:
Gene annotation was performed in the pipeline https://github.com/jingwanglab/Populus_genomic_prediction_climate_vulnerability/tree/main/1-Genome_analyses

1.3 TE analysis
For each individual, run the following steps:
(1) Extract non-homolog regions in MAF file.
$python 1.maf_extract_nonhomo.py -i ./maf_file/species_as_ref.maf -o ./non_homolog/species_nonhomolog
maf_extract_nonhomo.py - Script to extract non-homolog regions in MAF file.

(2) Merge non-homolog regions.
$sort -k1,1 -k2,2n species_nonhomolog | awk -F '\t' '{print $1"\t"$2"\t"$3"\t"($2+$3)"\t"$4}' | awk -F '\t' '{print $1"\t"$2"\t"$4"\t"$5}' | bedtools merge -i - > ./non_homolog/species_nonhomolog_sort.merge.bed

(3) Obtain homolog regions.
$bedtools complement -i ./non_homo/species_nonhomolog_sort.merge.bed -g ./genome_bed/species_genome.bed > ./homo/species_homolog.bed

(4) Obtain TEs present in individual-specific region.
$bedtools intersect -a ./non_homo/species_nonhomolog_sort.merge.bed -b species_TE.bed -wao > species_nonhomolog_TE.intersect

(5) Obtain TEs present in shared region.
$bedtools intersect -a ./homo/species_homolog.bed -b species_TE.bed -wao > species_homolog_TE.intersect

TE cluster:
cluster_TE_build.py - Script to generate TE files of pairwise genomes for cluster.

TE_vmatch_cluster.sh - Script to cluster syntentic fl-LTRs for pairwise genomes.

TE_analysis.sh - Script to integrate TE analysis.

1.4.Expression analysis
(1) gene_expression.sh - Script to obtain gene expression (FPKM) based on RNA-seq data.

(2) gene_TAU.py - Script to calculate tissue-specific expression for each of gene.

1.5.Methylation analysis
(1) methy_call.sh - Script for methylation calling.

(2) gene_methy.sh - Script to obtain methylation level in genes (body and flanking 2 kb regions).

(3) gene_methy_bin.sh - Script to obtain methylation level in genes (body and flanking 2 kb regions based on divided bins).

2.Pan-genome analysis
(1) Genome_pangenome.sh - Script to construct sequence-based pan-genome and align paired genome sequences using MUMMER. 

(2) wga.pangenome.py - Script to extract core, variable and specific region in each genome.

(3) Gene_pangenome.sh - Script to construct gene-based pan-genome.

(4) pan_gene_compare.py - Script to compare genome characteristics of different types of pan genes.

3.Synteny diversity analysis
(1) Use the coords file of the mummer alignment results to perform SYRI pairwise alignment.
mummer.sh - Script to obtain mummer alignment results.

(2) visualize SYRI results by plotsr.
syri.sh - Script to visualize the results of syri.

## Before caculate synteny diversity, please run all pairwise whole genome comparison using MUMmer and run SyRi to identify the syntenic and rearranged regions for each comparison.
Let's assume all the alignments in a folder like below:
	/xxx/SD_HR
	/xxx/SD_HR/ZL2020-SXHZ
	/xxx/SD_HR/ZL2020-SXHZ/ZSP191aL
	/xxx/SD_HR/ZL2020-SXHZ/ZSP192L
	...
	/xxx/SD_HR/ZL2020-SXHZ/ZSP196L
	...
	/xxx/SD_HR/ZSP196L/ZL2020-SXHZ
	/xxx/SD_HR/ZSP196L/ZSP191aL
	...
	
(3) Prepare the block file for the collinear area of genome pairwise alignment.
for i in {ZL2020-SXHZ,ZSP191aL,ZSP192L,ZSP194L,ZSP196L};do for j in {ZL2020-SXHZ,ZSP191aL,ZSP192L,ZSP194L,ZSP196L};do cd ../syri/$j/${i};grep SYNAL syri.out|cut -f 1-3,6-8,10 > ../SD_HR/$j/${i}/${i}.wga.syn.block.txt ; done ;done

(4) Calculate all collinear regions in each genome separately.
for i in {ZL2020-SXHZ,ZSP191aL,ZSP192L,ZSP194L,ZSP196L};do cd ../SD_HR/$i;/usr_storage/software/bedtools2/bin/multiIntersectBed -i ./*/*.wga.syn.block.txt -names ZSP191aL ZSP192L ZSP194L ZSP196L >ZL2020-SXHZ.syn.txt & done

(5) Calculate the alignment results on each chromosome of each genome.
for k in {ZSP191aL,ZSP192L,ZSP194L,ZSP196L}; do for j in {1..19}; do /usr_storage/software/mummer-4.0.0beta2/bin/show-aligns -r ../SD_HR/ZL2020-SXHZ/$k/out_m_i90_l100.delta Chr$j Chr$j > ../SD_HR/ZL2020-SXHZ/$k/out_m_i90_l100.Chr$j.aligns & done &done

(6) Prepare files in croods format by get.all.syn.coord.pl.
for i in {ZL2020-SXHZ,ZSP191aL,ZSP192L,ZSP194L,ZSP196L};do cd ../SD_HR/${i}; perl get.all.syn.coord.pl ./${i}.syn.all.txt ../SD_HR/${i}/ ./${i}.syn.all.coords.txt & done
get.all.syn.coord.pl - Script to get coordinates of syntenic regions in all genomes.

(7) Use run.cal.syn.div.pl to get the colinear diversity of each chromosome
for i in {ZL2020-SXHZ,ZSP191aL,ZSP192L,ZSP194L,ZSP196L};do cd ../SD_HR/ZSP191aL/${i}; perl ./run.cal.syn.div.pl ./${i}.syn.all.coords.txt ../SD_HR/${i}/ ../SD_HR/chrlen/ ./calculate.syn.diversity.pairwise.chr.pl ./splitChr/ > np.chr;done
for i in {ZL2020-SXHZ,ZSP191aL,ZSP192L,ZSP194L,ZSP196L};do for j in {1..19};do cd ../SD_HR/$i; perl ./syn.div.merge.pl ./splitChr/ Chr${j} ../SD_HR/Chrbed_5/ZL2020-SXHZ.leng.txt ./Chr${j}.syn.div.pos.txt ;done;done
run.cal.syn.div.pl - Script to caculate synteny diversity for every postion of each chromosome.
syn.div.merge.pl - Script to merge and calculate synteny diversity for every postion of genome.

(7) Calculate collinear diversity using calculate.syn.diversity.window.pl with window and step. 
for k in {1..19}; do perl calculate.syn.diversity.window.pl ./Chr$k.syn.div.pos.txt 10000 2000 ./Chr$k.syn.div.win10kb.step2kb.txt & done &
for k in {1..19}; do perl calculate.syn.diversity.window.pl ./Chr$k.syn.div.pos.txt 200000 200000 ./Chr$k.syn.div.win200kb.step200kb.txt & done &
calculate.syn.diversity.window.pl - Script to calculate synteny diversity in a given window	.

4. Variant-calling and population genetics
(1) trimmomatic.sh - Script to use Trimmomatic to filter raw reads.

(2) bwa&picard.sh - Script to use BWA to map reads to reference genome, use SAMtools to sort the alignment results and use Picard to mark PCR duplicate.

(3) GATK-HaplotypeCaller.sh - Script to perform SNP and Indel calling by GATK.

(4) SNP_INDEL_seperate.sh - Scripts to SNPs and indels.

(5) SNP_filter.sh - Scripts to filter SNPs.

(6) INDEL_GATK.sh - Scripts to filter indels by GATK. 

(7) INDEL_filter.sh - Scripts to filter indels by vcftools. 

(8) SV_cutesv_sniffle.sh - Scripts to calculates Svs by using cutesv and sniffle, and finally merges them using SURVIVOR.
All filtered SNPs, indels, and SVs datasets were merged using GATK for subsequent analysis.

(9) structure.sh - Script to caculate population admixture.

(10) PCA.sh - Script to caculate population PCA.

(11) PCA.R - Script to visualize the population PCA.


5.GEA (Genotype-Environment Association) analyses are conducted using LFMM, GEMMA and RDA.
(1) environment_cor.R - Script to calculate the correlation between environmental variables.

(2) GF.R - Script to calculate the importance of the environment using GF.

(3) lfmm.sh - Script to calculate the GEA by LFMM.

(4) GEMMA.sh - Script to calculate the GEA by GEMMA.

(5) RDA.R - Script to calculate the GEA by RDA.

(6) IBD_IBE.sh - Script to perform IBD, IBE, pIBD and pIBE analyses.


6.Genetic offset

(1)Local_Reverse_genetic_offset.R was used to calculate local genetic offset, which assumes in situ tolerance and reverse genetic offset to assess the maladaptation of populations when simultaneously considering the contributions of in situ adaptation and migration to future shifting climates.In addition, the unique contributions of SVs were calculated and projected onto the map in four equal parts from 0 to 1.

(2) Forwrad_genetic_offset.R was used for the forwad genetic offset by defining different geographic distances.

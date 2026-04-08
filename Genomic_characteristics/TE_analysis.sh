#For each species, run the following steps:
#Extract non-homolog regions in MAF file.
$python maf_extract_nonhomo.py -i ./maf_file/species_as_ref.maf -o ./non_homolog/species_nonhomolog

# Merge non-homolog regions.
$sort -k1,1 -k2,2n species_nonhomolog | awk -F '\t' '{print $1"\t"$2"\t"$3"\t"($2+$3)"\t"$4}' | awk -F '\t' '{print $1"\t"$2"\t"$4"\t"$5}' | bedtools merge -i - > ./non_homolog/species_nonhomolog_sort.merge.bed

#Obtain homolog regions.
$bedtools complement -i ./non_homo/species_nonhomolog_sort.merge.bed -g ./genome_bed/species_genome.bed > ./homo/species_homolog.bed

#Obtain TEs present in species-specific region.
$bedtools intersect -a ./non_homo/species_nonhomolog_sort.merge.bed -b species_TE.bed -wao > species_nonhomolog_TE.intersect

#Obtain TEs present in shared region.
$bedtools intersect -a ./homo/species_homolog.bed -b species_TE.bed -wao > species_homolog_TE.intersect


#Generate TE files of pairwise  genomes for cluster by cluster_TE_build.py .
#Cluster syntentic fl-LTRs for pairwise poplar genomes.
for dir in ./*
do
	cd $dir
	if [ -f "ltr_junction.fa" ];then
	mkdir 70
	mkdir 80
	mkdir 90
	mkvtree -db ltr_junction.fa -dna -pl -allout -v && vmatch -dbcluster 70 70 70/Cluster -s -identity 70 -exdrop 6 -seedlength 10 -d ltr_junction.fa && vmatch -dbcluster 80 80 80/Cluster -s -identity 80 -exdrop 5 -seedlength 15 -d ltr_junction.fa &&	vmatch -dbcluster 90 90 90/Cluster -s -identity 90 -exdrop 4 -seedlength 20 -d ltr_junction.fa
	cd ../
	else
	echo "$dir no file"
	cd ../
	fi
done
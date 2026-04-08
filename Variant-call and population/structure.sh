#### LD-pruned
vcftools --gzvcf pd_chr.vcf.gz --plink --out pd
plink --file pd --maf 0.05  --indep-pairwise 50 10 0.2 --allow-extra-chr --allow-no-sex --make-bed --out pd
sed -i -e '/Contig/d' pd.prune.in
plink --file pd --extract pd.prune.in --recode --allow-extra-chr --allow-no-sex --make-bed --out pd_ldpruned

#### use Admixture to perform structure analyses 
#!/bin/bash
   set -e
 file="/xxx/final_vcf/vcf/structure/pd_ldpruned.bed"
   outdir="/xxx/final_vcf/vcf/structure/K"
   AD="/data/apps/admixture/1.3.0/"
   for K in 1 2 3 4 5 6 7 8 9 10
   do
          echo "cd $outdir;$AD/admixture --cv $file $K|tee log${K}.out" >> /xxx/pd_structure.sh
   done

#Gene-based pan-genome computation of five genome using OrthFinder
/xxx/orthofinder -f /xxx/orthofinder/pep -S blast -t 10

#Extract pan-gene family and pan genes in Orthifinder output.
python /xxx/orthofinder/pan_gene_extract.py Orthogroups.tsv 

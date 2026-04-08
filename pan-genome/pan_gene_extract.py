

# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 11:06:55 2021

@author:Stt
"""

import pandas as pd


data= open('./Orthogroups.tsv').readlines()
data_m = data[1:]

core_lines = []
dis2_lines=[]
dis3_lines=[]
dis4_lines=[]
private_lines=[]


##### obtain pan-gene family
for line in data_m:
    xx = line.strip().split('\t')
    spe_i = 0
    for i in range(1,5):
#        ids = xx[i].split(',')
        if len(xx[i]) >0 and xx[i]!='':
            spe_i += 1
                
    if spe_i ==5:
        core_lines.append(line)
        continue
    elif spe_i ==1:
        private_lines.append(line)
        continue
    elif spe_i ==2:
        dis2_lines.append(line)
        continue
    elif spe_i ==3:
        dis3_lines.append(line)
        continue    
    elif spe_i ==4:
        dis4_lines.append(line)
        continue    


dis_lines =dis1_lines+ dis2_lines + dis3_lines + dis4_lines + dis5_lines + 

types = ['core','dis','private']

for ty in types:
    f_o = open('/xxx/pan_gene/'+ty+'/'+ty+'_og.tsv','w')
    for line in softcore_lines:
        f_o.write(line)
    f_o.close()
    
    
    ## extract gene ID of each species in pan-gene family
    ogs = pd.read_csv('/xxx/pan_gene/'+ty+'/'+ty+'_og.tsv', sep='\t')    
    spe = ogs.columns.tolist()
    spe = spe[1:]
      
    a_genes=[]
    for sp in spe:
        sp_df = ogs[ogs[sp].notna()]
        sp_genes = sp_df[sp].tolist()
        all_genes = []
        for ll in sp_genes:
            xx = ll.split(',')
            for i in xx:
                all_genes.append(i)
                a_genes.append(i)
        f = open('/xxx/pan_gene/'+ty+'/'+sp+'_gene.list','w')
        for gene in all_genes:
            f.write(gene+'\n')
        f.close()
         
    f = open('/xxx/pan_gene/'+ty+'/all_gene.list','w')
    for gene in a_genes:
        f.write(gene+'\n')
    f.close()
    
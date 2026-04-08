# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 11:57:55 2020

@author: Stt
"""

import itertools
import os

path = '/XX/TE_analysis/fl_LTR/'
files = os.listdir(path) 

combines = list(itertools.combinations(files,2))

for combine in combines:
    s1 = combine[0].split('_')[0]
    s2 = combine[1].split('_')[0]
    s = s1+'_'+s2
    mkdir_cmd = 'mkdir '+s
    cat_cmd = 'cat /XX/TE_analysis/fl_LTR/'+combine[0]+' /XX/TE_analysis/fl_LTR/'+combine[1]+' > '+s+'/ltr_juction.fa'
    os.system(mkdir_cmd)
    os.system(cat_cmd)
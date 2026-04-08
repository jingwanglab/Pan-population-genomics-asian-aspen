#!/bin/bash
#1. performing pca analysis
smartpca -p pd.txt > pca_out.log
#！！！！！！！！！！！！
#2. using the twstats to select the most significant PCs
twstats -t twtable -i pd.pca.eval -o pd.pca.out


#ragtag.py correct
for i in {ZSP191aL,ZSP192L,ZSP194L,ZSP196L};do ../bin/ragtag.py correct ../ZL2020-SXHZ.fasta ../$i.genome.fasta -o ../ragtag/$i -T corr ;done

#ragtag.py scaffold
for i in {ZSP191aL,ZSP192L,ZSP194L,ZSP196L};do ../bin/ragtag.py scaffold ../ZL2020-SXHZ.fasta ../ragtag/$i/$i.genome.corrected.fasta -o ../ragtag/$i ;done

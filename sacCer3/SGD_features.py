import sys
import os

from collections import defaultdict

types = ["ORF","rRNA_gene","ncRNA_gene","snRNA_gene","tRNA_gene","snoRNA_gene"]
chrom = [ None, "I" , "II", "III", "IV", "V", "VI", "VII","VIII","IX","X",
                "XI","XII","XIII","XIV","XV","XVI","XVII" ]

strand = { "W" : "+", "C" : "-" }

feature_types = set()
parents = defaultdict(set)
#coding = open("sacCer3.coding.bed",'w')
#repeat = open("sacCer3.repeat.bed",'w')
#rna = open("sacCer3.rna.bed",'w')

with open("SGD_features.tab") as f, open("SGD_features.bed",'w') as out:
    for line in f:
        line = line.strip('\n').split('\t')
#        parents[ line[1] ].add( line[6] )
        feature_types.add(  line[1] )
for item in feature_types:
    print item


#        if line[1] in types and line[2] == "Verified":
#            try:
#                my_chr = chrom[int(line[8])]
#            except ValueError:
#                my_chr = line[8]
#            start,end = line[9],line[10]
#            if int(start) > int(end):
#                end,start = start,end
#            outline = [ my_chr, start, end , line[3], "1", strand[line[11]] ]
#            out.write( "\t".join(outline)+"\n")

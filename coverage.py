import sys
import os

from collections import defaultdict
from operator import itemgetter

import numpy as np

# ==========================================================
#  ---- MAIN ----
# ==========================================================

pipeline = "/home/TeamGilbert/GroupData/CLIP_July2016/Pipeline/"
gene_bed = pipeline+"sacCer3/SGD_features.bed"

def main():
    
    hits = sys.argv[1]
    hits = gene_hits(gene_bed,hits)

    with open(pipeline+"sacCer3/Pus1targets.csv") as f:
        f.readline()
        for line in f:
            line = line.strip().split()
            try:
                print hits[ line[0] ]
            except KeyError:
                print "None"
#    for key,value in sorted(hits.items(),key=itemgetter(1)):
#        print key,value

def gene_hits(gene_bed,hitfile):
    features = load_features(gene_bed)
    all_hits = {}
    with open(hitfile) as f:
        while True:
            # Get all hits for a single chromosome
            chrom_hits = get_chromosome_hits(features,f)
            if not chrom_hits:
                break
            for gene,count in chrom_hits.items():
                all_hits[gene] = count
    return adjust_scores(all_hits)


def get_chromosome_hits(features,handle):
    chrom,hits = get_chromosome(handle)
    if not chrom:
        return None

    starts = features[chrom]['start']
    ends = features[chrom]['end']
    names = features[chrom]['name']

    start = 0
    chrom_hits = defaultdict(int)

#    for position,peak in hits:
    for i,name in enumerate(names):
        start,end = starts[i],ends[i]
        for position in range(start,end+1):
            chrom_hits[name] += hits[position]
                
    return chrom_hits

def adjust_scores(hits):
    mean = np.mean(hits.values())
    hits = { key : value/mean 
             for key,value in hits.items() }
    return hits


# =========================================================
#  ---- Get Data ----
# =========================================================

def load_features(gene_bed):
    """Load ORF annotations from bedfile to dictionary"""
    annotations = defaultdict(gene_dict)
    with open(gene_bed) as f:
        for line in f:
            line = line.strip().split()
            if not line:
                continue
            chrom,start,end,gene,_,strand = line
            chrom_anno = annotations[chrom]
            chrom_anno['start'].append( int(start) )
            chrom_anno['end'].append( int(end) )
            chrom_anno['name'].append(gene)
            chrom_anno['strand'].append(strand)
    return annotations


def get_chromosome(handle):
    """Get the coverage for each position in the current
       chromosome in the file"""
    hits = defaultdict(int)
    # The first line tells me what chromosome I'm looking at
    last = handle.tell()
    try:
        chrom = handle.readline().split()[0]
    except IndexError:
        return None,None
    handle.seek(last)
    # Iterate through the file, finding all hits for that chrom
    while True:
        line = handle.readline().strip().split()
        try:
            my_chrom,position,score = line
        except ValueError:
            break
        if my_chrom == chrom:
            if int(score) > 0:
                hits[ int(position) ] = int(score)
        else:
            handle.seek( last )
            break
        last = handle.tell()
    return chrom,hits


def gene_dict():
    return { 'start' : [],
             'end'   : [],
             'strand': [],
             'name'  : [] }



if __name__=="__main__":
    main()

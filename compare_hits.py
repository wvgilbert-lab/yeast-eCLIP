import sys
import os

from operator import itemgetter

from coverage import *

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

# ==========================================================
#  ---- MAIN ----
# ==========================================================

def main():
    pulldown,library,outfile = sys.argv[1:]

    pulldown_hits = gene_hits(gene_bed, pulldown)    
    input_hits    = gene_hits(gene_bed, library )

    pulldown_hits = rank_hits(pulldown_hits)
    input_hits    = rank_hits(input_hits)

    targets = []
    with open(pipeline+"sacCer3/Pus1targets.csv") as f:
        f.readline()
        for line in f:
            line = line.strip().split()
            try:
                targets.append( line[0] )
            except:
                continue

    plot_hits(pulldown_hits,
              input_hits,
              outfile,
              targets)


def rank_hits(hits):
    """Sort the dictionary as a list of tuples"""
    rank = sorted( hits.items(), 
                   reverse=True,
                   key=itemgetter(1) )
    return rank

def plot_hits(pulldown,library,outfile,targets=[]):
    """Plot the rank of each gene in two libraries"""
    sorted_hits = defaultdict(two_tuple)
    for i,hit in enumerate(library):
        sorted_hits[ hit[0] ][0] = i
    for j,hit in enumerate(pulldown):
        sorted_hits[ hit[0] ][1] = j
    gene_names = sorted_hits.keys()
    sorted_hits = [ tuple(sorted_hits[key]) 
                    for key in gene_names ]
    X,Y = zip(*sorted_hits)

    tgt_X,tgt_Y = [],[]
    for target in targets:
        try:
            position = gene_names.index(target)
        except ValueError:
            tgt_X.append( 0 )
            tgt_Y.append( 0 )
            continue
        tgt_X.append(X[position])
        tgt_Y.append(Y[position])


    better,worse = 0,0
    better_tgt,worse_tgt = 0,0

    for x,y in zip(X,Y):
        if y < x:
            better += 1
        if y >= x:
            worse += 1
    for x,y in zip(tgt_X,tgt_Y):
        if y < x:
            better_tgt += 1
        if y >= x:
            worse_tgt += 1
    print better,worse
    print better_tgt,worse_tgt

    sys.exit()
            
    fig,axis = plt.subplots()
    
    axis.plot(X,Y,"k.",alpha=0.3)
    axis.plot([0,6000],[0,6000],
              'b--',alpha=0.7,linewidth=2.5)
    axis.plot(tgt_X,tgt_Y,'r^')

    axis.set_xlabel("input ranks")
    axis.set_ylabel("CLIP ranks")
    
    fig.savefig(outfile)
    

def two_tuple():
    return [None,None]

if __name__=="__main__":
    main()

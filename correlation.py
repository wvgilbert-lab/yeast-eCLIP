import sys
import os
import argparse

import pandas as pd

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import seaborn as sns

from matplotlib_venn import venn2
import numpy as np


annotations = "/home/TeamGilbert/GroupData/CLIP_July2016/Pipeline/sacCer3/SGD_features.chr.bed"

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('rep1', type=str,
                        help='First replicate peaks')
    parser.add_argument('rep2', type=str,
                        help='Second replicate peaks')
    parser.add_argument('--pref', 
                        help='Output prefix')
    args = parser.parse_args()


    replicates = overlap(args.rep1,args.rep2,args.pref)

    data = pd.read_csv(replicates, sep='\t',
                       names=['r1_chr','r1_start','r1_end','r1_p','r1_r','r1_strand',
                              'r2_chr','r2_start','r2_end','r2_p','r2_r','r2_strand',
                              'overlap' ])
    plot_pvals(data,args.pref)
    plot_rvals(data,args.pref)


    gene_overlap(args.rep1,args.rep2,args.pref)


def gene_overlap(rep1,rep2,prefix):
    feats_1 = prefix+".genes1.bed"
    feats_2 = prefix+".genes2.bed"
    # get gene list for rep 1
    command_1 = [ "bedtools", "intersect", "-wa", "-s",
                  "-b", rep1, "-a", annotations, 
                  ">", feats_1  ]
    os.system( " ".join(command_1) )
    # get gene list for rep 2
    command_2 = [ "bedtools", "intersect", "-wa", "-s",
                  "-b", rep2, "-a", annotations, 
                  ">", feats_2  ]
    os.system( " ".join(command_2) )
    # overlap between the two?
    with open( feats_1 ) as f:
        genes_1 = set([ line.split('\t')[3]
                        for line in f ])
    with open( feats_2 ) as f:
        genes_2 = set([ line.split('\t')[3] 
                        for line in f ])
    genes = genes_1 & genes_2
    # write to file
    # plot overlap
    f,ax = plt.subplots()
    venn2( (len(genes_1),len(genes_2),len(genes)),
           set_labels=('Replicate 1','Replicate 2'),
           ax=ax )
    ax.set_title( prefix.split("/")[-1] )
    f.savefig( prefix+".genes.png")
    return 
    # Cleanup
    os.remove( feats_1 )
    os.remote( feats_2 )

def overlap(rep1,rep2,prefix):
    outfile = "overlap.bed"
    if prefix:
        outfile = prefix+"."+outfile

    command = [ "bedtools", "intersect", 
                "-a", rep1, "-b", rep2, "-wo", "-s",
                ">", outfile ]
    os.system( " ".join(command) )
    return outfile


def plot_pvals(data,prefix):
    R = np.corrcoef(data['r1_p'],data['r2_p'])
    f,ax = plt.subplots()
    sns.regplot(data=data,x="r1_p",y="r2_p",ax=ax)
    plt.figtext(0.25,0.75,"R="+str(R[0][1]),
                axes=ax)
    ax.set_xlabel('Replicate 1 peak p-values')
    ax.set_ylabel('Replicate 2 peak p-values')
    f.savefig(prefix+".pvals.png")

def plot_rvals(data,prefix):
    R = np.corrcoef(data['r1_r'],data['r2_r'])
    f,ax = plt.subplots()
    sns.regplot(data=data,x="r1_r",y="r2_r",ax=ax)
    plt.figtext(0.25,0.75,"R="+str(R[0][1]),
                axes=ax)
    ax.set_xlabel('Replicate 1 peaks')
    ax.set_ylabel('Replicate 2 peaks')
    f.savefig(prefix+".fc.png")



if __name__ == "__main__":
    main()

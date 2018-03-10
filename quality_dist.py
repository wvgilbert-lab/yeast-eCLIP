import sys
import os

from Bio import SeqIO
from numpy import median

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

def main():
    infile,outfile = sys.argv[1:]
    scores = get_read_qualities(infile)

    X = [ i for i,n in enumerate(Y) ]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.bar(X,scores)
    ax.set_xlabel('Read quality score (median)')
    ax.set_ylabel('Read counts')
    fig.savefig(outfile)
    
def get_read_qualities(fastq):
    scores_dist = [ 0 for i in range(40) ]
    with open(fastq) as f:
        for record in SeqIO.parse(f,'fastq'):
            scores = record.letter_annotations["phred_quality"] 
            scores_dist[ int(median(scores)) ] += 1
            #median_scores.append( median(scores) ) 
    return scores_dist


if __name__=="__main__":
    main()

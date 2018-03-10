""" Receives the two input-normalized replicates for an eCLIP
    experiment and 
       (1) Visualizes the distribution of peaks
       (2) Enforces a p-value and fold-enrichment threshold
       (3) Visualizes the fraction of stringent peaks 
"""

import sys
import os
import argparse
import csv

import pybedtools as bedtools
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.backends.backend_pdf as pdf_plot
from matplotlib import pyplot as plt
import seaborn as sns

sns.set_style('white')

def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("peaks_1", metavar="rep1.peaks.bed",
                        help="Input normalized peaks for replicate 1")
    parser.add_argument("peaks_2", metavar="rep2.peaks.bed")
    parser.add_argument("prefix",
                        help="prefix for output files, including path")
    parser.add_argument("-p",default=0.00001,type=float,
                        help="p-value threshold")
    parser.add_argument("-r",default=4,type=float,
                        help="enrichment threshold")
    args = parser.parse_args()

    pdf = pdf_plot.PdfPages( args.prefix+".quality.pdf" )
    
    peaks_1 = get_peak_scores(args.peaks_1)
    peaks_2 = get_peak_scores(args.peaks_2)

    ok_peaks = args.prefix+".single_rep_peaks.bed"
    write_single_peaks( peaks_1, peaks_2, ok_peaks,
                        r=args.r,p=args.p )
    
    f1,f2 = peaks_dist(peaks_1,peaks_2,args.r)
    pdf.savefig(f1,width=15,height=5)
    pdf.savefig(f2,width=15,height=10)
    f = threshold(args.peaks_1,args.peaks_2,args.prefix,args.r,args.p)
    pdf.savefig(f)
    
    pdf.close()


def write_single_peaks( peaks_1,peaks_2,outfile, r, p ):
    with open( outfile,'w' ) as out:
        for peaks in [ peaks_1, peaks_2 ]:
            for position,(pval,fc) in peaks.items():
                my_pval = 10**(-pval)
                my_fc   = 2**(fc)
                if my_pval<p and my_fc>r:
                    chrom,start,end,strand = position
                    outline = [ chrom, start, end, str(pval), str(fc), strand ]
                    out.write( "\t".join(outline)+"\n" )
                    
    
def get_peak_scores( peakfile ):
    """ Get the (-log_10 pvalue, log_2 fold-change) for each peak,
        store it in a dictionary """
    peaks = {}
    with open(peakfile) as f:
        for line in f:
            line = line.strip().split('\t')
            if line:
                my_id = ( line[0],line[1],line[2],line[5] ) #{}:{}-{}{})".format(line[0],line[1],line[2],line[5])
                peaks[ my_id ] = ( float(line[3]) , float(line[4]) )
    return peaks

                
def peaks_dist(rep1,rep2,pval):

    # Plot peaks
    peaks_f,ax = plt.subplots(2,sharex=True,sharey=True)
    for i,r in enumerate([rep1,rep2]):
        my_peaks = [ values[1] for key,values in r.items() ]
        sns.distplot( my_peaks, bins=100, kde=False, rug=True, ax=ax[i] )
    #ax[1].set_yscale("log")
    ax[1].set_xlabel("$\log_2$ fold-change (CLIP vs SM-input)")
    ax[0].set_title("Rep1"); ax[1].set_title("Rep2")
    # Plot peaks v pvalues
    
    bins = list( range(-15,15) )
    pvals_f,ax = plt.subplots(1, 2,sharex=True,sharey=True)
    for i,r in enumerate([rep1,rep2]):
        my_pvals,my_peaks = zip(*list(r.values()))
        sns.regplot( x=np.array(my_peaks), y=np.array(my_pvals),
                     fit_reg=False, ax=ax[i])
    ax[0].set_ylabel("$-\log_{10}$ p-value (CLIP vs SM-input)")
    ax[0].set_title("Rep1"); ax[1].set_title("Rep2")
    #f,ax = plt.subplots(2, sharex=True, sharey=True)
    #sns.distplot(peaks1, bins=30, kde=False, rug=True, ax=ax[0])
    #sns.distplot(peaks2, bins=30, kde=False, rug=True, ax=ax[1])
    #ax[0].axvline(2,color='green',alpha=0.6);ax[1].axvline(2,color='green',alpha=0.6)
    #ax[0].axvline(3,color='green',alpha=0.8);ax[1].axvline(3,color='green',alpha=0.8)
    #ax[1].set_yscale("log");
    
    return peaks_f,pvals_f

def load_peaks(bedfile):
    peaks = []
    with open(bedfile) as f:
        reader = csv.DictReader(f,delimiter='\t',
                                fieldnames=['chrom','start','end','peak','pval','strand'])
        peaks = [ float(row['peak']) 
                  for row in reader
                  if float(row['peak']) < 300 ]
    return peaks


def threshold(peaks1,peaks2,prefix,minR,maxP):

    p1 = bedtools.BedTool( peaks1 )
    p2 = bedtools.BedTool( peaks2 )

#    all_peaks  = open( prefix+".all_rep_peaks.bed", 'w' )
    good_peaks = open( prefix+".stringent_peaks.bed", 'w' )
    
    overlap = p1.intersect(p2,wo=True)

    print( overlap )
    scores = []
    for row in overlap:
        score1 = 2**float(row[3]) > minR and 10**(-float(row[4])) < maxP
        score2 = 2**float(row[9]) > minR and 10**(-float(row[10]))< maxP
        try:
            peak_line = [ row[0], row[1], row[8], row[3], row[9], row[5] ]
#            if score1 or score2:
#                all_peaks.write( "\t".join(peak_line)+"\n" )
            if score1 and score2:
                good_peaks.write( "\t".join(peak_line)+"\n" )
        except:
            print( row )
        scores.append( (score1,score2) )

    good_frac = len([ s for s in scores if s[0] and s[1] ])/ float(len(scores))
    ok_frac   = len([ s for s in scores if s[0] or s[1] ]) / float(len(scores))
    rest_frac = 1 - ok_frac
    ok_frac = ok_frac - good_frac

    print( good_frac,ok_frac,rest_frac)
    
    f,ax = plt.subplots()
    ax.pie([ rest_frac,ok_frac,good_frac ],
           colors = ['yellowgreen', 'gold', 'lightskyblue'],
           labels=["rest of CLIP peaks","meets threshold in only one replicate",
                   "significant CLIP peaks"])
    ax.axis('equal')
    return f


if __name__=="__main__":
    main()


import sys
import os
import argparse

from collections import defaultdict

import numpy as np

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import seaborn as sns
sns.set(style='ticks')

annotations = "/home/cassandra/TeamGilbert/GroupData/CLIP_July2016/Pipeline/sacCer3/SGD_features.tab"

# ==========================================================
#  ---- MAIN / Testing ----
# ==========================================================

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('peaks',
                        help='normalized CLIP peaks (bed)')
    parser.add_argument('targets',
                        help='bedfile or structure file with positions'
                             'of targets sites in CDS or genome')
    parser.add_argument('outfile',
                        help="output plot")
    parser.add_argument('--format',choices=['structure','bed','CDS'],
                        default="CDS",
                        help="targets file format. default:CDS")
    parser.add_argument('--background',default=None,
                        help="bedfile with list of regions to use as a background control")
    args = parser.parse_args()


    structure = load_structure(args.targets)
    starts = CDS_to_genomic(structure)
    peaks = load_peaks(args.peaks)
    site_context = site_peaks( starts,peaks )

    pref = args.outfile[:-4]

    for i,site in enumerate(site_context):
        if site:
            start = starts[i]
            struct= structure[i]

            name = "{}_{}".format(struct['gene'],struct['position'])
            seq  = [ n for n in struct['sequence'] ]
            seq[10] = "Y"
            seq = "".join(seq)
            colors = site[100-10:100+50]

            
            call_VARNA(seq,struct['structure'],colors,name+"."+pref+".png")
     
    metaplot(site_context, args.outfile)


# ==========================================================
#  ---- Load Data ----
# ==========================================================

def load_structure(infile):
    data = list()
    with open(infile) as f:
        while True:
            try:
                name = f.readline().strip()[1:]
            except ValueError:
                break
            # Read the sequence
            sequence = f.readline().strip()
            if not sequence:
                break
            # Read line and parse structure and MFE
            structure = f.readline().strip().split()
            structure,MFE = structure[0],structure[-1]
            structure.strip()
            MFE = float(MFE.strip("()"))
            # Parse out gene name 
            gene,pos,_ = name.split("_")
            # Store data in list of dicts
            data.append( { 'gene'     : gene, 
                           'position' : int(pos[1:]),
                           'sequence' : sequence,
                           'structure': structure,
                           'MFE'      : MFE } )

    return data

def load_peaks(peakfile):
    """ Load CLIP peaks from a bedfile into a dictionary 
        of { chr : { position : hits }}"""
    peaks = defaultdict(lambda:defaultdict(int))
    with open(peakfile) as f:
        for line in f:
            chrom,start,end,pval_R,fold_R,strand = line.strip().split()
            for i in range( int(start), int(end)+1):
                peaks[chrom][i] = float(fold_R)
    return peaks


def load_bed(bedfile):
    starts = defaultdict(list)
    ends   = defaultdict(list)
    with open(bedfile) as f:
        for line in f:
            line = line.strip().split('\t')
            chrom,start,end = line[:3]
            starts[chrom].append( int(start) )
            end[chrom].append( int(end) )
    return starts,ends

# ==========================================================
#  Genomic coordinates 
# =========================================================            

def get_CDS_start(genelist):
    """ Receives a list of relative CDS positions and 
        converts it to absolute genomic position """
    chrom = [ None, "I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", 
                   "XI","XII","XIII","XIV","XV","XVI","XVII","XVIII","XIX","XX", ]
    strand = { "W" : 1, "C" : -1}
    CDS_start = {}
    with open(annotations) as f:
        for line in f:
            line = line.strip().split('\t')
            if line[1] == "CDS" and line[6] in genelist:
                my_chr = "chr"+chrom[int(line[8])]
                CDS_start[line[6]] = ( my_chr,
                                       int(line[9]),
                                       strand[line[11]] )
    return CDS_start

def CDS_to_genomic(positions):
    target_positions = []
    genelist = [ record['gene'] for record in positions ]
    CDS_starts = get_CDS_start(genelist)
    for position in positions:
        try:
            chrom,start,strand = CDS_starts[ position['gene'] ]
        except KeyError:
            continue
        target_positions.append( (chrom, strand*position['position']+start, strand) )
    return target_positions


def site_peaks(positions,peaks,r=100):
    site_context = []
    for position in positions:
        chrom,site,strand = position
        my_peaks = [ peaks[chrom][i]
                     for i in range(site-r,site+r+1) ]
        if strand == -1:
            my_peaks = [ p for p in reversed(my_peaks) ]
        if len(set(my_peaks)) == 1:
            site_context.append( None )
        else:
            site_context.append( my_peaks )
    return site_context

# ==========================================================
#  ---- PLOTS ----
# ==========================================================

def metaplot(site_peaks,outfile):

    density = [ 0 for i in range(201) ]

    for peak in site_peaks:
       if peak:
           for i,val in enumerate(peak):
               if val > 0 :
                   density[i] += 1

    f,ax = plt.subplots(2,sharex=True)
    sns.regplot( x=np.array(range(len(density))), 
                 y=np.array(density),
                 fit_reg=False,
                 color='k',ax=ax[0] )
    
    for Y in site_peaks:
        if Y:
            ax[1].plot( range(len(density)), Y)
        #sns.regplot(range(len(density)),Y,ax=ax[1],marker='-')
    ax[0].axvline(100,color="green",linewidth=2)
    ax[1].axvline(100,color="green",linewidth=2)
    f.savefig(outfile)


def call_VARNA(sequence,structure,peaks,outfile):

    colormap = ";".join([ str(p) for p in peaks ])
    command = [ "java", "-cp", 
                "/home/cassandra/Programs/VARNA/VARNAv3-93.jar", 
                "fr.orsay.lri.varna.applications.VARNAcmd", 
                # Sequence and structure
                "-structureDBN", '"{}"'.format(structure),
                "-sequenceDBN", '"{}"'.format(sequence),
                # Color settings
                "-colorMap", '"{}"'.format(colormap),
                "-colorMapStyle green", 
                "-colorMapMax", '"10.0"', "-colorMapMin", '"0.0"',
                "-nsBasesColor", '"#961F43"',
                # Nucleotide styles
                "-bpStyle line", 
                # Output
                "-o", outfile, 
                "-resolution", '"11.0"']
    os.system( " ".join(command) )
     


# ==========================================================
#  Actually run things
# ==========================================================
if __name__=="__main__":
    main()

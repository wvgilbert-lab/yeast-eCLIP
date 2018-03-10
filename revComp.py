""" Take the reverse complement of the sequences
    in a .fastq or .fasta file. """

import sys
import os
import gzip

from Bio import SeqIO
from Bio.Seq import reverse_complement

def revComp(infile,outfile):

    # Figure out format 
    in_gz = infile.endswith(".gz")
    my_format = infile.split(".")    
    if in_gz:
        my_format = my_format[-2]
    else:
        my_format = my_format[-1]

    if in_gz:
        f = gzip.open(infile,'rb')
    else:
        f = open(infile)

    if outfile.endswith(".gz"):
        out = gzip.open(outfile,'wb')
    else:
        out = open(outfile,'w')
        
    
    if my_format == 'fasta':
        for record in SeqIO.parse(f,'fasta'):
            my_seq = str(record.seq.reverse_complement())
            out.write( ">"+record.id+"\n")
            out.write( my_seq + "\n" )

    elif my_format == 'fastq':
        while True:
            header = f.readline()
            seq = f.readline().strip()
            header2 = f.readline()
            quality = f.readline().strip()
            if not seq:
                break
            out.write( header + reverse_complement(seq) +"\n" )
            out.write( header2 + quality[::-1]+"\n")

if __name__=="__main__":
    infile,outfile = sys.argv[1:]
    revComp(infile,outfile)

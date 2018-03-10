import sys
import os
import gzip
import argparse

from collections import defaultdict

# ==============================================================================

def main():

    # Parser for command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--fastq1",
                        help="fastq file for the first read pair")
    parser.add_argument("--fastq2",
                        help="fastq file for the second read pair")
    parser.add_argument("--prefix",
                        help="output file prefix,including path")
    parser.add_argument("-N",default=10,type=int,
                        help="random-mer length")
    parser.add_argument("--metrics",
                        help="output file for randommer counts")
    args = parser.parse_args()
    randomer = args.N

    # Open input files
    f1 = fastq_file_handle( args.fastq1 )
    f2 = fastq_file_handle( args.fastq2 )

    # Open output files 
    out1 = gzip.open( args.prefix+".r1.fastq.gz",'w' )
    out2 = gzip.open( args.prefix+".r2.fastq.gz",'w' )

    # Initialize a counter of randomer instances
    counter = defaultdict( int )
    
    # Iterate over fastq records to reformat reads
    while True:
        r1_record = get_record( f1 )
        r2_record = get_record( f2 )
        if not (r1_record and r2_record):
            break

        # Find and count the randomer 
        n_mer = r2_record[1][:randomer]
        counter[ n_mer ] += 1

        # Reformat both records
        header_1 = "@{}:{}".format( n_mer,
                                    r1_record[0][1:] )
        header_2 = "@{}:{}".format( n_mer,
                                    r2_record[0][1:])
        sequence_2 = r2_record[1][randomer:]
        quality_2 = r2_record[3][randomer:]

        # Write the new records to the output files
        out1.write("\n".join([ header_1,
                               r1_record[1],
                               "+",
                               r2_record[3] ])+"\n")
        out2.write("\n".join([ header_2,
                               sequence_2,
                               "+",
                               quality_2 ])+"\n")

    # Close all open file handles
    for handle in [ f1, f2 , out1, out2 ]:
        handle.close()

    # Output nmer counts to metrics file
    with open( args.metrics,'w' ) as out:
        for nmer,count in counter.items():
            out.write( "{}\t{}\n".format( nmer,count ) )

# ==============================================================================

def fastq_file_handle( infile ):
    """ Receives an input file and checks whether it's gzipped,
        then opens it accordingly and returns the file handle
    """
    if infile.endswith(".gz"):
        f = gzip.open( infile )
    else:
        f = open( infile )
    return f


def get_record( f ):
    """ Get a fastq record from an open file handle.
        Return as tuple (header, sequence, header2, quality). """
    header = f.readline().strip()
    sequence = f.readline().strip()
    header2 = f.readline().strip()
    quality = f.readline().strip()
    if ( header and sequence and header2 and quality ):
        if header.endswith( '/1' ) or header.endswith( '/2' ):
            header = header[:-2]
        return (header , sequence,
                header2, quality ) 
    else:
        return None

# ==============================================================================

if __name__=="__main__":
    main()

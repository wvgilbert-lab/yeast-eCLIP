import sys
import os
import argparse
import gzip
import pysam


# ==========================================================
#  ---- MAIN ----
# ==========================================================

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('reads_1')
    parser.add_argument('reads_2')
    parser.add_argument('outdir')
    parser.add_argument('outfile')
    args = parser.parse_args()

    library_quality( args.reads_1, args.reads_2,
                     args.outdir, args.outfile )

    
def library_quality(reads1, reads2, outdir, outfile):
    with open(outfile,'w') as out:

        header = [ "Pipeline_step", "total_reads", "reads_discarded"]
        out.write( "\t".join(header)+"\n" )

        # --- Initial reads
        total_reads = count_fastq_reads( reads1,reads2 )
        outline = [ "Initial_reads", str(total_reads), "0" ]
        out.write( "\t".join(outline)+"\n" )
        
        # --- cutadapt round1
        trimmed = [ outdir+"/intermediates/"+fname
                    for fname in os.listdir(outdir+"/intermediates")
                    if ".trimmed.fastq" in fname ]
        trimmed_reads = count_fastq_reads( trimmed[0],trimmed[1] )
        outline = [ "cutadapt_1",
                    str(trimmed_reads),
                    str(total_reads-trimmed_reads) ]
        out.write( "\t".join(outline)+"\n")
        
        # --- cutadapt round2
        trimmed = [ outdir+"/intermediates/"+fname
                    for fname in os.listdir(outdir+"/intermediates/")
                    if "round2.fastq" in fname ]
        trimmed_reads_2 = count_fastq_reads( trimmed[0],trimmed[1] )
        outline = [ "cutadapt_2",
                    str(trimmed_reads_2),
                    str(trimmed_reads - trimmed_reads_2) ]
        out.write( "\t".join(outline)+"\n")
        
        # --- Repeat masking
        unmasked = [ outdir+"/mapped/"+fname
                     for fname in os.listdir(outdir+"/mapped/")
                     if "Unmapped.out.mate" in fname and "genome_mapped" not in fname ]
        masked_reads = count_unmapped_reads(unmasked[0],unmasked[1])
        outline = [ "Repeat_masking",
                    str(masked_reads),
                    str(trimmed_reads_2 - masked_reads) ]
        out.write( "\t".join(outline)+"\n" )
        
        # --- Read mapping
        unmasked = [ outdir+"/mapped/"+fname
                     for fname in os.listdir(outdir)
                     if "Unmapped.out.mate" in fname and "genome_mapped" in fname ]
        unmapped_reads = count_unmapped_reads( unmasked[0],unmasked[1] )
        outline = [ "Read_mapping",
                    str(masked_reads - unmapped_reads),
                    str(unmapped_reads) ]
        out.write( "\t".join(outline)+"\n" )
        
        # --- Dereplicating
        derep = [ outdir+fname
                  for fname in outdir
                  if fname.endswith(".rmDup.r2.bam") ]
        command = [ "samtools view -F 0x40", derep[0], "| cut -f1 | sort | uniq | wc -l" ]
#        os.system( " ".join(command) )
        
# ==========================================================
#  ---- Counting functions ----
# ==========================================================

def count_fastq_reads(fastq1,fastq2=None):
    """ Count the number of reads in a single fastq file 
        or in a pair of fastq's for paired-end reads """
    reads1 = count_single_fastq(fastq1)
    if fastq2:
        reads2 = count_single_fastq(fastq2)
        if reads1==reads2:
            return reads1
        else:
            sys.exit("Read pair fastq files don't match")
    return reads1

def count_single_fastq(fastq):
    if fastq.endswith(".gz"):
        f = gzip.open(fastq)
    else:
        f = open(fastq)
    for i,line in enumerate(f):
        continue
    return (i+1)/4


def count_unmapped_reads(mate1,mate2):
    return count_fastq_reads(mate1,mate2)

def count_bam_reads(bamfile):
    """ Count the number of unique reads in a bamfile """
    # Generate the name of the sorted file
    sorted = bamfile.split(".")
    sorted[-1] = ".sorted.bam"
    sorted = ".".join(sorted)
    # If the sorted file does not exist, create it
    if not os.path.isfile( sorted ):
        os.system( " ".join(["samtools","sort",bamfile,"-o",sorted]) )        
    # Create an index file if it doesn't already exist
    index = sorted+".bai"
    if not os.path.isfile( index ):
        os.system( " ".join(["samtools","index",sorted]) )
    # Use samtools idxtools to count the number of reads

# ==========================================================
#  ---- RUN THINGS ----
# ==========================================================

if __name__=="__main__":
    main()

import sys
import os
import argparse

from gscripts.clipseq import barcode_collapse_pe

#printOnly = True
printOnly = False
intermediates = []

project_dir = "/home/TeamGilbert/GroupData/CLIP_July2016"
genome_index = project_dir+"/Pipeline/sacCer3/sacCer3_STAR_index/"
annotations  = project_dir+"/Pipeline/sacCer3/annotation/Saccharomyces_cerevisiae.R64-1-1.77.fullName.gtf"

picard = project_dir+"/Pipeline/picard-2.6.0/picard.jar"

# ==========================================================
#  ---- MAIN ----
# ==========================================================
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('reads_1', metavar='r1.fastq.gz',
                        help="First read file from paired-end sequencing")
    parser.add_argument('reads_2', metavar='r2.fastq.gz',
                        help="Second read file from paired-end sequencing")
    parser.add_argument('outpref', metavar="output_prefix",
                        help="Output prefix, including ")
    parser.add_argument('--n',type=int,default=10,
                        help="length of random n-mer (default: 10)")
    parser.add_argument('--adapt5', nargs='+',
                        default=["TTACTATGCCGCTGGTGGCTCTAGATGTGCAAGTCTCAAGATGTCAGGCT","GATGTGCAAGTCTCAAGATGTCAGGCT"],
                        help="list of 5' adapters")
    parser.add_argument('--adapt3', nargs='+',
                        default=["GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"],
                        help="list of 3' adapters")
    args = parser.parse_args()

    logfile = args.outpref+".log"
    
    # Two rounds of adapter trimming, first on 5' and 3',
    # second round just 3'
    with open(logfile,'a') as out:
        out.write(" ========================== \n")
        out.write("      Adapter trimming \n ")
        out.write(" ========================== \n\n" )

    trim1,trim2 = trim_adapters(args.reads_1,args.reads_2,args.outpref,
                                adapt_5 = args.adapt5,
                                adapt_3 = args.adapt3)

    trim1,trim2 = trim_adapters(trim1,trim2,args.outpref+".trimmed",
                                adapt_3=args.adapt3)
                                
    # Map reads to the genome without removing repeats
    # Criteria should be fairly lax
    with open(logfile,'a') as out:
        out.write("\n ====================== \n")
        out.write("      Read mapping \n ")
        out.write(" ====================== \n\n" )
    mapped = map_reads(trim1,trim2,args.outpref)

    # Use picard to determine the insert length distribution
    # before and after collapsing duplicates
    outfile = args.outpref+".sorted.bam"
    sort_sam(mapped,outfile,args.outpref)
    insert_length_dist(outfile,args.outpref)

    # Collapse barcodes
    derep = args.outpref+".rmDup.bam"
    if not printOnly:
        barcode_collapse_pe.barcode_collapse(mapped,derep)

    insert_length_dist(derep,args.outpref+".rmDup")



# ==========================================================
#  ---- Functions ----
# ==========================================================
    
def trim_adapters(r1,r2,outpref,
                  adapt_5=[],adapt_3=[],log=None):
    """ Use cutadapt to remove the adapters from the read 
        pairs """
    # First generate all of the options for the command
    command = [ "cutadapt", "-f", "fastq", 
                "--match-read-wildcards",  
                "--times","1" ,"-e","0.1", "-O","1",
                "--quality-cutoff", "6",  "-m 18"]

    # Then add all of the adapters to trim
    for adapt in adapt_3:
        command.extend(["-A",adapt])
    for adapt in adapt_5:
        command.extend(["-g",adapt])

    # Finally add all the input and output filenames
    out1 = outpref+ ".1.fastq.gz"
    out2 = outpref+ ".2.fastq.gz"
    # Add the filenames to a list of intermediate files
    intermediates.extend([ out1,out2 ])

    # Add the rest of the command options
    command.extend([ "-o" , out1, "-p", out2, 
                     r1, r2 ,
                     ">>" , outpref+".log" ])
    if log:
        command[-1] = log
    # Run the command
    run_command( command )
    return out1,out2


def map_reads(read1,read2,outpref):
    outfile = outpref+".mapped.bam"
    intermediates.append( outfile )
    command = [ "STAR",  
                "--runMode alignReads",  
                "--runThreadN 4",  
                "--genomeDir", genome_index, 
                "--genomeLoad NoSharedMemory",
                "--readFilesIn", read1, read2,
                "--readFilesCommand", "zcat",
                "--outSAMunmapped Within",  
                "--outFilterMultimapNmax 20",  
                "--outFileNamePrefix", outfile,
                "--outSAMattributes All",  
                "--outStd BAM_Unsorted",  
                "--outSAMtype BAM Unsorted",  
                "--outFilterType BySJout",  
                "--outReadsUnmapped Fastx",  
                "--outFilterScoreMin 10",  
                "--outSAMattrRGline ID:foo",  
                "--alignEndsType EndToEnd", 
                ">", outpref+".bam" ] 
    run_command( command )
    return outfile

def sort_sam(bamfile,outfile,outpref):
    # Sort a given bamfile by 
    command = [ "java -jar", picard, "SortSam",
                "INPUT=" + bamfile,
                "OUTPUT=" + outfile,
                "SO=coordinate", 
                "VALIDATION_STRINGENCY=SILENT",
                ">>", outpref+".log" ]
    run_command(command )
    # Create an index for the new bamfile
    command = [ "samtools index", outfile ]
    run_command( command )
    return outfile

def insert_length_dist(mapped_reads,outpref):
    command = [ "java -jar", picard, "CollectInsertSizeMetrics",
                "INPUT="+mapped_reads,
                "OUTPUT="+outpref+".sizes.txt",
                "HISTOGRAM_FILE="+outpref+".sizes.pdf"]
    run_command( command )
    return

def run_command( command ):
    if printOnly:
        print( " ".join(command)+"\n" )
    else:
        print( " ".join(command)+"\n" )
        os.system( " ".join(command) )


def cleanup( outpref ):
    for fname in intermediates:
        os.remove( fname )



# ==========================================================
#  ---- Actually run things ----
# ==========================================================

if __name__=="__main__":
    main()

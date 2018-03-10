import sys
import os
import argparse

from Bio import SeqIO

from revComp import revComp

repeat_index = "/home/TeamGilbert/GroupData/CLIP_July2016/Pipeline/sacCer3/sacCer3_repeats/"
genome_index = "/home/TeamGilbert/GroupData/CLIP_July2016/Pipeline/sacCer3/sacCer3_STAR_index/"


# -----------------
#       MAIN
# -----------------

def main():

    # Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('reads_1', metavar='r1.fastq.gz',
                        help="First read file from paired-end sequencing")
    parser.add_argument('reads_2', metavar='r2.fastq.gz',
                        help="Second read file from paired-end sequencing")
    parser.add_argument('outdir', metavar="output_dir",
                        help="Output directory")
    parser.add_argument('--n',type=int,default=10,
                        help="length of random n-mer (default: 10)")
    parser.add_argument('--adapt5', nargs='+',
                        help="list of 5' adapters")
    parser.add_argument('--adapt3', nargs='+',
                        help="list of 3' adapters")
    parser.add_argument("--strand", default="both"
                        help="What read pair are we mapping?"
                        "'both','forward', or 'reverse'")
    args = parser.parse_args()


    # Files, filenames, directories
    make_output_directory(args.outdir)

    # Fastqc round 1
    fastqc(args.reads_1,args.outdir)
    fastqc(args.reads_2,args.outdir)

    # First round of cutadapt
    trim1,trim2 = cutadapt(args.reads_1, args.reads_2, args.outdir,
                       adapt_5=args.adapt5, adapt_3=args.adapt3 )
    # Second round of cutadapt deals with double-ligation on 3'
    trim1,trim2 = cutadapt(trim1, trim2, args.outdir, 
                           adapt_3=args.adapt3, suffix="round2")


    if args.strand == "both":
        pass
    elif args.strand == "forward":
        STAR_map_reads(trim1,args.outdir)
    elif args.strand == "reverse":
        trim2_rc = trim2.rsplit(".",1)[0]+".revComp.fastq"
        revComp(trim2,trim2_rc)
        STAR_map_reads(trim2_rc,args.outdir) 
    
   


# ----------------------------------------------------------
#       Wrapper functions for each step in the pipeline 
# ----------------------------------------------------------

def fastqc(infile,outdir):
    prefix = get_prefix(infile)
    outfile = outdir+"/logs/"+prefix+".fastqc.log"
    command = [ "fastqc", infile, "-q", 
                "-o" , outdir+"/fastqc/" ,
                ">>" , outfile ]
    print( " ".join(command)+"\n" )
    run_command( command )
    return outfile


def cutadapt(r1,r2,outdir,adapt_5=[],adapt_3=[],suffix=None):
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
    prefix1 = get_prefix(r1)
    prefix2 = get_prefix(r2)
    if not suffix:
        suffix = "trimmed"
    out1 = "{}/intermediates/{}.{}.fastq.gz".format(outdir,
                                                    prefix1,
                                                    suffix)
    out2 = "{}/intermediates/{}.{}.fastq.gz".format(outdir,
                                                    prefix2,
                                                    suffix)
    command.extend([ "-o" , out1, "-p", out2, 
                     r1, r2 ,
                     ">>" , outdir+"/logs/cutadapt.log" ])
    print( " ".join(command)+"\n" )
    run_command( command )
    # Run the command
    return out1,out2

def STAR_map_reads(read1,outdir):
    prefix1 = get_prefix(read1)
    command = [ "STAR", 
                "--genomeDir", genome_index, 
                "--readFilesIn", read1, 
                "--outFileNamePrefix", outdir+"/mapped/"+prefix1,
                ">>", outdir+"/logs/star.repeats.txt" ]
#   command = [ "STAR", "--runMode", "alignReads",
 #               "--runThreadN 8",
 #               "--genomeDir", genome_index,
 #               "--genomeLoad", "NoSharedMemory",  
 #               "--readFilesIn", read1,
 #               "--readFilesCommand", "gunzip",
            # Include unmapped reads
 #               "--outSAMunmapped Within",
            # Limit multiple mapping
 #               "--outFilterMultimapNmax 1",  
 #               "--outFilterMultimapScoreRange 1",
            # Output file options
 #               "--outFileNamePrefix","{}/mapped/{}".format(outdir,
  #                                                          prefix1),
            # Output as unsorted BAM to stdout
   #             "--outSAMattributes All",  "--outStd BAM_Unsorted",  
   #             "--outSAMtype BAM Unsorted", "--outFilterType BySJout",  
   #             "--outReadsUnmapped Fastx", "--outFilterScoreMin 10",
   #             "--outSAMattrRGline ID:foo", "--alignEndsType EndToEnd", 
   #             ">", "{}/mapped/{}.bam".format(outdir,prefix1) ]
    print( " ".join(command)+"\n" )
    run_command( command )


def STAR_map_reads(read1,read2,outdir):
    prefix1 = get_prefix(read1)
    command = [ "STAR", 
                "--genomeDir", genome_index, 
                "--readFilesIn", read1, read2,
                "--outFileNamePrefix", outdir+"/mapped/"+prefix1,
                ">>", outdir+"/logs/star.repeats.txt" ]
    run_command(" ".join(command)+"\n")

def run_command(command_list):
    command = " ".join(command_list)
    os.system(command)

# ---------------------------
#       Other functions
# ---------------------------

def make_output_directory(outdir):
    """Make sure the directory exits 
    and make subdirectories"""
    # Outdir exits
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    # Make subdirs
    for subdir in ["fastqc","mapped","logs","intermediates"]:
        this_dir = "{}/{}".format(outdir,subdir) 
        if not os.path.isdir( this_dir ):
            os.makedirs(this_dir)

def get_prefix(fname):
    """Remove path and suffixes from filename
       to get a prefix """
    fname = fname.split("/")[-1]
    if fname.endswith(".gz"):
        prefix = fname.split(".")[:-2]
        return ".".join(prefix)
    else:
        return fname.rsplit(".",1)[0]
    

def count_blank(fastq):
    blank = 0
    with open(fastq) as f:
        for i,record in enumerate(SeqIO.parse(f,'fastq')):
            my_seq = str( record.seq )
            if len(my_seq) == my_seq.count("N") + my_seq.count("G"):
                blank += 1
    return float(blank)/i+1

# --------------------------------
#       Actually run things 
# --------------------------------

if __name__=="__main__":
    main()

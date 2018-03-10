import sys
import os
import argparse

import multiprocessing

#from library_quality import library_quality
#print_only = True

#print_only = False
print_only = True

project_dir = "/home/ubuntu/"

genome_index = project_dir+"/Pipeline/sacCer3/sacCer3_STAR_index/"
annotations  = project_dir+"/Pipeline/sacCer3/annotation/Saccharomyces_cerevisiae.R64-1-1.77.fullName.gtf"

gscripts = project_dir+"Pipeline/gscripts-1.1/gscripts/"
count_aligned = gscripts+"/general/count_aligned_from_sam.py"
fix = gscripts+"/clipseq/fix_scores.py"
picard = project_dir+"/picard/picard.jar"

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
    parser.add_argument('--species', default='sacCer3',
                        help="reference genome to map to. (sacCer3 or Sigma1278b)")
    parser.add_argument('--prefix', 
                        help="prefix for output files")
    args = parser.parse_args()


    # Files, filenames, directories
    make_output_directory(args.outdir)

    reads1,reads2 = demultiplex( args.reads_1,args.reads_2,
                                 args.outdir,args.prefix )

    # Fastqc round 1
    fastqc(reads1,args.outdir)
    fastqc(reads2,args.outdir)

    # First round of cutadapt
    trim1,trim2 = cutadapt(reads1,reads2, args.outdir,
                           adapt_5=args.adapt5, adapt_3=args.adapt3,
                           cut_round=1)
    # Second round of cutadapt deals with double-ligation on 3'
    trim1,trim2 = cutadapt(trim1, trim2, args.outdir, 
                           adapt_3=args.adapt3,
                           cut_round=2,suffix="round2")
    
    # Fastqc round 2
    fastqc(trim1,args.outdir)
    fastqc(trim2,args.outdir)
    # Now sort and map the reads 
#    mate1,mate2 = sort_reads(trim1,trim2,args.outdir)

    mapped = STAR_map_reads(trim1,trim2,args.outdir+"/"+args.prefix+".genome_mapped")
    unique = args.outdir+"/"+args.prefix+".genome_mapped.unique.bam"
    get_unique(mapped,unique)

    # Collapse barcodes
    derep = args.outdir+"/mapped/genome_mapped.rmDup.bam"
    metrics = args.outdir+"/logs/rmDup.log"
    collapse(unique,derep,metrics)
    # unique replaced mapped

    # Sort the mapped reads by coordinate
    mapped = args.outdir+"/"+args.prefix+".rmDup.bam"
    sorted_reads = sort_sam(derep,mapped,args.outdir)
    
    sorted_read2 = samtools_view(sorted_reads)

    peaks_pref = sorted_read2.rsplit(".",1)[0]    
    peaks = clipper(sorted_read2,peaks_pref)
    fix_scores(peaks,args.outdir+"/"+args.prefix)
    


# ----------------------------------------------------------
#       Wrapper functions for each step in the pipeline 
# ----------------------------------------------------------
def demultiplex(fastq1,fastq2,outdir,prefix):
    command = [ "python", project_dir+"/Pipeline/demultiplex.py",
                "--fastq1", fastq1, "--fastq2", fastq2,
                "--prefix", outdir+"/intermediates/"+prefix, 
                "--metrics", outdir+"/logs/demux.metrics.txt" ]
    run_command( command )
    return (outdir+"/intermediates/"+prefix+".r1.fastq.gz",
            outdir+"/intermediates/"+prefix+".r2.fastq.gz")

def fastqc(infile,outdir):
    prefix = get_prefix(infile)
    outfile = outdir+"/logs/"+prefix+".fastqc.log"
    command = [ "fastqc", infile, "-q", 
                "-o" , outdir+"/fastqc/" ,
                ">>" , outfile ]
    run_command( command )
    return outfile


def cutadapt(r1,r2,outdir,adapt_5=[],adapt_3=[],suffix=None,cut_round=1):
    # First generate all of the options for the command
    if cut_round == 1:
        overlap = 1
    else:
        overlap = 5
    overlap = str(overlap)
    command = [ "cutadapt", "-f", "fastq", 
                "--match-read-wildcards",  
                "--times","1" ,"-e","0.1", "-O",overlap,
                "--quality-cutoff", "6",  "-m 18"]
    # Then add all of the adapters to trim
    if cut_round == 1:
        command.extend([ "-a", 'NNNNNGATCGTCGGACTGTAGAACTCTGAACGTGTAG',
                         "-g", "GCACCCGAGAATTCCA",                     
                         "-g", "GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
                         "-g", "CAAGCAGAAGACGGCATACGAGATAACTTG" ])
    command.extend([ "-a", 'TGGAATTCTCGGG',
                     "-a", 'CCAAGGAACTCCAGTCAC',
                     "-a", 'TGGAATTCTCGGGTGCCAAGGCACTCCAGTCACCAAGTTATCTCGTATGCCGTCTTCTGCTTG'])

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
    run_command( command )
    # Run the command
    return out1,out2


def count_aligned_reads(bamfile,outdir):
    command = [ "samtools", "view", bamfile,
                "|", 
                count_aligned,
                ">", outdir+"/mapped/repeat_metrics.txt" ]
    run_command( command )

                
def sort_reads(mate1,mate2,outdir):
    # Unzip the two input fastqs in parallel
    run_command( ["gunzip", mate1, "&"] )
    run_command( ["gunzip", mate2, "&"] )
    run_command( ["wait"] )
    # Renome .gz from stored filename
    mate1 = mate1.rsplit(".",1)[0]
    mate2 = mate2.rsplit(".",1)[0]
    # Now sort the fastqs
    sorted_mate1 = outdir+"/intermediates/sorted.mate1.fastq"
    sorted_mate2 = outdir+"/intermediates/sorted.mate2.fastq"
    command = [ "fastq-sort", "--id", mate1, 
                ">", sorted_mate1,
                "&& ",
                "fastq-sort", "--id", mate2, 
                ">", sorted_mate2 ]
    run_command( command )
    return sorted_mate1,sorted_mate2


def STAR_map_reads(read1,read2,prefix):
    prefix1 = get_prefix(read1)
    command = [ "STAR",  
                # Job parameters
                "--runMode alignReads",  
                "--runThreadN",str(multiprocessing.cpu_count()),  
                "--outFilterMismatchNoverLmax 0.09",
                # Reference genome
                "--genomeDir", genome_index, 
                "--genomeLoad NoSharedMemory",
                # Input files
                "--readFilesIn", read1, read2,
                "--readFilesCommand", "zcat",
                # Output format
                "--outSAMunmapped Within",  
                "--outFilterMultimapNmax 20",  
                "--outFilterMultimapScoreRange 1",  
                "--outFileNamePrefix", prefix+".bam", 
                "--outSAMattributes All",  
                "--outStd BAM_Unsorted",  
                "--outSAMtype BAM Unsorted",  
                "--outFilterType BySJout",  
                "--outReadsUnmapped Fastx",  
                "--outFilterScoreMin 10",  
                "--outSAMattrRGline ID:foo",  
                "--alignEndsType EndToEnd", 
                ">", prefix+".bam" ] 
    run_command( command )
    return prefix+".bam"

def get_unique(mapped,outfile):
    command = [ "samtools view",
                "-q", "200",
                "-hb", mapped,
                ">", outfile]
    run_command( command )

def collapse(mapped,derep,metrics):
    command = [ "python", gscripts+"clipseq/barcode_collapse_pe.py",
                "--bam", mapped,
                "--out_file", derep,
                "--metrics_file", metrics ]
    run_command( command )
    return derep

def sort_sam(derep_bam,outfile,outdir):
    # Sort a given bamfile by 
#    outfile = derep_bam.rsplit(".",1)[0]+".sorted.bam"
    command = [ "java -jar", picard, "SortSam",
                "INPUT=" + derep_bam,
                "OUTPUT=" + outfile,
                "SO=coordinate", 
                "VALIDATION_STRINGENCY=SILENT",
                ">", outdir+"/logs/SortSam.log" ]
    run_command(command )
    # Create an index for the new bamfile
    command = [ "samtools index", outfile ]
    run_command( command )
    return outfile

def samtools_view(infile):
    outfile = infile.rsplit(".",1)[0]+".r2.bam"
    command = [ "samtools view", 
                "-hb", 
                "-f 128",  
                infile,
                ">", outfile ]
    run_command( command )
    return outfile

def clipper(infile,prefix):
    outbed = prefix+".peaks.bed"
    command = [ "clipper",
                # Inputs
                "-b", infile, 
                "-s", 'sacCer3',
                # Output
                "-o", outbed ]
    run_command( command )
    return outbed 

def fix_scores(infile,prefix):
    command = [ "python", fix,
                "--bed", infile,
                "--out_file", prefix+".peaks.fixed.bed"]
    run_command( command )

def bigBed(infile,outdir):
    command = [ "bedToBigBed", infile, 
                chrom_size,
                outdir+"/peaks.fixed.bb",
                "-type=bed6+4"]
    run_command( command )



# ---------------------------
#       Other functions
# ---------------------------

def run_command(command_list):
    command = " ".join(command_list)
    if not print_only:
        os.system(command)
    print( command+"\n" )

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

    
def cleanup(outdir):
    for fname in os.listdir(outdir+"/mapped/"):
        pass


def rev_comp(sequence):
    comp = {"A":"T","T":"A",
            "C":"G","G":"C",
            "N":"N"}
    rc = [ comp[nuc] for nuc in reversed(sequence) ]
    return "".join( rc )

# --------------------------------
#       Actually run things 
# --------------------------------


if __name__=="__main__":
    main()

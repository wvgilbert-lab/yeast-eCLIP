import sys, os
import csv

projectdir = "/nobackup1/axolotl/Yeast_eCLIP/"

def make_slurmfile( infile1, infile2, outdir, prefix):
    options = { "job-name"   : prefix+".eCLIP" ,
                "output"     : prefix+".log.txt",
                "export"     : "ALL",
                "mem-per-cpu": "8000",
                "nodes"      : "16",
                "partition"  : "newnodes,sched_mit_hill", 
                "time"       : "12:00:000",
                "mail-type"  : "ALL",
                "mail-user"  : "axolotl@mit.edu" }
    
    command = [ "python", projectdir+"/Pipeline/runPipeline.py",
                infile1, infile2, outdir, 
                "--prefix", prefix ,
                ">", projectdir+"batch_scripts/"+prefix+".sh"]

    slurmfile = projectdir+"batch_scripts/"+prefix+".submit.sh"

    with open(slurmfile,'w') as out:
        out.write("#!/bin/bash\n#\n")
    
        for option,value in options.items():
            out.write( "#SBATCH --{}={}\n".format(option,value) )

        out.write("\n#use the python 2.7 environment for this job \n")
        out.write("setpkgs -a anaconda3\n")
        out.write("source activate py27\n")

        out.write("\n#Run the pipeline\n")
        out.write( " ".join(command)+"\n" )
        out.write( "bash "+projectdir+"batch_scripts/"+prefix+".sh" )
    return slurmfile


if __name__=="__main__":
    
    infile = sys.argv[1]
    batchfiles = []

    with open(infile) as f:
        for row in csv.DictReader(f,delimiter=','):
            batch = make_slurmfile( projectdir+"raw_reads/"+row['read1'],
                                    projectdir+"raw_reads/"+row['read2'],
                                    row['outdir'],
                                    row['prefix'] )
            batchfiles.append( batch )

    for batch in batchfiles:
        os.system( "sbatch "+batch )
            

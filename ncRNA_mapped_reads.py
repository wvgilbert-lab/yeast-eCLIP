import sys
import os

from collections import defaultdict

# ==========================================================
#  ---- MAIN ----
# ==========================================================

def main():
    infile,outfile = sys.argv[1:]
    maps = count_reads( infile )

    with open(outfile,'w') as out:
        out.write("# ======================\n"+
                  "#        SUMMARY\n" +
                  "# ======================\n")
        out.write("tRNA-mapping reads: "+str(sum([ v 
                                                    for k,v in maps.items() 
                                                    if k[0]=="t"])) + 
                  "\n")
        out.write( "rRNA-mapping reads: "+str(sum([ v
                                                    for k,v in maps.items()
                                                    if k.startswith("RDN")]))+ 
                   "\n")
        out.write("other ncRNA: "+str(sum([ v for k,v in maps.items()
                                            if k[0]!="t" and not k.startswith("RDN")]))+
                  "\n\n")
        out.write("# ======================\n"+
                  "#       Detailed\n" +
                  "# ======================\n")
        for key,value in maps.items():
            out.write("{}\t{}\n".format(key,value))
    
def count_reads(infile):

    mapped = defaultdict(int)

    with open(infile) as f:
        prev_name = None
        targets = []
        while True:
            # read a line,break if EOF
            line = f.readline().strip()
            if not line:
                break
            # Skip to next line if it's a header
            if line[0] == "@":
                continue
            # Otherwise, it's an alignment line
            else:
                # Get the name of the read and the target sequence
                line = line.split()
                query = line[0]
                reference = line[2]
                # If it's the first alignment with this read,
                # store the information for the previous read and 
                # start over with this one
                if query != prev_name:
                    for target in targets:
                        mapped[ target ] += float(1)/len(targets)
                    prev_name = query
                    targets = [ reference ]
                else:
                    targets.append( reference )

    return mapped


if __name__=="__main__":
    main()

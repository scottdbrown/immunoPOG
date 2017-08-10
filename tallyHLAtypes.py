'''
Tally HLA types

Input:
    - directory containing directories of OptiType results. Directory names POG###_P####
    - File with POG mutation info (POG ID \t SNV \t INDEL)
Output:
    - For each POGid, the HLA type which has the highest read support

Date: August 10, 2017
@author: sbrown
'''

## Import Libraries
import sys
import argparse
import os

DEBUG = False
VERB = False


if __name__ == "__main__":

    ## Deal with command line arguments
    parser = argparse.ArgumentParser(description = "Tally HLA types")
    parser.add_argument("result_directory", help = "Directory containing the results", type = str)
    parser.add_argument("pog_mutation_fof", help = "File with POG mutation files", type = str)
    parser.add_argument("pog_output_fof", help = "File to write, will be input for 'prepPeptideBindingPredictions...'", type = str)
    parser.add_argument("-d", "--debug", action = "store_true", dest = "DEBUG", help = "Flag for setting debug/test state.")
    parser.add_argument("-v", "--verbose", action = "store_true", dest = "VERB", help = "Flag for setting verbose output.")
    args = parser.parse_args()

    ## Set Global Vars
    DEBUG = args.DEBUG
    VERB = args.VERB

    res = {}

    ## Go through directories, get POG ids, then go through those directories and get result files.
    for root, dirs, files in os.walk(args.result_directory):
        for d in dirs:
            if d.startswith("POG"):
                pid = d.split("_")[0]
                if pid not in res:
                    res[pid] = [0,""]
                for droot, ddirs, dfiles in os.walk(os.path.join(root, d)):
                    for f in dfiles:
                        if f.endswith("result.tsv"):
                            HEADER = True
                            for line in open(os.path.join(droot, f), "r"):
                                if HEADER:
                                    HEADER = False
                                else:
                                    if line.startswith("0"):
                                        ## only look at the first entry
                                        line = line.split("\t")
                                        hlas = ",HLA-".join(line[1:7])
                                        hlas = hlas.replace("*","")
                                        hlas = "HLA-" + hlas
                                        reads = float(line[7])

                                        ## replace what is already stored if this has more read support
                                        if reads > res[pid][0]:
                                            res[pid] = [reads, hlas]

    ##
    out = open(args.pog_output_fof, "w")
    for line in open(args.pog_mutation_fof, "r"):
        line = line.rstrip().split("\t")
        pid = line[0].split("-")[0]
        if pid in res and res[pid][0] != 0:
            ## have HLA data
            out.write("{}\n".format("\t".join([line[0], res[pid][1], line[1], line[2]])))
        else:
            print("No HLA data for {}. Skipping.".format(pid))
    out.close()

    print("done.")

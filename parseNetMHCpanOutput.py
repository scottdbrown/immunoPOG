'''
Parse NetMHCpan Output

Input:
    - Directory of NetMHCpan result files/directories
    - Minimal peptide reference file
Output:
    - Summary file of all pMHC containing mutations

Date: August 9, 2017
@author: sbrown
'''

## Import Libraries
import sys
import argparse
import os

DEBUG = False
VERB = False

resI = {}       ## hold results
pepInfo = {}    ## hold submitted 21(ish)mer info.


def parsePeptideInfo(filename):
    HEADER = True
    for line in open(filename, "r"):
        if HEADER:
            HEADER = False
        else:
            tumor, minPeptideID, ENST, mutation, varPos, ref_mut, seq, refProtSeq = line.rstrip().split("\t")
            pepInfo[int(minPeptideID)] = [ref_mut, ENST, mutation, varPos, refProtSeq]
            

def parseClassI(filename):
    sample = filename.split("/")[-2]    ## second to last is the sample name folder

    if sample not in resI:
        resI[sample] = {}

    IN_PREDICTIONS = False
    for line in open(filename, "r"):
        if IN_PREDICTIONS:
            if line.strip() == "":
                ## finished this chunk
                IN_PREDICTIONS = False
            elif not line.startswith("---"):    ## visual formatting line of output file.
                try:
                    line = line.strip().rstrip().split()
                    ## extract columns of interest.
                    pos = int(line[0])
                    hla = line[1]
                    peptide = line[2]
                    # core = line[3]
                    # Of = line[4]
                    # Gp = line[5]
                    # Gl = line[6]
                    # Ip = line[7]
                    # Il = line[8]
                    # Icore = line[9]
                    identity = int(line[10])
                    # score = line[11]
                    ic50 = float(line[12])
                    # rank = line[13]
                except:
                    print("Error reading line")
                    print("File: {}".format(filename))
                    print("Line: {}".format(line))
                    sys.exit()

                ## replace X in peptide with U
                peptide = peptide.replace("X","U")


                ## Look up identity.
                ref_mut, ENST, mut, varPos, refProtSeq = pepInfo[identity]
                
                ## make spot for data
                if ENST not in resI[sample]:
                    resI[sample][ENST] = {}
                if mut not in resI[sample][ENST]:
                    resI[sample][ENST][mut] = {}
                if hla not in resI[sample][ENST][mut]:
                    resI[sample][ENST][mut][hla] = {}


                ## make sure peptide is not present in the wildtype sequence
                if peptide not in refProtSeq:


                    ## have variant position in submitted peptide, and short peptide position. Need to calculate if mutant position is in each short peptide.

                    if mut.split("-")[0] == "SNV":
                        pepVarPos = int(varPos) - pos + 1
                        pepLen = len(peptide)

                        if pepVarPos > 0 and pepVarPos <= pepLen:
                            ## variant is within peptide.

                            ## check that amino acid is correct (control for Selenocysteine U/Sec converting to X)
                            if (ref_mut == "ref" and peptide[pepVarPos-1] != mut[0] and mut[0] != "U") or (ref_mut == "mut" and peptide[pepVarPos-1] != mut[-1] and mut[-1] != "U"):
                                print("Mutation does not match peptide!")
                                print("NetMHC output: {}".format(line))
                                print("PepInfo: {}".format(pepInfo[identity]))
                                sys.exit()

                            if pos not in resI[sample][ENST][mut][hla]:
                                resI[sample][ENST][mut][hla][pos] = {}
                            if pepLen not in resI[sample][ENST][mut][hla][pos]:
                                resI[sample][ENST][mut][hla][pos][pepLen] = {}
                            resI[sample][ENST][mut][hla][pos][pepLen][ref_mut] = [peptide, pepVarPos, ic50]

                    elif mut.split("-")[0] == "INDEL":
                        if varPos.endswith("*"):
                            ## frameshift.
                            pepVarPos = int(varPos.split("-")[0]) - pos + 1
                            pepLen = len(peptide)

                            if pepVarPos <= pepLen:
                                ## peptide contains variant site or is downstream of variant
                                if pos not in resI[sample][ENST][mut][hla]:
                                    resI[sample][ENST][mut][hla][pos] = {}
                                if pepLen not in resI[sample][ENST][mut][hla][pos]:
                                    resI[sample][ENST][mut][hla][pos][pepLen] = {}
                                resI[sample][ENST][mut][hla][pos][pepLen][ref_mut] = [peptide, pepVarPos, ic50]

                        else:
                            ## no frameshift
                            if "ins" in mut.split("-")[2]:
                                ## in frame insertion
                                pepVarPos1 = int(varPos.split("-")[0]) - pos + 1
                                pepVarPos2 = int(varPos.split("-")[1]) - pos + 1

                                if pepVarPos1 <= pepLen and pepVarPos2 > 0:
                                    ## peptide contains at least one edge of the insertion
                                    if pos not in resI[sample][ENST][mut][hla]:
                                        resI[sample][ENST][mut][hla][pos] = {}
                                    if pepLen not in resI[sample][ENST][mut][hla][pos]:
                                        resI[sample][ENST][mut][hla][pos][pepLen] = {}
                                    resI[sample][ENST][mut][hla][pos][pepLen][ref_mut] = [peptide, pepVarPos1, ic50]

                            elif "del" in mut.split("-")[2]:
                                ## in frame deletion
                                pepVarPos = float(varPos) - pos + 1

                                if pepVarPos <= pepLen and pepVarPos >= 1:
                                    ## peptide contains the novel junction site made by the deletion.
                                    if pos not in resI[sample][ENST][mut][hla]:
                                        resI[sample][ENST][mut][hla][pos] = {}
                                    if pepLen not in resI[sample][ENST][mut][hla][pos]:
                                        resI[sample][ENST][mut][hla][pos][pepLen] = {}
                                    resI[sample][ENST][mut][hla][pos][pepLen][ref_mut] = [peptide, pepVarPos, ic50]

                            else:
                                sys.exit("Non-handled mutation type. Unexpected. {}".format(line))



        elif line.strip().startswith("Pos"):
            ## now in the section of the file that has predictions
            IN_PREDICTIONS = True

def writeOutput(out):
    for samp in resI:
        for ENST in resI[samp]:
            for mut in resI[samp][ENST]:
                for hla in resI[samp][ENST][mut]:
                    for pos in resI[samp][ENST][mut][hla]:
                        for pepLen in resI[samp][ENST][mut][hla][pos]:
                            mutPepSeq = resI[samp][ENST][mut][hla][pos][pepLen]["mut"][0]
                            pepVarPos = resI[samp][ENST][mut][hla][pos][pepLen]["mut"][1]
                            mutIC50 = resI[samp][ENST][mut][hla][pos][pepLen]["mut"][2]

                            out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(samp, ENST, mut, hla, mutPepSeq, pepVarPos, mutIC50))

if __name__ == "__main__":

    ## Deal with command line arguments
    parser = argparse.ArgumentParser(description = "Parse NetMHCpan Output")
    parser.add_argument("--peptideInfoFiles", help = "File(s) containing details for submitted peptides", type = str, nargs = "+")
    parser.add_argument("--result_dirs", help = "Path to results", type = str, nargs = "+")
    parser.add_argument("outputFile", help = "File to write output to", type = str)
    parser.add_argument("-d", "--debug", action = "store_true", dest = "DEBUG", help = "Flag for setting debug/test state.")
    parser.add_argument("-v", "--verbose", action = "store_true", dest = "VERB", help = "Flag for setting verbose output.")
    args = parser.parse_args()

    ## Set Global Vars
    DEBUG = args.DEBUG
    VERB = args.VERB

    out = open(args.outputFile, "w")
    out.write("barcode\tENST\tmutation\thla\tmut_peptide\tpeptide_variant_position\tmut_ic50\n")

    for peptideInfoFile, result_dir in zip(args.peptideInfoFiles, args.result_dirs):
        ## reset these for each variable
        resI = {}
        pepInfo = {}

        ## Parse peptide info file
        if VERB: print("Parsing peptide info file...")
        parsePeptideInfo(peptideInfoFile)

        ## Parse NetMHCpan output and write
        out = open(args.outputFile, "a")


        if VERB: print("Parsing NetMHCpan results...")
        for root, dirs, filenames in os.walk(result_dir):
            for f in filenames:
                if f.endswith(".pMHC"):
                    if VERB: print("I: Parsing {}".format(os.path.join(root,f)))
                    parseClassI(os.path.join(root,f))

                    ## Write results
                    if VERB: print("Writing output...")
                    writeOutput(out)

                    ## clear resholder
                    resI = {}

        out.close()
    
    print("done.")
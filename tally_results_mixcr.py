'''
Tally results
Go through folders of output from MiTCR run and summarize.

Date: September 22, 2016
@author: sbrown

Edited April 6, 2017:
    - Tally MiXCR data

Edited April 10, 2017:
    - Assigning chain based on J gene, not V gene, as some V alpha genes are also V delta genes...
'''

## Import Libraries
import sys
import argparse
import os

DEBUG = False
VERB = False


if __name__ == "__main__":

    ## Deal with command line arguments
    parser = argparse.ArgumentParser(description = "Tally results")
    ## add_argument("name", "(names)", metavar="exampleOfValue - best for optional", type=int, nargs="+", choices=[allowed,values], dest="nameOfVariableInArgsToSaveAs")
    parser.add_argument("resultsDir", help = "Directory containing result folders", type = str)
    parser.add_argument("outputDir", help = "Directory to save output", type = str)
    parser.add_argument("-d", "--debug", action = "store_true", dest = "DEBUG", help = "Flag for setting debug/test state.")
    parser.add_argument("-v", "--verbose", action = "store_true", dest = "VERB", help = "Flag for setting verbose output.")
    args = parser.parse_args()

    ## Set Global Vars
    DEBUG = args.DEBUG
    VERB = args.VERB

    cdr3 = open(os.path.join(args.outputDir, "mixcrCDR3.tsv"), "w")
    samp = open(os.path.join(args.outputDir, "mixcrSamp.tsv"), "w")

    cdr3.write("sample\tchain\taaSeq\tnucSeq\tabundance\tvallele\tdallele\tjallele\n")
    samp.write("sample\treadLen\tnumReads\n")

    for root, dirs, files in os.walk(args.resultsDir):
        for d in dirs:
            readLen = open(os.path.join(root, d, "readLength.txt"), "r").readlines()[0].rstrip()
            numreads = open(os.path.join(root, d, "numReads.txt"), "r").readlines()[0].rstrip()

            samp.write("{}\t{}\t{}\n".format(d, readLen, numreads))

            ## expected result file name:
            resFile = "alignments_rescued_2_extended.clns.TCR.txt"
        
            HEADER = True

            for line in open(os.path.join(root, d, resFile), "r"):
                ## line 1:
                #[0]cloneId [1]cloneCount      [2]cloneFraction   [3]clonalSequence  [4]clonalSequenceQuality   [5]allVHitsWithScore       [6]allDHitsWithScore       [7]allJHitsWithScore [8]allCHitsWithScore       [9]allVAlignments  [10]allDAlignments  [11]allJAlignments  [12]allCAlignments  [13]nSeqFR1 [14]minQualFR1      [15]nSeqCDR1        [16]minQualCDR1     [17]nSeqFR2 [18]minQualFR2      [19]nSeqCDR2        [20]minQualCDR2     [21]nSeqFR3 [22]minQualFR3      [23]nSeqCDR3        [24]minQualCDR3     [25]nSeqFR4 [26]minQualFR4      [27]aaSeqFR1        [28]aaSeqCDR1       [29]aaSeqFR2        [30]aaSeqCDR2       [31]aaSeqFR3        [32]aaSeqCDR3       [33]aaSeqFR4        [34]refPoints
                if HEADER:
                    HEADER = False
                else:
                    line = line.rstrip().split("\t")
                    chain = line[7][:3]
                    abund = line[1]
                    nucSeq = line[23]
                    aaSeq = line[32]
                    vallele = line[5]
                    #vseg = line[7]
                    jallele = line[7]
                    #jseg = line[9]
                    dallele = line[6]
                    #dseg = line[11]

                    #cdr3.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(d, chain, aaSeq, nucSeq, abund, vallele, dallele, jallele))
                    cdr3.write("\t".join([d, chain, aaSeq, nucSeq, abund, vallele, dallele, jallele]))
                    cdr3.write("\n")
    cdr3.close()
    samp.close()

    print("done.")

'''
Prepare MiXCR run
Given list of files, prepare MiXCR run.

Date: May 29, 2017
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
    parser = argparse.ArgumentParser(description = "Prepare MiXCR run")
    ## add_argument("name", "(names)", metavar="exampleOfValue - best for optional", type=int, nargs="+", choices=[allowed,values], dest="nameOfVariableInArgsToSaveAs")
    parser.add_argument("sourceFof", help = "File with files to analyze", type = str)
    parser.add_argument("outputBasename", help = "File basename for output files (to be used by clusterTAS)", type = str)
    parser.add_argument("-d", "--debug", action = "store_true", dest = "DEBUG", help = "Flag for setting debug/test state.")
    parser.add_argument("-v", "--verbose", action = "store_true", dest = "VERB", help = "Flag for setting verbose output.")
    args = parser.parse_args()

    ## Set Global Vars
    DEBUG = args.DEBUG
    VERB = args.VERB

    fof = open("{}.fof".format(args.outputBasename), "w")
    script = open("{}.scripts".format(args.outputBasename), "w")

    for line in open(args.sourceFof, "r"):
        line = line.rstrip().split("\t")
        runid = line[0]
        rnas = line[1:]

        if len(rnas) == 1 and rnas[0].endswith(".bam"):
            isBam = True
        else:
            isBam = False
            if len(rnas) > 2:
                print("WARNING: More than two sequence files given for {}.".format(runid))
                print("!! Will only use the first two files !!")
            elif len(rnas) == 1:
                print("WARNING: Expects two fastq files per sample, only found 1 ({}). Unexpected results may occur.".format(runid))

        if isBam:
            fof.write("{}\n".format(rnas[0]))

            # MiXCR v2 params
            script.write("{}\t/gsc/software/linux-x86_64-centos5/samtools-1.3.1/bin/samtools view {} | head -1 | awk '{{print $10}}' | wc -c > readLength.txt; /gsc/software/linux-x86_64-centos5/samtools-1.3.1/bin/samtools view {} | wc -l > numReads.txt; /gsc/software/linux-x86_64-centos5/samtools-1.3.1/bin/samtools bam2fq -1 reads_1.fastq -2 reads_2.fastq {}; /home/sbrown/bin/mixcr-2.1/mixcr align -p rna-seq -OallowPartialAlignments=true reads_1.fastq reads_2.fastq alignments.vdjca; /home/sbrown/bin/mixcr-2.1/mixcr assemblePartial -p alignments.vdjca alignments_rescued_1.vdjca; /home/sbrown/bin/mixcr-2.1/mixcr assemblePartial -p alignments_rescued_1.vdjca alignments_rescued_2.vdjca; /home/sbrown/bin/mixcr-2.1/mixcr extendAlignments alignments_rescued_2.vdjca alignments_rescued_2_extended.vdjca; /home/sbrown/bin/mixcr-2.1/mixcr assemble -ObadQualityThreshold=0 alignments_rescued_2_extended.vdjca alignments_rescued_2_extended.clns; /home/sbrown/bin/mixcr-2.1/mixcr exportClones -c TCR alignments_rescued_2_extended.clns alignments_rescued_2_extended.clns.TCR.txt; rm reads*.fastq;\n".format(runid, os.path.basename(bam), os.path.basename(bam), os.path.basename(bam)))
        else:
            fof.write("{}\n".format("\t".join(rnas)))

            #  for numReads, dividing by two instead of 4 because there are two (paired end) files.
            script.write("{}\tzcat {} | head -2 | tail -1 | wc -c > readLength.txt; zcat {} | wc -l | awk '{{print $1/2}}' > numReads.txt; /home/sbrown/bin/mixcr-2.1/mixcr align -p rna-seq -OallowPartialAlignments=true {} {} alignments.vdjca; /home/sbrown/bin/mixcr-2.1/mixcr assemblePartial -p alignments.vdjca alignments_rescued_1.vdjca; /home/sbrown/bin/mixcr-2.1/mixcr assemblePartial -p alignments_rescued_1.vdjca alignments_rescued_2.vdjca; /home/sbrown/bin/mixcr-2.1/mixcr extendAlignments alignments_rescued_2.vdjca alignments_rescued_2_extended.vdjca; /home/sbrown/bin/mixcr-2.1/mixcr assemble -ObadQualityThreshold=0 alignments_rescued_2_extended.vdjca alignments_rescued_2_extended.clns; /home/sbrown/bin/mixcr-2.1/mixcr exportClones -c TCR alignments_rescued_2_extended.clns alignments_rescued_2_extended.clns.TCR.txt;\n".format(runid, os.path.basename(rnas[0]), os.path.basename(rnas[0]), os.path.basename(rnas[0]), os.path.basename(rnas[1])))


    fof.close()
    script.close()
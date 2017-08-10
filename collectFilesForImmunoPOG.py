'''
Collect Files for immunoPOG
Given a file of POG IDs, will collect all files required for immunoPOG analysis

Date: April 21, 2017
@author: sbrown

Edited June 5, 2017:
    - Also accept list of RNA-seq fastq files.
'''

## Import Libraries
import sys
import argparse
import os
import pprint
import glob

DEBUG = False
VERB = False

## Expected file locations (glob format):
pogDir = "/projects/POG/POG_data/"
flatFilePathTemplate = "{pog_id}/flatfile/{pog_id}.tab"
rnaBamPathTemplate = "{pog_id}/rna/{biop_num_t}_t_{lib_id_t}/reposition/hg19a_jg-e69_bwa-*/"
rnaBamEnd = "*_withJunctionsOnGenome_dupsFlagged.bam"
snvFilePathTemplate = "{pog_id}/wgs/{biop_num_t}_t_{lib_id_t}_{biop_num_n}_n_{lib_id_n}/reviewed/snv/{pog_id}[._]{lib_id_t}_{lib_id_n}.passed.somatic.coding.snvs.merged.annotated.txt"
indelFilePathTemplate = "{pog_id}/wgs/{biop_num_t}_t_{lib_id_t}_{biop_num_n}_n_{lib_id_n}/reviewed/indel/{pog_id}_{lib_id_t}_{lib_id_n}.passed.somatic.coding.indels.merged_annotated_collapsed.txt"


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def askChoice(choices):
    '''Ask user which of the list they want'''
    print("Multiple files, which one should be used? Enter '-1' for none.")
    for i in range(0,len(choices)):
        print("[{}]\t{}".format(i,choices[i]))
    response = "null"
    while not is_number(response) or response > len(choices):
        response = input("Enter number: ")
        response = int(response)

    return response


if __name__ == "__main__":

    ## Deal with command line arguments
    parser = argparse.ArgumentParser(description = "Collect Files for immunoPOG")
    parser.add_argument("samplesFile", help = "File listing each POG ID to fetch", type = str)
    parser.add_argument("outputBase", help = "Output file base name", type = str)
    parser.add_argument("-fq", "--fastq_fof", dest = "fqfof", help = "Optional rnaseq fastq fof for requested POG ids.")
    parser.add_argument("-d", "--debug", action = "store_true", dest = "DEBUG", help = "Flag for setting debug/test state.")
    parser.add_argument("-v", "--verbose", action = "store_true", dest = "VERB", help = "Flag for setting verbose output.")
    args = parser.parse_args()

    ## Set Global Vars
    DEBUG = args.DEBUG
    VERB = args.VERB

    lib_info = {}   ## hold all libraries for all POG of interest

    filesToUse = {}

    fqs = []

    if args.fqfof:
        ## if fastq fof provided.
        for line in open(args.fqfof, "r"):
            fqs.append(line.rstrip())
    if DEBUG: print(fqs)


    for POG in open(args.samplesFile, "r"):
        POG = POG.rstrip()
        if VERB: print("Collecting files for {}".format(POG))

        lib_info[POG] = {}
        lib_info[POG]["normal"] = []
        lib_info[POG]["diseased"] = {}

        HEADER = True
        for line in open(os.path.join(pogDir, flatFilePathTemplate.format(pog_id=POG)), "r"):
            if HEADER:
                HEADER = False
            else:
                line = line.rstrip().split("\t")
                lib = line[3]   ## library id
                dn = line[9]    ## diseased/normal
                prot = line[16] ## protocol (wgs/rna)
                samp = line[42] ## sample prefix (biop1, blood1, etc)

                if dn == "Normal" and "WGS" in prot:
                    lib_info[POG]["normal"].append((lib,samp))
                elif dn == "Diseased" and "RNA" in prot:
                    if samp not in lib_info[POG]["diseased"]:
                        lib_info[POG]["diseased"][samp] = {"rna":[], "wgs":[]}
                    lib_info[POG]["diseased"][samp]["rna"].append(lib)
                elif dn == "Diseased" and "WGS" in prot:
                    if samp not in lib_info[POG]["diseased"]:
                        lib_info[POG]["diseased"][samp] = {"rna":[], "wgs":[]}
                    lib_info[POG]["diseased"][samp]["wgs"].append(lib)



        if len(lib_info[POG]["normal"]) > 1:
            print("Multiple normals identified for {}.".format(POG))
            print("Check flat file: {}".format(os.path.join(pogDir, flatFilePathTemplate.format(pog_id=POG))))
        elif len(lib_info[POG]["normal"]) == 0:
            print("No normal found for {}.".format(POG))

        for biop_d in lib_info[POG]["diseased"]:
            for lib_d in lib_info[POG]["diseased"][biop_d]["rna"]:
                ## find RNAseq
                files = []
                if args.fqfof:
                    ## lookup in fq list
                    for f in fqs:
                        if lib_d in f:
                            files.append(f)
                    
                    if "{}-{}".format(POG,biop_d) not in filesToUse:
                        filesToUse["{}-{}".format(POG,biop_d)] = {"rna":"","snv":"","indel":""}

                    filesToUse["{}-{}".format(POG,biop_d)]["rna"] = files

                else:
                    ## lookup bam based on expected path
                    files = glob.glob(os.path.join(pogDir,rnaBamPathTemplate.format(pog_id=POG, biop_num_t=biop_d, lib_id_t=lib_d), rnaBamEnd))
                    if len(files) == 0:
                        print("Missing RNA-seq file for {}, {}, {}.".format(POG, biop_d, lib_d))
                        if glob.glob(os.path.join(pogDir,rnaBamPathTemplate.format(pog_id=POG, biop_num_t=biop_d, lib_id_t=lib_d))):
                            print("Choose file to use from directory...")
                            files = glob.glob(os.path.dirname(os.path.join(pogDir,rnaBamPathTemplate.format(pog_id=POG, biop_num_t=biop_d, lib_id_t=lib_d),"*")))
                            ind = askChoice(files)
                        else:
                            ind = -1

                    if len(files) > 1:
                        ind = askChoice(files)
                    else:
                        ind = 0

                    if ind >= 0:
                        file = files[ind]
                        if "{}-{}".format(POG,biop_d) not in filesToUse:
                            filesToUse["{}-{}".format(POG,biop_d)] = {"rna":"","snv":"","indel":""}

                        filesToUse["{}-{}".format(POG,biop_d)]["rna"] = [file]

            ## find SNV & INDEL
            if "{}-{}".format(POG,biop_d) in filesToUse: ## we have RNA-seq
                for lib_d in lib_info[POG]["diseased"][biop_d]["wgs"]:
                    for lib_n, biop_n in lib_info[POG]["normal"]:
                        snvFiles = glob.glob(os.path.join(pogDir, snvFilePathTemplate.format(pog_id=POG, biop_num_t=biop_d, lib_id_t=lib_d, biop_num_n=biop_n, lib_id_n=lib_n)))
                        if len(snvFiles) == 0:
                            print("No SNV file with expected name.")
                            if os.path.exists(os.path.dirname(os.path.join(pogDir, snvFilePathTemplate.format(pog_id=POG, biop_num_t=biop_d, lib_id_t=lib_d, biop_num_n=biop_n, lib_id_n=lib_n)))):
                                print("Choose file to use from directory...")
                                snvFiles = glob.glob(os.path.join(os.path.dirname(os.path.join(pogDir, snvFilePathTemplate.format(pog_id=POG, biop_num_t=biop_d, lib_id_t=lib_d, biop_num_n=biop_n, lib_id_n=lib_n))),"*"))
                                ind = askChoice(snvFiles)
                            else:
                                ind = -1

                        elif len(snvFiles) > 1:
                            ind = askChoice(snvFiles)
                        else:
                            ind = 0

                        if ind >= 0:
                            snvFile = snvFiles[ind]
                        else:
                            snvFile = None

                        filesToUse["{}-{}".format(POG,biop_d)]["snv"] = snvFile

                        indFiles = glob.glob(os.path.join(pogDir, indelFilePathTemplate.format(pog_id=POG, biop_num_t=biop_d, lib_id_t=lib_d, biop_num_n=biop_n, lib_id_n=lib_n)))
                        if len(indFiles) == 0:
                            print("No INDEL file with expected name.")
                            if os.path.exists(os.path.dirname(os.path.join(pogDir, indelFilePathTemplate.format(pog_id=POG, biop_num_t=biop_d, lib_id_t=lib_d, biop_num_n=biop_n, lib_id_n=lib_n)))):
                                print("Choose file to use from directory...")
                                indFiles = glob.glob(os.path.join(os.path.dirname(os.path.join(pogDir, indelFilePathTemplate.format(pog_id=POG, biop_num_t=biop_d, lib_id_t=lib_d, biop_num_n=biop_n, lib_id_n=lib_n))),"*"))
                                ind = askChoice(indFiles)
                            else:
                                print("Expected directory ({}) does not exist.".format(os.path.dirname(os.path.join(pogDir, indelFilePathTemplate.format(pog_id=POG, biop_num_t=biop_d, lib_id_t=lib_d, biop_num_n=biop_n, lib_id_n=lib_n)))))
                                ind = -1

                        elif len(indFiles) > 1:
                            ind = askChoice(indFiles)
                        else:
                            ind = 0

                        if ind >= 0:
                            indFile = indFiles[ind]
                        else:
                            indFile = None

                        filesToUse["{}-{}".format(POG,biop_d)]["indel"] = indFile


    
    if VERB: pprint.pprint(lib_info)
    if VERB: print("\n\n")
    if VERB: pprint.pprint(filesToUse)

    ## MiXCR requires fastq.
    ## OptiType requires fastq.
    ## NetMHCpan requires the HLA results and the SNV/indels info

    ## think about how I want to set this up. Write lists of files for each probably, with some sort of id.

    
    rnaOut = open("{}_rnaseq.tsv".format(args.outputBase), "w")
    mutOut = open("{}_mutations.tsv".format(args.outputBase), "w")

    for pid in filesToUse:
        if filesToUse[pid]["rna"] is not None and filesToUse[pid]["snv"] is not None and filesToUse[pid]["indel"] is not None:

            rnaOut.write("{}\t{}\n".format(pid, "\t".join(filesToUse[pid]["rna"])))
            mutOut.write("{}\t{}\t{}\n".format(pid, filesToUse[pid]["snv"], filesToUse[pid]["indel"]))

    rnaOut.close()
    mutOut.close()

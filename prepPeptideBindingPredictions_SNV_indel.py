'''
Prepare peptide binding predictions

Input:
    - fof: tab separated file with a sample per line. Fields are sample name, comma-separated HLA alleles (in NetMHCpan format), snv file, indel file.
    - reference transcriptome: which matches genome used for mutation calling
    - reference proteome: which matches genome used for mutation calling
Output:
    - directories for each sample containing the peptides to run predictions on.
    - master list of all "minimal" peptides that contain mutations
    - scripts to run to get peptide-MHC binding predictions for each sample (designed to be run on Genesis)

Date: August 8, 2017
@author: sbrown
'''

## Import Libraries
import sys
import argparse
import os
from Bio.Seq import Seq
from Bio.Alphabet import generic_rna
import math

NETMHCPAN_BIN = "/home/sbrown/bin/netMHCpan-3.0/netMHCpan"  ## on genesis

DEBUG = False
VERB = False

SEQ_DICT = {}
HLA_DICT = {}
enst2ensp = {}

MIN_PEPTIDE_LEN = 8
MAX_PEPTIDE_LEN = 11

numMutationsSkipped = 0
numMutationsComplete = 0


## since ensp reference sequences may contain *
def processProtRefWithStop(protseq):
    protArray = protseq.split("*")
    curBestLen = 0
    bestProt = ""
    for p in protArray:
        if p.startswith("M") and len(p) > curBestLen:
            bestProt = p
            curBestLen = len(p)

    return bestProt


## since enst reference sequences may contain 5'UTR
def getORFandStart(sequence, protseq, geneName):
    '''Given a CDNA sequence and prot seq, find and return the position of the prot (nucleotide position of the A in the start ATG codon) and position of the stop codon'''

    startPos_best = 0
    protein_best = ""
    
    if "U" in protseq:
        protseq = protseq.split("U")
    
        ## For sequence[0:], sequence[1:], sequence[2:] (three different reading frames)
        for i in range(0,3):
            seq = sequence[i:]
            curPos = i + 1

            while len(seq) % 3 != 0:
                seq += "N"

            ## Translate including stops.
            prot = "".join(Seq(seq, generic_rna).translate(to_stop=False))

            ## if adding "N" to the end resulted in an X being the last amino acid, remove it.
            while prot.endswith("X"):
                prot = prot[:-1]

            ## Parse sequence into strings delimited by "*" (stop codon). Keep track of the position of the first codon of each, and the length.
            protArray = prot.split("*")
            for j in range(0,len(protArray)):
                p = protArray[j]
                if p.endswith(protseq[0]):
                    allChunksMatch = True
                    for k in range(1,len(protseq)):
                        if len(protseq[k]) > 0 and protseq[k] == protArray[j+k]:
                            allChunksMatch = allChunksMatch and True
                        elif len(protseq[k]) > 0:
                            ## then there is a mismatch
                            allChunksMatch = False

                    if allChunksMatch:
                        protein_best = p
                        startPos_best = curPos
                
                ## shift cursor to beginning of next protein that will be checked
                curPos += 3*len(p) + 3 #plus 3 for the stop codon

        ## Return the position of the longest.
        if not protein_best.endswith(protseq[0]):
            if DEBUG: print("Reference transcript is: {}".format(sequence))
            if DEBUG: print("Reference protein sequence is: {}".format(protseq))
            print("Warning: No ORF ends with the provided reference protein sequence for {}. SKIPPING.".format(geneName))
            return [None, None]

        else:
            ## if protSeq shorter than protein_best, shift the start position up to the start of the protein.
            startPos_best += 3*(len(protein_best) - len(protseq[0]))
            endPos_best = startPos_best + 3*len("U".join(protseq))

    else:
        ## For sequence[0:], sequence[1:], sequence[2:] (three different reading frames)
        for i in range(0,3):
            seq = sequence[i:]
            curPos = i + 1

            while len(seq) % 3 != 0:
                seq += "N"

            ## Translate including stops.
            prot = "".join(Seq(seq, generic_rna).translate(to_stop=False))

            ## if adding "N" to the end resulted in an X being the last amino acid, remove it.
            while prot.endswith("X"):
                prot = prot[:-1]

            #if DEBUG: print("Potential protein is: {}".format(prot))

            ## Parse sequence into strings delimited by "*" (stop codon). Keep track of the position of the first codon of each, and the length.
            protArray = prot.split("*")
            for p in protArray:
                #if DEBUG: print("startPos_best = {}".format(startPos_best))
                #if len(p) > len(protein_best):
                if p.endswith(protseq):
                    protein_best = p
                    startPos_best = curPos
                    #if DEBUG: print("protein_best = {}".format(p))
                    #if DEBUG: print("refprot = {}".format(protseq))
                    #if DEBUG: print("protein_best start = {}".format(startPos_best))
                
                ## shift cursor to beginning of next protein that will be checked
                curPos += 3*len(p) + 3 #plus 3 for the stop codon

        ## Return the position of the longest.
        if not protein_best.endswith(protseq):
            if DEBUG: print("Reference transcript is: {}".format(sequence))
            if DEBUG: print("Reference protein sequence is: {}".format(protseq))
            print("Warning: No ORF ends with the provided reference protein sequence for {}. SKIPPING.".format(geneName))
            return [None, None]

        else:
            ## if protSeq shorter than protein_best, shift the start position up to the start of the protein.
            startPos_best += 3*(len(protein_best) - len(protseq))
            endPos_best = startPos_best + 3*len(protseq)

    #if DEBUG: print("startPos_best = {}\nendPos_best = {}".format(startPos_best, endPos_best))
    return [startPos_best, endPos_best]



if __name__ == "__main__":

    ## Deal with command line arguments
    parser = argparse.ArgumentParser(description = "Prepare peptide binding predictions")
    parser.add_argument("infiles", help = "FOF with sample name, hla alleles, snv file, and indel file (tab sep)", type = str)
    parser.add_argument("ref_trans", help = "Reference transcriptome (matches genome version for mutation calls)", type = str)
    parser.add_argument("ref_prot", help = "Reference proteome (matches genome version for mutation calls)", type = str)
    parser.add_argument("workingDir", help = "Path to directory to write files.", type = str)
    parser.add_argument("-d", "--debug", action = "store_true", dest = "DEBUG", help = "Flag for setting debug/test state.")
    parser.add_argument("-v", "--verbose", action = "store_true", dest = "VERB", help = "Flag for setting verbose output.")
    args = parser.parse_args()

    ## Set Global Vars
    DEBUG = args.DEBUG
    VERB = args.VERB

    print("Reading in sample file...")
    ## sample name \t comma-separated HLA alleles \t snv file \t indel file \n
    files = []
    for line in open(args.infiles, "r"):
        line = line.rstrip().split("\t")
        files.append([line[0], line[1], line[2], line[3]])

    ## Read in reference transcriptome
    ## {ENST: [cDNA_SEQUENCE,prot_SEQUENCE]}
    print("Reading in reference transcriptome...")
    for line in open(args.ref_trans, "r"):
        if line.startswith(">"):
            header = line.rstrip()
            trans = header.split(" ")[0]
            ## cleave off ">"
            trans = trans[1:]
            ## remove ".#" at end
            trans = trans.split(".")[0]
            if trans not in SEQ_DICT:
                SEQ_DICT[trans] = ["",""]
            else:
                sys.exit("Error: Transcript {} already has sequence:\n\n{}\n\nExiting.".format(trans, SEQ_DICT[trans]))
        else:
            SEQ_DICT[trans][0] += line.rstrip()

    ## Read in reference proteome
    print("Reading in reference proteome...")
    for line in open(args.ref_prot, "r"):
        if line.startswith(">"):
            header = line.rstrip()
            trans = header.split(" ")[4].split(":")[1]
            ## remove ".#" at end
            trans = trans.split(".")[0]
            
        elif trans in SEQ_DICT: ## we have cDNA seq
            SEQ_DICT[trans][1] += line.rstrip()

    


    ## snv header (may not be exact):
    ## [0]chromosome [1]pos [2]ref [3]alt [4]dbsnp [5]cosmic_id [6]effect_type [7]functional_class [8]impact [9]hgvs_protein [10]aa_change [11]hgvs_cds [12]gene_symbol [13]ensembl_gene_id [14]ensembl_transcript_id [15]mutationseq_score [16]Strelka [17]P01598_Alt_Count [18]P01598_Ref_Count [19]P01606_Alt_Count [20]P01606_Ref_Count [21]P01608_Alt_Count [22]P01608_Ref_Count [23]is_real [24]comments [25]is_best_transcript [26]picked_by [27]other_transcripts [28]other_genes [29]cDNA_pos [30]Sift [31]Polyphen [32]Existing_variations [33]Domain_(type,desc,interpro_id) [34]type

    ## indel header (may not be exact):
    ## [0]chromosome [1]pos [2]ref [3]alt [4]type [5]dbsnp [6]cosmic_id [7]effect_type [8]functional_class [9]impact [10]hgvs_protein [11]aa_change [12]hgvs_cds [13]gene_symbol [14]ensembl_gene_id [15]ensembl_transcript_id [16]key [17]gv_effect [18]evidence [19]is_real [20]comments [21]P01598_Alt_Count [22]P01598_Total_Count [23]P01606_Alt_Count [24]P01606_Total_Count [25]P01608_Alt_Count [26]P01608_Total_Count [27]variant_source [28]called_by [29]is_best_transcript [30]picked_by [31]other_transcripts [32]other_genes

    print("Reading in mutation information...")

    peps = {}
    pepsInfo = {}
    PEP_NUM = 0

    for tumour, hlas, snvFile, indelFile in files:
        print("Processing {}...".format(tumour))
        HLA_DICT[tumour] = hlas

        ## SNV
        for line in open(snvFile, "r"):
            if not line.startswith("#"):

                VALID = True

                line = line.rstrip().split("\t")
                
                ## If missense variant
                ## Note there is also "stop_lost" which might be worth predicting, as an OOF-indel. (new sequence generated)
                ## But this is not currently implemented.
                if "missense_variant" in line[6]:
                    #if DEBUG: print(line)
                    hgvs = line[11]
                    aaChange = line[10]
                    snvPos = int(hgvs[2:-3])
                    ref = hgvs[-3:].split(">")[0]
                    mut = hgvs[-3:].split(">")[1]
                    hugo = line[12]
                    enst = line[14]


                    ## get reference sequence
                    enst_seq = SEQ_DICT[enst][0]
                    ensp_seq = SEQ_DICT[enst][1]

                    ## cleave off leading X if it is present.
                    while ensp_seq.startswith("X"):
                        ensp_seq = ensp_seq[1:]

                    if not ensp_seq.startswith("M"):
                        print("Warning: Reference protein sequence for {} ({}) does not start with 'M'. SKIPPING".format(enst, hugo))
                        numMutationsSkipped += 1
                        VALID = False

                    ## if "*" in ensp_seq, extract the longest sequence
                    if VALID and "*" in ensp_seq:
                        ensp_seq = processProtRefWithStop(ensp_seq)

                        if len(ensp_seq) == 0:
                            ## no appropriate sequence
                            print("Warning: Reference protein sequence for {} ({}) contains '*', and no segment starts with 'M'. SKIPPING".format(enst, hugo))
                            numMutationsSkipped += 1
                            VALID = False

                    if VALID:
                        #if DEBUG: print("Getting start codon position for {}, indel {}".format(enst, hgvs))
                        if DEBUG: print("Getting protein position in {}.".format(enst))
                        orfPos, stopPos = getORFandStart(enst_seq, ensp_seq, "{} ({})".format(enst, hugo))

                    if VALID and orfPos != None:

                        #if DEBUG: print("Coding start position is {}: {}...".format(orfPos, enst_seq[orfPos-1:orfPos+5]))


                        ## check the reference matches what mutation caller said was reference
                        refSeq = enst_seq[orfPos-1:]

                        noStartStart = False
                            
                        if not refSeq.startswith("ATG"):
                            print('Warning: Reference sequence for {} ({}) starting at position {} (based on reference protein matching ORF) does not start with a start codon.'.format(enst, hugo, orfPos))
                            noStartStart = True
                            numMutationsSkipped += 1
                            print("SKIPPING")

                        ## Move this outside the else if non-Start codon starts are fine to include.
                        else:
                        
                            transRefNuc = refSeq[snvPos-1:snvPos]
                            if ref != transRefNuc:
                                print("Warning: Reference allele does not match reference sequence: {}: {}/{}".format(enst, ref, transRefNuc))
                                print(line)
                                numMutationsSkipped += 1
                                print("SKIPPING")
                            else:

                                ## Get protein codon/aa
                                protMutPos = int(aaChange[1:-1])
                                
                                ## Translate reference
                                ## if length is not multiple of 3, pad the end with N
                                while len(refSeq)%3 != 0:
                                    refSeq += "N"

                                ## refSeq is from start codon onwards.
                                ## position of stop codon is known. Any stop codons before that point need to be replaced with "U"
                                refSeqTrans = "".join(Seq(refSeq, generic_rna).translate(to_stop=False)).split("*")

                                i = 0
                                refProtein = refSeqTrans[i]
                                ## while this does not reach the stop codon
                                while 3*len(refProtein) < stopPos - orfPos:
                                    ## add subsequent chunks delimited by "U"
                                    i += 1
                                    refProtein += "U"
                                    refProtein += refSeqTrans[i]


                                ## mutate reference cdna
                                mutCDNA = refSeq[:snvPos-1] + mut + refSeq[snvPos:]
                                while len(mutCDNA)%3 != 0:
                                    mutCDNA += "N"

                                if mutCDNA.startswith("ATG"):
                                    ## start codon is not broken by mutation

                                    ## translate mutant
                                    mutCDNAtrans = "".join(Seq(mutCDNA, generic_rna).translate(to_stop=False)).split("*")

                                    i = 0
                                    mutProtein = mutCDNAtrans[i]
                                    while 3*len(mutProtein) < stopPos - orfPos:
                                        i += 1
                                        mutProtein += "U"
                                        mutProtein += mutCDNAtrans[i]


                                    if DEBUG and noStartStart: print("Mutation line: {}".format(line))
                                    if DEBUG and noStartStart: print("Mutation: {}".format(aaChange))
                                    if DEBUG and noStartStart: print("Reference protein: {}.....{}.....{}".format(refProtein[:protMutPos-1], refProtein[protMutPos-1:protMutPos], refProtein[protMutPos:]))
                                    if DEBUG and noStartStart: print("Mutant protein: {}.....{}.....{}".format(mutProtein[:protMutPos-1], mutProtein[protMutPos-1:protMutPos], mutProtein[protMutPos:]))



                                    ## get minimal peptides
                                    ## bounds calculated here are 0-indexed.
                                    if protMutPos - MAX_PEPTIDE_LEN <= 0:
                                        ## the MAX_PEPTIDE_LEN-mer with mutation at the C-terminal end starts before the start codon
                                        lowerBound = 0
                                    else:
                                        lowerBound = protMutPos - MAX_PEPTIDE_LEN

                                    if protMutPos + MAX_PEPTIDE_LEN - 1 > len(mutProtein):
                                        ## the MAX_PEPTIDE_LEN-mer with the mutation at the N-terminal end ends after the stop codon
                                        upperBound = len(mutProtein) - 1
                                    else:
                                        upperBound = protMutPos + MAX_PEPTIDE_LEN - 2

                                    ## extract the peptide
                                    refPeptide = refProtein[lowerBound:upperBound+1]
                                    mutPeptide = mutProtein[lowerBound:upperBound+1]
                                    ## variant position within the peptide:
                                    varPos = protMutPos - lowerBound

                                    if len(refPeptide) < MIN_PEPTIDE_LEN or len(mutPeptide) < MIN_PEPTIDE_LEN:
                                        print("Warning: Variant {} does not fall within the coding sequence of {}. SKIPPING".format(hgvs, enst))
                                        numMutationsSkipped += 1
                                        print("refProtein: {}\nrefProtein length: {}\nmutProtein: {}".format(refProtein, len(refProtein), mutProtein))
                                        print("Expected protein mutation: {}".format(aaChange))
                                    
                                    else:
                                        ## add peptide to list.
                                        if tumour not in peps:
                                            peps[tumour] = {}
                                            pepsInfo[tumour] = {}

                                        mut_type = "SNV"
                                        
                                        ## Reference
                                        ## UNCOMMENT THESE LINES IF YOU WANT TO RUN PREDICTIONS ON REFERENCE SEQUENCE
                                        #PEP_NUM += 1
                                        #peps[tumour][PEP_NUM] = refPeptide
                                        #pepsInfo[tumour][PEP_NUM] = "{}\t{}\t{}\t{}\t{}\t{}\n".format(enst, "{}-{}-{}".format(mut_type,hgvs,aaChange), varPos, "ref", refPeptide, refProtein)

                                        ## Mutant
                                        PEP_NUM += 1
                                        peps[tumour][PEP_NUM] = mutPeptide
                                        pepsInfo[tumour][PEP_NUM] = "{}\t{}\t{}\t{}\t{}\t{}\n".format(enst, "{}-{}-{}".format(mut_type,hgvs,aaChange), varPos, "mut", mutPeptide, refProtein)

                        

        ## Process INDEL file
        for line in open(indelFile, "r"):
            if not line.startswith("#"):

                VALID = True

                line = line.rstrip().split("\t")

                hgvs = line[12]
                #if hgvs != "NA":
                if "NA" not in hgvs and "splice" not in line[7]:
                    mutPeptide = "temp" ## this may get set to None later on if deemed not able to be predicted.
                    ## hgvs_cds = c.#[_#][del/dup/ins]*
                    hgvs = hgvs[2:] ## remove leading "c."

                    if "del" in hgvs:
                        pos, seq = hgvs.split("del")
                    elif "dup" in hgvs:
                        pos, seq = hgvs.split("dup")
                    elif "ins" in hgvs:
                        pos, seq = hgvs.split("ins")
                    else:
                        print("Unknown variant type: {}".format(hgvs))
                        sys.exit()

                    pos = pos.split("_")
                    indelStartPos = int(pos[0])
                    
                    if len(pos) > 1:
                        if "*" in pos[1]:
                            ## affects the stop codon
                            if "del" in hgvs:
                                indelEndPos = indelStartPos + len(seq)
                            else:
                                indelEndPos = indelStartPos + 1
                        else:
                            indelEndPos = int(pos[1])
                    else:
                        indelEndPos = indelStartPos

                    if "dup" in hgvs:
                        ## generalize to "ins"
                        dupStartPos = indelStartPos
                        indelStartPos = indelEndPos
                        indelEndPos = indelStartPos + 1

                    if indelStartPos < 1:
                        print("WARNING: INDEL {} is (at least partially) positioned upstream of start codon. SKIPPING".format(hgvs))
                        numMutationsSkipped += 1
                    else:

                        enst = line[15]
                        # get ORF position
                        enst_seq = SEQ_DICT[enst][0]
                        ensp_seq = SEQ_DICT[enst][1]

                        ## cleave off leading X if it is present.
                        while ensp_seq.startswith("X"):
                            ensp_seq = ensp_seq[1:]

                        if not ensp_seq.startswith("M"):
                            print("Warning: Reference protein sequence for {} ({}) does not start with 'M'. SKIPPING".format(enst, hugo))
                            numMutationsSkipped += 1
                            VALID = False

                        ## if "*" in ensp_seq, extract the longest sequence
                        if VALID and "*" in ensp_seq:
                            ensp_seq = processProtRefWithStop(ensp_seq)

                            if len(ensp_seq) == 0:
                                ## no appropriate sequence
                                print("Warning: Reference protein sequence for {} ({}) contains '*', and no segment starts with 'M'. SKIPPING".format(enst, hugo))
                                numMutationsSkipped += 1
                                VALID = False

                        if VALID:

                            #if DEBUG: print("Getting start codon position for {}, indel {}".format(enst, hgvs))
                            orfPos, stopPos = getORFandStart(enst_seq, ensp_seq, "{} ({})".format(enst, hugo))

                        if VALID and orfPos != None:
                        
                            ## check the reference matches what mutation caller said was reference
                            refSeq = enst_seq[orfPos-1:]
                            while len(refSeq)%3 != 0:
                                refSeq += "N"
                            if "del" in hgvs:
                                if not refSeq[indelStartPos-1:].startswith(seq):
                                    print("Warning: deletion sequence does not match reference! SKIPPING")
                                    numMutationsSkipped += 1
                                    mutPeptide = None
                            if "dup" in hgvs:
                                if not refSeq[dupStartPos-1:].startswith(seq):
                                    print("Warning: duplication base does not match reference! SKIPPING")
                                    numMutationsSkipped += 1
                                    mutPeptide = None

                            ## I think this if/elif is redundant now. Should be handled above.
                            if not ensp_seq.startswith("M"):
                                print("Warning: Reference protein sequence for {} ({}) does not start with 'M'. SKIPPING".format(enst, hugo))
                                numMutationsSkipped += 1
                            elif not refSeq.startswith("ATG"):
                                print('Warning: Reference sequence for {} ({}) starting at position {} (based on longest ORF) does not start with a start codon.'.format(enst, hugo, orfPos))
                                numMutationsSkipped += 1
                                #print("enst_seq: {}".format(enst_seq))
                                #print("ensp_seq: {}".format(ensp_seq))
                                #print("refSeq: {}".format(refSeq))
                                print("SKIPPING")
                            else:

                                ## make the mutation
                                mutSeq = ""

                                stopPosMut = stopPos

                                if "del" in hgvs:
                                    ## remove the deletion sequence
                                    mutSeq = refSeq[:indelStartPos-1] + refSeq[indelEndPos:]

                                    ## shift the stop codon back
                                    stopPosMut -= (indelEndPos - indelStartPos + 1)

                                # Converting dup to ins.
                                elif "ins" in hgvs or "dup" in hgvs:
                                    ## insert the new or duplicated sequence
                                    mutSeq = refSeq[:indelStartPos] + seq + refSeq[indelEndPos-1:]

                                    ## shift the stop codon forward
                                    stopPosMut += len(seq)

                                    
                                if not mutSeq.startswith("ATG"):
                                    print("Warning: Start codon broken by mutation: {}. SKIPPING".format(hgvs))
                                    numMutationsSkipped += 1
                                    mutPeptide = None
                                if len(refSeq)%3 != len(mutSeq)%3 and "frameshift" not in line[7]:
                                    print("Warning: Applying indel caused frameshift, but indel not classified as frameshift: {} - {}. SKIPPING".format(hgvs, line[7]))
                                    numMutationsSkipped += 1
                                    mutPeptide = None
                                if len(refSeq)%3 == len(mutSeq)%3 and "inframe" not in line[7]:
                                    print("Warning: Applying indel did not cause frameshift, but indel classified as frameshift: {} - {}. SKIPPING".format(hgvs, line[7]))
                                    numMutationsSkipped += 1
                                    mutPeptide = None

                                while len(mutSeq) % 3 != 0:
                                    mutSeq += "N"
                                

                                ## manual checking that mut is happening correctly.
                                if DEBUG:
                                    print("INDEL line: {}".format(line))
                                    if "del" in hgvs:
                                        print("Reference: {}".format(refSeq[:indelStartPos-1] + "....." + refSeq[indelStartPos-1:indelEndPos] + "....." + refSeq[indelEndPos:]))
                                        print("Mutant: {}".format(refSeq[:indelStartPos-1] + "....." + refSeq[indelEndPos:]))
                                    elif "dup" in hgvs:
                                        print("Reference: {}".format(refSeq[:indelStartPos] + "....." + refSeq[indelEndPos:]))
                                        print("Mutant: {}".format(refSeq[:indelStartPos] + "....." + seq + "....." + refSeq[indelEndPos:]))
                                    elif "ins" in hgvs:
                                        print("Reference: {}".format(refSeq[:indelStartPos] + "....." + refSeq[indelEndPos-1:]))
                                        print("Mutant: {}".format(refSeq[:indelStartPos] + "....." + seq + "....." + refSeq[indelEndPos-1:]))

                                    else:
                                        print("Unknown mutation type - INVESTIGATE!")

                                    input("Press enter to continue...")



                                ## get minimal peptide

                                ## frameshift_variant, [disruptive_]inframe_deletion, [disruptive_]inframe_insertion
                                if mutPeptide and "inframe" in line[7]:
                                    if "del" in hgvs:
                                        ## variant will be novel aa junction (if codon(s) deleted), or novel aa(s) (if oof triplet(s) deleted)
                                        if indelStartPos%3 == 1:
                                            ## non-disruptive
                                            ## novel aa junction.
                                            ## variantPos is between the aa.
                                            firstIndelCodon = math.floor((indelStartPos+2)/3)
                                            if firstIndelCodon < 1:
                                                print("Warning: INDEL begins 5' of coding sequence. SKIPPING")
                                                numMutationsSkipped += 1
                                                mutPeptide = None
                                            else:
                                                lastIndelCodon = math.floor((indelEndPos+2)/3)
                                                varPos = firstIndelCodon-1 + 0.5   ## between aa just before deletion and new codon after deletion
                                                #firstCodon = firstIndelCodon - 10
                                                firstCodon = math.ceil(varPos) - MAX_PEPTIDE_LEN + 1
                                                if firstCodon < 1:
                                                    firstCodon = 1
                                                
                                                pepVarPos = varPos - firstCodon + 1
                                                
                                                lastCodon = math.floor(varPos) + MAX_PEPTIDE_LEN - 1
                                                



                                                mutProteinArray = "".join(Seq(mutSeq, generic_rna).translate(to_stop=False)).split("*")

                                                i = 0
                                                mutProtein = mutProteinArray[i]
                                                ## while this does not reach the stop codon
                                                while 3*len(mutProtein) < stopPosMut - orfPos:
                                                    ## add subsequent chunks delimited by "U"
                                                    i += 1
                                                    mutProtein += "U"
                                                    mutProtein += mutProteinArray[i]

                                                mutPeptide = mutProtein[firstCodon-1:lastCodon]

                                                ## check that variant site is within peptide
                                                ## (if mutation generates a stop codon at the mutation site, peptide will end before the pepVarPos)
                                                if pepVarPos > len(mutPeptide):
                                                    ## outside, skip this one.
                                                    numMutationsSkipped += 1
                                                    mutPeptide = None

                                                refProteinArray = "".join(Seq(refSeq, generic_rna).translate(to_stop=False)).split("*")

                                                i = 0
                                                refProtein = refProteinArray[i]
                                                ## while this does not reach the stop codon
                                                while 3*len(refProtein) < stopPos - orfPos:
                                                    ## add subsequent chunks delimited by "U"
                                                    i += 1
                                                    refProtein += "U"
                                                    refProtein += refProteinArray[i]
                                                #refPeptide = refProtein[firstCodon-1:lastCodon+int((len(seq)+2)/3)]

                                                aaChange = "del{}_{}".format(firstIndelCodon, lastIndelCodon)

                                        else:
                                            ## disruptive, but result still in frame.
                                            ## new aa and deleted aa(s)
                                            firstIndelCodon = math.floor((indelStartPos+2)/3)
                                            if firstIndelCodon < 1:
                                                print("Warning: INDEL begins 5' of coding sequence. SKIPPING")
                                                numMutationsSkipped += 1
                                                mutPeptide = None
                                            else:
                                                lastIndelCodon = math.floor((indelEndPos+2)/3)
                                                varPos = firstIndelCodon    ## will be the modified/new aa
                                                firstCodon = varPos - MAX_PEPTIDE_LEN + 1
                                                if firstCodon < 1:
                                                    firstCodon = 1

                                                pepVarPos = varPos - firstCodon + 1

                                                lastCodon = varPos + MAX_PEPTIDE_LEN - 1

                                                mutProteinArray = "".join(Seq(mutSeq, generic_rna).translate(to_stop=False)).split("*")

                                                i = 0
                                                mutProtein = mutProteinArray[i]
                                                ## while this does not reach the stop codon
                                                while 3*len(mutProtein) < stopPosMut - orfPos:
                                                    ## add subsequent chunks delimited by "U"
                                                    i += 1
                                                    mutProtein += "U"
                                                    mutProtein += mutProteinArray[i]

                                                mutPeptide = mutProtein[firstCodon-1:lastCodon]

                                                ## check that variant site is within peptide
                                                ## (if mutation generates a stop codon at the mutation site, peptide will end before the pepVarPos)
                                                if pepVarPos > len(mutPeptide):
                                                    ## outside, skip this one.
                                                    numMutationsSkipped += 1
                                                    mutPeptide = None

                                                refProteinArray = "".join(Seq(refSeq, generic_rna).translate(to_stop=False)).split("*")

                                                i = 0
                                                refProtein = refProteinArray[i]
                                                ## while this does not reach the stop codon
                                                while 3*len(refProtein) < stopPos - orfPos:
                                                    ## add subsequent chunks delimited by "U"
                                                    i += 1
                                                    refProtein += "U"
                                                    refProtein += refProteinArray[i]

                                                aaChange = "disdel{}_{}".format(firstIndelCodon, lastIndelCodon)



                                    elif "ins" in hgvs or "dup" in hgvs:
                                        ## minimal peptide will be indelStartPos - 10aa to indelEndPos + 10aa
                                        if indelStartPos%3 == 0:
                                            ## non-disruptive
                                            ## novel aa(s)
                                            ## variantPos is the aa(s).
                                            firstIndelCodon = math.ceil((indelStartPos+2)/3)
                                            if firstIndelCodon < 1:
                                                print("Warning: INDEL begins 5' of coding sequence. SKIPPING")
                                                numMutationsSkipped += 1
                                                mutPeptide = None
                                            else:
                                                lastIndelCodon = firstIndelCodon + int(len(seq)/3) - 1 ## first codon plus number of amino acids added - 1
                                                varPosLeft = firstIndelCodon
                                                varPosRight = lastIndelCodon
                                                
                                                firstCodon = varPosLeft - MAX_PEPTIDE_LEN + 1
                                                if firstCodon < 1:
                                                    firstCodon = 1
                                                
                                                pepVarPosLeft = varPosLeft - firstCodon + 1
                                                pepVarPosRight = varPosRight - firstCodon + 1
                                                
                                                lastCodon = varPosRight + MAX_PEPTIDE_LEN - 1

                                                mutProteinArray = "".join(Seq(mutSeq, generic_rna).translate(to_stop=False)).split("*")

                                                i = 0
                                                mutProtein = mutProteinArray[i]
                                                ## while this does not reach the stop codon
                                                while 3*len(mutProtein) < stopPosMut - orfPos:
                                                    ## add subsequent chunks delimited by "U"
                                                    i += 1
                                                    mutProtein += "U"
                                                    mutProtein += mutProteinArray[i]

                                                mutPeptide = mutProtein[firstCodon-1:lastCodon]

                                                ## check that variant site is within peptide
                                                ## (if mutation generates a stop codon at the mutation site, peptide will end before the pepVarPos)
                                                if pepVarPosLeft > len(mutPeptide):
                                                    ## outside, skip this one.
                                                    numMutationsSkipped += 1
                                                    mutPeptide = None

                                                refProteinArray = "".join(Seq(refSeq, generic_rna).translate(to_stop=False)).split("*")

                                                i = 0
                                                refProtein = refProteinArray[i]
                                                ## while this does not reach the stop codon
                                                while 3*len(refProtein) < stopPos - orfPos:
                                                    ## add subsequent chunks delimited by "U"
                                                    i += 1
                                                    refProtein += "U"
                                                    refProtein += refProteinArray[i]

                                                aaChange = "ins{}_{}".format(firstIndelCodon, lastIndelCodon)
                                                pepVarPos = "{}-{}".format(pepVarPosLeft, pepVarPosRight)

                                        else:
                                            ## disruptive, but result still in frame.
                                            ## new aa(s) and modified aa(s)
                                            firstIndelCodon = math.floor((indelStartPos+2)/3)
                                            if firstIndelCodon < 1:
                                                print("Warning: INDEL begins 5' of coding sequence. SKIPPING")
                                                numMutationsSkipped += 1
                                                mutPeptide = None
                                            else:
                                                #lastIndelCodon = firstIndelCodon + int(len(seq)/3) - 1 ## first codon plus number of amino acids added - 1
                                                lastIndelCodon = math.floor((indelStartPos + len(seq) + 2)/3)
                                                varPosLeft = firstIndelCodon
                                                varPosRight = lastIndelCodon
                                                
                                                firstCodon = math.ceil(varPosLeft) - MAX_PEPTIDE_LEN + 1
                                                if firstCodon < 1:
                                                    firstCodon = 1
                                                
                                                pepVarPosLeft = varPosLeft - firstCodon + 1
                                                pepVarPosRight = varPosRight - firstCodon + 1
                                                
                                                lastCodon = math.floor(varPosRight) + MAX_PEPTIDE_LEN - 1
                                                
                                                mutProteinArray = "".join(Seq(mutSeq, generic_rna).translate(to_stop=False)).split("*")

                                                i = 0
                                                mutProtein = mutProteinArray[i]
                                                ## while this does not reach the stop codon
                                                while 3*len(mutProtein) < stopPosMut - orfPos:
                                                    ## add subsequent chunks delimited by "U"
                                                    i += 1
                                                    mutProtein += "U"
                                                    mutProtein += mutProteinArray[i]

                                                mutPeptide = mutProtein[firstCodon-1:lastCodon]

                                                ## check that variant site is within peptide
                                                ## (if mutation generates a stop codon at the mutation site, peptide will end before the pepVarPos)
                                                if pepVarPosLeft > len(mutPeptide):
                                                    ## outside, skip this one.
                                                    numMutationsSkipped += 1
                                                    mutPeptide = None

                                                refProteinArray = "".join(Seq(refSeq, generic_rna).translate(to_stop=False)).split("*")

                                                i = 0
                                                refProtein = refProteinArray[i]
                                                ## while this does not reach the stop codon
                                                while 3*len(refProtein) < stopPos - orfPos:
                                                    ## add subsequent chunks delimited by "U"
                                                    i += 1
                                                    refProtein += "U"
                                                    refProtein += refProteinArray[i]

                                                aaChange = "disins{}_{}".format(firstIndelCodon, lastIndelCodon)
                                                pepVarPos = "{}-{}".format(pepVarPosLeft, pepVarPosRight)



                                elif mutPeptide and "frameshift" in line[7]:
                                    aaChange = "fs"

                                    if "del" in hgvs:
                                        ## variant will be novel aas
                                        ## doesn't matter where in codon it is affected, 5' end onwards will be affected.

                                        firstIndelCodon = math.floor((indelStartPos+2)/3)
                                        if firstIndelCodon < 1:
                                            print("Warning: INDEL begins 5' of coding sequence. SKIPPING")
                                            numMutationsSkipped += 1
                                            mutPeptide = None
                                        else:
                                            varPos = firstIndelCodon
                                            firstCodon = math.ceil(varPos) - MAX_PEPTIDE_LEN + 1
                                            if firstCodon < 1:
                                                firstCodon = 1
                                            
                                            pepVarPos = varPos - firstCodon + 1
                                            
                                            mutProteinArray = "".join(Seq(mutSeq, generic_rna).translate(to_stop=False)).split("*")

                                            i = 0
                                            mutProtein = mutProteinArray[i]
                                            ## while this does not reach the first codon affected by the indel (any UGA after indel treated as Stop)
                                            while len(mutProtein) < firstIndelCodon-1:
                                                ## add subsequent chunks delimited by "U"
                                                i += 1
                                                mutProtein += "U"
                                                mutProtein += mutProteinArray[i]

                                            mutPeptide = mutProtein[firstCodon-1:]

                                            ## check that variant site is within peptide
                                            ## (if mutation generates a stop codon at the mutation site, peptide will end before the pepVarPos)
                                            if pepVarPos > len(mutPeptide):
                                                ## outside, skip this one.
                                                numMutationsSkipped += 1
                                                mutPeptide = None

                                            pepVarPos = "{}-*".format(pepVarPos)

                                            refProteinArray = "".join(Seq(refSeq, generic_rna).translate(to_stop=False)).split("*")

                                            i = 0
                                            refProtein = refProteinArray[i]
                                            ## while this does not reach the stop codon
                                            while 3*len(refProtein) < stopPos - orfPos:
                                                ## add subsequent chunks delimited by "U"
                                                i += 1
                                                refProtein += "U"
                                                refProtein += refProteinArray[i]


                                    elif "ins" in hgvs or "dup" in hgvs:
                                        ## novel aa(s)
                                        ## variantPos is the aa(s)
                                        firstIndelCodon = math.ceil((indelStartPos+2)/3)
                                        if firstIndelCodon < 1:
                                            print("Warning: INDEL begins 5' of coding sequence. SKIPPING")
                                            numMutationsSkipped += 1
                                            mutPeptide = None
                                        else:
                                            varPosLeft = firstIndelCodon
                                            #varPosRight = "*"
                                            
                                            firstCodon = math.ceil(varPosLeft) - MAX_PEPTIDE_LEN + 1
                                            if firstCodon < 1:
                                                firstCodon = 1
                                            
                                            pepVarPosLeft = varPosLeft - firstCodon + 1
                                            pepVarPosRight = "*"
                                            
                                            mutProteinArray = "".join(Seq(mutSeq, generic_rna).translate(to_stop=False)).split("*")

                                            i = 0
                                            mutProtein = mutProteinArray[i]
                                            ## while this does not reach the first codon affected by the indel (any UGA after indel treated as Stop)
                                            while len(mutProtein) < firstIndelCodon-1:
                                                ## add subsequent chunks delimited by "U"
                                                i += 1
                                                mutProtein += "U"
                                                mutProtein += mutProteinArray[i]
                                                
                                            mutPeptide = mutProtein[firstCodon-1:]

                                            ## check that variant site is within peptide
                                            ## (if mutation generates a stop codon at the mutation site, peptide will end before the pepVarPos)
                                            if pepVarPosLeft > len(mutPeptide):
                                                ## outside, skip this one.
                                                numMutationsSkipped += 1
                                                mutPeptide = None

                                            refProteinArray = "".join(Seq(refSeq, generic_rna).translate(to_stop=False)).split("*")

                                            i = 0
                                            refProtein = refProteinArray[i]
                                            ## while this does not reach the stop codon
                                            while 3*len(refProtein) < stopPos - orfPos:
                                                ## add subsequent chunks delimited by "U"
                                                i += 1
                                                refProtein += "U"
                                                refProtein += refProteinArray[i]

                                            pepVarPos = "{}-{}".format(pepVarPosLeft, pepVarPosRight)

                                    
                                

                                if mutPeptide and mutPeptide.endswith("X"):
                                    ## mutation was near 3' end and required "N" bases to be added to end of transcript, resulting in X amino acid.
                                    ## as this is the end of the transcript, end translation here prior to "X" codon.
                                    mutPeptide = mutPeptide[:-1]

                                ## Manually check the minimal peptide
                                if DEBUG:
                                    print("Reference minimal peptide:\t{}".format(refPeptide))
                                    print("Mutant minimal peptide:\t\t{}".format(mutPeptide))

                                    input("Press enter to continue...")



                                ## add peptide to list.
                                if mutPeptide:
                                    if tumour not in peps:
                                        peps[tumour] = {}
                                        pepsInfo[tumour] = {}

                                    mut_type = "INDEL"
                                    #varPos = "{}-{}".format(varPosLeft,varPosRight)
                                    
                                    ## Reference
                                    ## UNCOMMENT THESE LINES IF YOU WANT TO RUN PREDICTIONS ON REFERENCE SEQUENCE
                                    ## HARDER TO INTERPRET "reference" FOR AN INDEL
                                    #PEP_NUM += 1
                                    #peps[tumour][PEP_NUM] = refPeptide
                                    #pepsInfo[tumour][PEP_NUM] = "{}\t{}\t{}\t{}\t{}\n".format(enst, "{}-{}-{}".format(mut_type,hgvs,aaChange), varPos, "ref", refPeptide)

                                    ## Mutant
                                    PEP_NUM += 1
                                    peps[tumour][PEP_NUM] = mutPeptide
                                    pepsInfo[tumour][PEP_NUM] = "{}\t{}\t{}\t{}\t{}\t{}\n".format(enst, "{}-{}-{}".format(mut_type,hgvs,aaChange), pepVarPos, "mut", mutPeptide, refProtein)





   
    ## make directories and write files
    print("Writing directories and files...")
    transferfof = open(os.path.join(args.workingDir, "filesToTransfer.fof"), "w")
    analysisfile = open(os.path.join(args.workingDir, "analysisToRun.sh"), "w")

    for tum in peps:
        os.makedirs(os.path.join(args.workingDir, tum))
        
        ## get HLA
        hlas = HLA_DICT[tum]

        ## write peptides
        out = open(os.path.join(args.workingDir, tum, "peptides.fa"), "w")
        toWrite = ""
        for k in peps[tum]:
            numMutationsComplete += 1
            toWrite += ">{}\n{}\n".format(k, peps[tum][k])
        out.write(toWrite)
        out.close()

        transferfof.write("{}\n".format(os.path.join(args.workingDir, tum, "peptides.fa")))
        analysisfile.write("{}\tsource /home/sbrown/bin/pythonvenv/python3/bin/activate;{} -tdir tmpdirXXXXXX -a {} -f peptides.fa > bindingRes.pMHC;\n".format(tum, NETMHCPAN_BIN, hlas))

    transferfof.close()
    analysisfile.close()


    idsfile = open(os.path.join(args.workingDir, "peptideIDs.tsv"), "w")
    idsfile.write("tumour\tminPeptideID\tENST\tmutation\tvarPos\tref_mut\tpepSequence\trefProtSeq\n")
    for tum in pepsInfo:
        for k in pepsInfo[tum]:
            idsfile.write("{}\t{}\t{}".format(tum,k,pepsInfo[tum][k]))
    idsfile.close()

    print("Number of mutations skipped: {}".format(numMutationsSkipped))
    print("Number of mutations processed successfully: {}".format(numMutationsComplete))

    print("done.")

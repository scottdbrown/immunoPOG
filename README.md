# immunoPOG
immunoPOG scripts

August 10, 2017
---------------

Scripts have been in development over the past couple of months.
First commit represents a "final" version of them, cleaned up.

Usage:
------

A python virtual environment is used from here: `/projects/sbrown_prj/bin/pythonvenv/python3`.

These scripts are designed to be able to be run on sets of samples. Required input is a text file with POG ids, one per line:
```
POG001
POG002
...
```

# Step 1: Generate file of files for sequence data and mutation data
```
    (python3) [sbrown@gphost03 /projects/sbrown_prj/POG/immunoPOG/data]
    λ python ../scripts/collectFilesForImmunoPOG.py -fq immunoPOGs_170421_fastqs.tsv immunoPOGs_170421.txt immunoPOGs_170609
```
* immunoPOGs_170421.txt is the file of POG ids explained above.
* immunoPOGs_170421_fastqs.tsv is a list of .fastq.gz files for the POGs. If not provided, the script will find the .bam files.

This script generates "_mutations.tsv" and "_rnaseq.tsv" files.

# Step 2: Run OptiType on fastq files
This is covered elsewhere. This should result in a directory containing directories for each sample.

# Step 3: Generate a file containing POG sample id, HLA type, SNV, and INDEL files.
```
    (python3) [sbrown@gphost03 /projects/sbrown_prj/POG/immunoPOG/analysis/mutPepBind]
    λ python ../../scripts/tallyHLAtypes.py /projects/trans_scratch/validations/optitype/POGs ../../data/immunoPOGs_170609_mutations.tsv immunoPOG19_toRun.tsv
```
* `/projects/trans_scratch/validations/optitype/POGs` is the path to the directory containing OptiType result directories.

If there are multiple result directories for a POG id (from multiple sequence files), this will select the first choice HLA type from the file with the highest number of supporting reads.

# Step 4: For each mutation, get the minimal peptide that needs to be submitted to NetMHCpan 3.0
```
    (python3) [sbrown@gphost03 /projects/sbrown_prj/POG/immunoPOG/analysis/mutPepBind/immunoPOG19_170810] 
    λ python ../../../scripts/prepPeptideBindingPredictions_SNV_indel.py ../immunoPOG19_toRun.tsv ../../../data/Homo_sapiens.GRCh37.69.cdna.all.fa ../../../data/Homo_sapiens.GRCh37.69.pep.all.fa `pwd`
```
* `Homo_sapiens.GRCh37.69.cdna.all.fa` is the reference transcriptome and `Homo_sapiens.GRCh37.69.pep.all.fa` is the reference proteome that match the mutation calling.

*_Note_*: Mutations occuring in proteins which do not begin with "M" will be skipped, as will mutations which occur before the start codon, interrupt the start codon, do not match the reference sequence, or are splice variants.

This script sets up files needed for batch submission to the genesis.phage.bcgsc.ca cluster.

# Step 5: Running NetMHCpan predictions
This uses the utility `clusterTAS`, available here: [https://github.com/scottdbrown/bcgsc-scripts](https://github.com/scottdbrown/bcgsc-scripts) - [clusterTAS](https://github.com/scottdbrown/bcgsc-scripts/blob/master/clusterTAS). This is not required, but the output files and downstream parsing have been set up to fit with the input and output of `clusterTAS`.
Briefly, clusterTAS will create directories on Genesis for each job, and copy over any necessary files (using Apollo if requested). It will only submit `x` number of jobs, depending either on a user defined max number of jobs, or the amount of space you have available on genesis. It will monitor the jobs on Genesis, and copy back result directories when complete, and submit more jobs.

*_Note_*: When running NetMHCpan 3.0 in parallel (multiple executables at once), make sure to specify the argument `-tdir tmpdirXXXXXX` in order to avoid tmpdirs from separate executables from overwriting eachother.

```
    (python3) [sbrown@gphost03 /projects/sbrown_prj/POG/immunoPOG/analysis/mutPepBind/immunoPOG19_170810] 
    λ clusterTAS --genMem 1 --maxJobs 50 --noApollo -c filesToTransfer.fof analysisToRun.sh /genesis/extscratch/sbrown/immunoPOG/ /projects/sbrown_prj/POG/immunoPOG/analysis/mutPepBind/immunoPOG19_170810/results/
```

# Step 6: Parsing the results
```
    (python3) [sbrown@gphost03 /projects/sbrown_prj/POG/immunoPOG/analysis/mutPepBind/immunoPOG19_170810] 
    λ python ../../../scripts/parseNetMHCpanOutput.py immunoPOG19_pMHC_parsed.tsv --peptideInfoFiles peptideIDs.tsv --result_dirs results/analysis_results/
```

This will parse through the NetMHCpan output files, and cross-reference back to the `peptideIDs.tsv` file to only keep the peptide-MHC predictions which contain novel peptide sequence. It does this in two ways: 1) being aware of where the variant is in the submitted peptide and only keeping 8-11mers which contain the variant, and 2) checking to make sure each peptide is not found in the reference sequence.

* `immunoPOG19_pMHC_parsed.tsv` is the output file containing a line for each pMHC prediction, with all relevant information.
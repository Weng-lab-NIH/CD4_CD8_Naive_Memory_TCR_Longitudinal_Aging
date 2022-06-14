# DESCRIPTION

This directory contains the scripts needed to run the TCR Data Processing scripts from the large FASTQ files generated from the Illumina Sequencing Machine to the human subject-specific FASTQ file used for analysis. Here is a brief summary of the purpose of each step as designated by each subdirectory:

1. **01_Filler_Match**: These scripts identify the barcode sequence in the FASTQ files and matches it to a corresponding subject whose barcode ID is known.

2. **02_MiGEC_Reformat**: These scripts reformat the FASTQ file so that they can be read by MiGEC to identify V-J genes, CDR3 AA, and CDR3 NT, along with their corresponding UMI counts. 

3. **03_Insert_PCR_Time_UMI**: These scripts inserts a 3-letter nucleic acid tag into each UMI to designate when the sample was PCR amplified. This tag is later used for UMI Grouping and contamination removal purposes. *Note*: FASTQ files used as input files for this script need a "Batch Number" written in the file name. 

4. **04_Combine_FASTQ**: These script(s) combines all of the FASTQs associated with each human subject so that all of each subject's respective TCR sequences can be eventually read by MiGEC. 

5. **05_UMI_Grouping**: These scripts cluster similar sequences (Hamming distance criteria) together based on when they were PCR amplified. Similar sequences get assigned the earliest PCR time-tag or a new UMI tag. 

6. **06_MiGEC**: This directory contains the compiled MiGEC .JAR file (v. 1.2.7) used to read and identify the FASTQ sequences.

7. **07_Separate_Alpha_Beta**: These scripts(s) create separate files for TCRA and TCRB sequences and writes the respective TCR information to each file.

8. **08_Merge_MiGEC_Data**: This script will read the assembled consensus FASTQ file generated from MiGEC and the separated Alpha/Beta files to combine all relevant information into one file (one for alpha, one for beta).

9. **09_Contamination**: Contains the scripts needed for the Contamination process. 

   * **01_Cut_Fix_Merged_Files**: Removes TCR sequences with stop codons and less than 3 reads.

   * **02_Compact_Merged_Files**: Collapses consensus 150bp TCR sequence to the first 100 bp and collapses UMI to 12bp while identifying PCR time tag. 

   * **03_Overlap**: Identifies which TCR sequences (100 bp sequence, 12bp UMI, and CDR3AA) are found in multiple donors and writes them to the output file. These sequences are known as the "contaminant sequences".

   * **04_Find_Source**: Identifies the source of the contaminant sequences (based on the PCR time tag which tells when the sequence was amplified relative to other sequences) and writes their corresponding information to the output file. 

   * **05_Remove_Source**: Re-writes the overlap files that don't contain any of the source contaminant TCR sequences. 

   * **06_Clean_Merged_File**: Reads the clean overlap file and removes the overlap sequences from the original merged file. 

10. **10_Clean_Summarize_FASTQ**: Summarizes the Merged File (can work on clean one and original one) to collapse the data into unique V-CDR3AA-J or V-CDR3AA-CDR3NT-J sequences. Also can summarize this data by counting number of unique TCR sequences, total UMIs, and reads. 

# DISCLAIMER(S): 

1. Some of the scripts aren't properly annotated due to when they were written and/or code writing style of the author.

2. Some of these scripts were written with Python 2.x.xx and will need to be adjusted if using Python 3.x.xx. A giveaway as to whether they were written with Python 2.x.xx or 3.x.xx should be the `print` statement if it exists. 
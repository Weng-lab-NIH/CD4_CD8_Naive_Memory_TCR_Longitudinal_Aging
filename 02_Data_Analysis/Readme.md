
# DESCRIPTION
This directory contains the scripts needed to analyze the processed TCR data. Here is a brief summary of the purpose of each step as designated by each subdirectory:

1.  **01_DivE**: These scripts will reformat the TCR files (V-CDR3AA-J, UMI, reads) to be compatible for DivE. These scripts also will read in a CSV file that lists the number of cells to project the TCR data to as well as the R script to run DivE(v. 1.0) with. 

2.  **02_Overlap_Estimation_Calculations**: These scripts will calculate the percentage overlap of 2 samples of the same human subject (1st visit vs. 2nd visit, Naive vs. Memory, CD4 vs. CD8) as well as for the reproducibility data.[^1]

3. **03_Simpson_Index**: This script will calculate the Inverse Simpson's Index for each human subject's sample. [^2]

4. **04_Find_Common_TCRs**: These scripts will calculate a normalized UMI frequency for each unique TCR (V-CDR3AA-J) in all of the human subject samples as well as calculate the corresponding number of public/private TCRs for each human subject. [^3]

5. **05_Estimate_abTCR_Diversity**: These scripts fit equations to tables with TCRa, TCRb, and abTCR richness values from multiple public single-cell sources. These equations will be used to estimate abTCR diversity based on DivE estimats of TCRa and TCRb richnesses, developed in **01_DivE**.


# DISCLAIMER(S):

[^1]: The script will also calculate the Morisita-Horn and Bhattacharyya Index, but these values aren't subsequently used or displayed in the paper. 

[^2]: The script will also calculate the normal Simpson's Index and Complementary Simpson's Index, but these values aren't subsequently used or displayed in the paper. 

[^3]: 2 sets of output files will have to be generated: 1) Ouput files that ignore if a particular TCR was from the 1st or 2nd visit (these files are used for Figure 7A-7B) and 2) Output files that consider if a particular TCR was from the 1st or 2nd visit (this is important as this information is used for Figure 7C-7D in the paper). 

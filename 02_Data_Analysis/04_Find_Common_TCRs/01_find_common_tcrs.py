"""
FILE NAME:      01_find_common_tcrs.py

COMPLEMENTARY
SCRIPT(S)/
FILE(S):        submit_find_common_tcrs.sh

DESCRIPTION:    This script will list each TCR found in the pool of samples analyzed and calculate the log10 UMI
                Percentage (normalized and adjusted for each sample's repertoire) as well as list how many samples
                have that TCR.

INPUT
FILES(S):       1) Cut Fixed Files (either on CDR3AA or CDR3NT).

OUTPUT
FILES(S):       1) TCR Matrix (1 for alpha & 1 for beta): each row is a unique TCR while each column represents
                each sample's adjusted frequency for that TCR.

INPUT
ARGUMENT(S):    1) Input Directory: The directory that contains all of the fixed cut files needed to run this program.
                2) Output Directory: The directory that will contain all of the TCR matrices.
                3) TCR Type: Specify whether the TCRs are based on CDR3AA or CDR3NT.
                4) File Name: Provide the file name prefix used to name the output files.

CREATED:        07Apr19

MODIFICATIONS
LOG:
16Aug19         Adjusted and Normalized the TCR frequency to each donor's repertoire size. Present the TCR
                as the median frequency instead of the mean frequency.
16Aug19         Optimization of creation of large TCR dataframe.


LAST MODIFIED
BY:             Thomas Nguyen

PYTHON
VERSION USED
TO WRITE
SCRIPT:         3.7

VERSION:        2.0

AUTHOR(S):      Thomas Nguyen

STATUS:         Working

TO DO LIST:     None at this moment.

DISCLAIMER(S):

NOTE(S):
"""
import argparse,os, glob, re
import pandas as pd
import numpy as np
from collections import defaultdict

def make_dir(output_directory):
    """
    This function makes the output directory if it already doesn't exist
    """
    if not os.path.exists(output_directory):
        try:
            os.makedirs(output_directory)
        except OSError as exc:
            if exc.errno != os.errno.EEXIST:
                raise
            pass


def find_build_tcr_df(files,tcr_type):
    """
    This function reads in all of the fixed cut files and creates list of all possible TCRs.
    """
    all_tcr_donor_df = pd.DataFrame()
    for file in files:
        match_all = re.search(r"P7M1S-none-(\d{3,4})-(alpha|beta)combinedfixed_cut.txt",file)
        if match_all:
            donor = match_all.group(1)
        tcr_count_dict = defaultdict(int)
        with open(file) as input:
            for line in input:
                info = line.rstrip().split('\t')
                if tcr_type == 'CDR3_AA':
                    v_gene,cdr3_aa,j_gene = info[:3]
                    tcr = '{}|{}|{}'.format(v_gene,cdr3_aa,j_gene)
                    umi = int(info[3])
                elif tcr_type == 'CDR3_NT':
                    v_gene, cdr3_aa, j_gene,cdr3_nt = info[:4]
                    tcr = '{}|{}|{}|{}'.format(v_gene, cdr3_aa, j_gene,cdr3_nt)
                    umi = int(info[4])
                tcr_count_dict[tcr] += umi
        file_df = pd.DataFrame.from_dict(tcr_count_dict,orient = "index",columns = [donor])
        all_tcr_donor_df = all_tcr_donor_df.combine_first(file_df)
    return all_tcr_donor_df


def additional_processing(tcr_df):
    """
    This function builds the CSV file that lists the UMI count for each TCR from each donor.
    This file will eventually be used to generate the Manhattan/Strip/Jitter Plots in R.
    """
    tcr_df['Number of Common Donors'] = (tcr_df > 0).sum(axis=1)
    for donor in tcr_df.columns.values.tolist()[:-2]:
        tcr_df[donor] = np.log10(float(100)*tcr_df[donor]/(tcr_df[donor].sum()))
    tcr_df = tcr_df.replace(0,np.nan)
    tcr_df['Median Value'] = tcr_df.iloc[:,:-2].median(axis = 1)
    return tcr_df



def Main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input_directory",
                        help="Specify the directory that contains all of the fixed cut files that need to be analyzed.")
    parser.add_argument("output_directory",
                        help="Specify the directory where the output files will reside in"
                             "(if it doesn't exist it will be created).")
    parser.add_argument("file_type",help = "Specify what type of fixed cut file needs to be analyzed",
                        choices = ['CDR3_AA','CDR3_NT'],type = str)
    parser.add_argument("file_name",help = "Specify the File Name Prefix you want attached to the results file",
                        type = str)

    args = parser.parse_args()

    make_dir(args.output_directory)

    input_alpha_files = sorted(glob.glob(args.input_directory + '*alphacombinedfixed_cut.txt'))
    input_beta_files = sorted(glob.glob(args.input_directory + '*betacombinedfixed_cut.txt'))

    alpha_sharing_matrix = find_build_tcr_df(input_alpha_files,args.file_type)
    alpha_sharing_matrix = additional_processing(alpha_sharing_matrix)
    #Sharing Matrix
    alpha_output_name1 = '{}_Alpha_Common_TCRs.csv'.format(args.file_name)
    alpha_output_name1 = os.path.join(args.output_directory,alpha_output_name1)
    alpha_sharing_matrix.to_csv(alpha_output_name1)


    beta_sharing_matrix = find_build_tcr_df(input_beta_files, args.file_type)
    beta_sharing_matrix = additional_processing(beta_sharing_matrix)
    #Sharing Matrix
    beta_output_name1 = '{}_Beta_Common_TCRs.csv'.format(args.file_name)
    beta_output_name1 = os.path.join(args.output_directory, beta_output_name1)
    beta_sharing_matrix.to_csv(beta_output_name1)


if __name__ == '__main__':
    Main()


"""
FILE NAME:      Calculate_Overlap_Statistics_01.py

COMPLEMENTARY
SCRIPT(S)/
FILE(S):        compile_overlap_statistics_results_02.py; generate_calculate_overlap_statistics_01a.py

COMPATIBLE
PIPELINE(S):    Human

DESCRIPTION:    This script will calculate the CDR3 AA and CDR3 NT percentage overlap by TCR/UMI, Bhattacharyya Index,
                and Morisita-Horn Index. For regular Human TCR donors, it will calculate the following types of overlap:

                1) Visit A vs. Visit B (within CD4M, CD4N, CD8M, CD8N, CD4NM, CD8NM)
                2) for CD4 or CD8, Naive vs. Memory (within Visit A or Visit B)
                3) for Naive or Memory or Naive-Memory combined CD4 vs. CD8 (within Visit A vs. Visit B)


INPUT
FILES(S):       1) Fixed Cut CDR3 NT File:s 6-column tab-delineated file that lists V-Gene, CDR3 AA, J-Gene,
                CDR3 NT, UMI counts, and Reads counts.

OUTPUT
FILES(S):       1) A .csv file that lists the files that were compared, the type of comparison, and the Percent Overlap,
                Morisita-Horn Overlap Index, and Bhattacharyya Overlap Index.

INPUT
ARGUMENT(S):    1) Input Directory: The directory that houses the Fixed Cut CDR3 NT Files.
                2) Output Directory: The directory where the results .csv file will be written to.

CREATED:        05Apr19

MODIFICATIONS
LOG:
19Aug19         Optimization of the creation of the overlap TCR dataframe.

PYTHON
VERSION USED
TO WRITE
SCRIPT:         3.7

VERSION:        2.0

AUTHOR(S):      Thomas Nguyen

STATUS:         Working

TO DO LIST:     None at this moment
"""
import argparse, re, os, glob
from collections import defaultdict
from itertools import combinations
import pandas as pd
import numpy as np

def make_dir(output_directory):
    """
    This function creates the output directory if it doesn't exist and prevents
    concurrent directory creations.
    """
    if not os.path.exists(output_directory):
        try:
            os.makedirs(output_directory)
        except OSError as exc:
            if exc.errno != os.errno.EEXIST:
                raise
            pass

def parse_filename(filename,file_type):
    """
    This function parses the name of a file to obtain the donor, timepoint, cell type, chain, and cell subset.
    The regular expression search was made and tested through https://regex101.com/r/HcLY7j/4.
    """
    if file_type != 'reproducibility':
        match = re.search(r"P7M1S-(.*)-(\d{3,4})([AB])CD(4|8|48)(N|M|NM)-(alpha|beta)combinedfixed_cut.txt",filename)
        if match:
            donor,timepoint,cell_type,cell_subset,chain = match.group(2),\
                                          match.group(3),\
                                          'CD{}'.format(match.group(4)),\
                                          match.group(5),match.group(6)
    elif file_type == 'reproducibility':
        match = re.search(r"P7M1S-(.*)-(.*)-(alpha|beta)combinedfixed_cut.txt",filename)
        donor,timepoint,cell_type,cell_subset,chain = match.group(2),np.nan, np.nan, np.nan, match.group(3)
    return donor,timepoint,cell_type,cell_subset,chain


def sort_files_by_donor(input_directory,file_type):
    """
    This function creates a combination of all possible valid comparisons (Same donor and same chain for
    non-reproducibility; same chain for reproducibility) for overlap analysis.
    """
    all_file_comparisons = list()
    input_files = sorted(glob.glob(input_directory + '*combinedfixed_cut.txt'))
    combination_list = list(combinations(input_files,2))
    for combo in combination_list:
        donor_1, timepoint_1, cell_type_1, subset_1, chain_1 = parse_filename(combo[0],file_type)
        donor_2, timepoint_2, cell_type_2, subset_2, chain_2 = parse_filename(combo[1],file_type)
        if file_type != 'reproducibility':
            if donor_1 == donor_2 and chain_1 == chain_2:
                file1 = combo[0]
                file2 = combo[1]
                all_file_comparisons.append((file1,file2))
        else:
            if chain_1  == chain_2:
                file1 = combo[0]
                file2 = combo[1]
                all_file_comparisons.append((file1, file2))
    return all_file_comparisons


def read_file(file):
    """
    This function parses a Fixed Cut File (based on CDR3 NT) and creates a dictionary of all TCR's by either CDR3-AA or
    CDR3-NT and tallies their UMI counts.
    """
    total_umi_counts = 0
    tcr_cdr3_aa = defaultdict(int)
    tcr_cdr3_nt = defaultdict(int)
    with open(file) as input:
        for line in input:
            info = line.rstrip().split('\t')
            v_gene,cdr3_aa,j_gene,cdr3_nt = info[:4]
            umi = int(info[4])
            aa = '{}|{}|{}'.format(v_gene,cdr3_aa,j_gene)
            nt = '{}|{}|{}|{}'.format(v_gene,cdr3_aa,cdr3_nt,j_gene)
            tcr_cdr3_aa[aa] += umi
            tcr_cdr3_nt[nt] += umi
            total_umi_counts += umi
    return tcr_cdr3_aa,tcr_cdr3_nt,total_umi_counts


def find_comparison(file_1,file_2,file_type):
    """
    This function determines the type of overlap comparison being conducted.
    """
    donor_1, timepoint_1, cell_type_1, subset_1, chain_1 = parse_filename(file_1,file_type)
    donor_2, timepoint_2, cell_type_2, subset_2, chain_2 = parse_filename(file_2,file_type)
    if file_type != 'reproducibility':
        # Timepoint comparison of same cell type & subset
        if timepoint_1 != timepoint_2 and cell_type_1 == cell_type_2 and subset_1 == subset_2:
            comparison = 'CD{}{} A vs. B'.format(cell_type_1,subset_1)
        # Comparison of Naive vs. Memory within CD4/CD8
        elif timepoint_1 == timepoint_2 and cell_type_1 == cell_type_2 and subset_1 != subset_2:
            comparison = 'CD{} N vs. M'.format(cell_type_1)
        # Comparison of CD4 vs. CD8 within Naive/Memory
        elif timepoint_1 == timepoint_2 and cell_type_1 != cell_type_2 and subset_1 == subset_2:
            comparison = '{} CD4 vs. CD8'.format(subset_1)
        else:
            comparison = None
    elif file_type == 'reproducibility':
        comparison = 'Reproducibility Comparison'
    return comparison


def build_tcr_ol_dataframe(tcr_dict_1,tcr_dict_2):
    """
    This function builds the TCR Dataframe that lists the UMI counts associated with each TCR found
    in the file pairs.
    """
    tcr_df_1 = pd.DataFrame.from_dict(tcr_dict_1,orient = 'index',columns = ['TCR Count 1'])
    tcr_df_2 = pd.DataFrame.from_dict(tcr_dict_2,orient = 'index',columns = ['TCR Count 2'])
    tcr_df = tcr_df_1.combine_first(tcr_df_2)
    return tcr_df


def calculate_percentage_overlap(tcr_umi_df,total_umi_count_1,total_umi_count_2):
    """
    This function calculates the percentage overlap (based on unique TCRs and UMI Count) between
    the 2 samples based on the following definition:

    OL_TCR% = 100*TCR_OL/(TCR_1+TCR_2-TCR_OL)
    OL_UMI% = 100*(UMI_OL_1 + UMI_OL_2)/(UMI_1+UMI_2)
    """
    ol_tcrs_df = tcr_umi_df[(tcr_umi_df['TCR Count 1'] > 0)&
                              (tcr_umi_df['TCR Count 2'] > 0)]
    tcr_1_df = tcr_umi_df[tcr_umi_df['TCR Count 1'] > 0]
    tcr_2_df = tcr_umi_df[tcr_umi_df['TCR Count 2'] > 0]

    ol_tcr_count = ol_tcrs_df.shape[0]
    ol_umi_1 = ol_tcrs_df['TCR Count 1'].sum()
    ol_umi_2 = ol_tcrs_df['TCR Count 2'].sum()

    ol_tcr_per = 100*float(ol_tcr_count) / (tcr_1_df.shape[0] + tcr_2_df.shape[0]-ol_tcr_count)
    ol_umi_per = 100*float(ol_umi_1+ol_umi_2) / (total_umi_count_1+total_umi_count_2)
    return ol_tcr_per,ol_umi_per


def morisita_horn(overlap_dataframe,umi_tot_1,umi_tot_2):
    """
    This function calculates the Morisita-Horn Similarity Index based on 2 Samples.
    This function is based on overlapping UMI percentage.
    Source of Equation:
    Horn, H.S. (1966). Measurement of "Overlap" in comparative ecological studies.
    The American Naturalist 100:419-424.
    """
    ####### Method 1 ##########
    mh_index_df = overlap_dataframe
    mh_index_df['Numerator-MH'] = (mh_index_df['TCR Count 1']/umi_tot_1)  * (mh_index_df['TCR Count 2']/umi_tot_2)
    mh_index_df['Denominator-MH'] = ((mh_index_df['TCR Count 1'] / umi_tot_1)**2) + ((mh_index_df['TCR Count 2'] / umi_tot_2)**2)
    numerator = mh_index_df['Numerator-MH'].sum()
    denominator = mh_index_df['Denominator-MH'].sum()
    ####### Method 2 #########
    # numerator = 0
    # denominator = 0
    # for row in overlap_dataframe.itertuples():
    #     tcr,umi_count_1,umi_count_2 = row[:3]
    #     x1 = float(umi_count_1)/umi_tot_1
    #     y1 = float(umi_count_2)/umi_tot_2
    #     numerator += x1*y1
    #     denominator ((x1**2) + (y1**2))
    try:
        mh_index = 2*numerator/denominator
    except ZeroDivisionError:
        mh_index = 0
    return mh_index

def bhattacharyya(overlap_dataframe,umi_tot_1,umi_tot_2):
    """
        This function calculates the Bhattacharyya Similarity Coefficient based on 2 samples;
        this function is based on overlapping UMI percentages.
        Source of Equation:
        Bhattacharyya, A. (1943). "On a measure of divergence between two statistical populations
        defined by their probability distributions". Bulletin of the Calcutta Mathematical Society.
        35: 99-109.
        """
    bhat_df = overlap_dataframe

    ##### Method 1 ############
    bhat_df['Bhat-SQRT'] = ((bhat_df['TCR Count 1']/umi_tot_1) * (bhat_df['TCR Count 2']/umi_tot_2)).pow(1./2)
    bhat_coeff = bhat_df['Bhat-SQRT'].sum()

    ###### Method 2 #########
    # bhat_coeff = 0
    # for row in bhat_df.itertuples():
    #     tcr,umi_count_1,umi_count_2 = row[:3]
    #     x1 = float(umi_count_1)/umi_tot_1
    #     y1 = float(umi_count_2)/umi_tot_2
    #     geo_mean = math.sqrt(x1*y1)
    #     bhat_coeff += geo_mean
    return bhat_coeff


def Main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input_directory",
                        help="The directory that contains the Fixed Cut Files (V-CDR3 AA-J-CDR3 NT)")
    parser.add_argument("output_directory",
                        help="Parent Directory where summary Files (1 for each donor) will reside in")
    parser.add_argument("file_type",
                        help = 'Specify whether the type of Fixed Cut files to be analyzed',
                        choices = ['donor','reproducibility'])

    args = parser.parse_args()
    make_dir(args.output_directory)
    matched_files = sort_files_by_donor(args.input_directory,args.file_type)
    all_data_df = pd.DataFrame()
    comparison_columns = ['Sample 1 Full Name', 'Sample 1 Donor', 'Sample 1 Timepoint', 'Sample 1 Cell Type',
                          'Sample 1 Subset', 'Sample 1 Chain', 'Sample 2 Full Name', 'Sample 2 Donor',
                          'Sample 2 Timepoint', 'Sample 2 Cell Type', 'Sample 2 Subset', 'Sample 2 Chain',
                          'Comparison Type', 'Sample 1 TCR Count by AA', 'Sample 1 TCR Count by NT',
                          'Sample 1 UMI Count', 'Sample 2 TCR Count by AA', 'Sample 2 TCR Count by NT',
                          'Sample 2 UMI Count','Percentage Overlap by AA-TCR', 'Percentage Overlap by AA-UMI',
                          'Percentage Overlap by NT-TCR','Percentage Overlap by NT-UMI',
                          'Morisita-Horn Index by AA','Morisita-Horn Index by NT',
                          'Bhattacharyya Index by AA', 'Bhattacharyya Index by NT']

    for file_pair in matched_files:
        donor_1, timepoint_1, cell_type_1, subset_1, chain_1 = parse_filename(file_pair[0], args.file_type)
        donor_2, timepoint_2, cell_type_2, subset_2, chain_2 = parse_filename(file_pair[1], args.file_type)
        file_comparison = find_comparison(file_pair[0],file_pair[1],args.file_type)
        if not file_comparison:
            continue

        #obtain a dictionary of unique TCR's by CDR3 AA and CDR3 NT, total UMI counts
        f1_tcr_cdr3_aa,f1_tcr_cdr3_nt,f1_total_umi = read_file(file_pair[0])
        f2_tcr_cdr3_aa, f2_tcr_cdr3_nt, f2_total_umi = read_file(file_pair[1])
        ol_dataframe_aa = build_tcr_ol_dataframe(f1_tcr_cdr3_aa,f2_tcr_cdr3_aa)
        ol_dataframe_nt = build_tcr_ol_dataframe(f1_tcr_cdr3_nt,f2_tcr_cdr3_nt)

        ############# Percentage Overlap ##################
        ol_tcr_aa,ol_umi_aa = calculate_percentage_overlap(ol_dataframe_aa,f1_total_umi,f2_total_umi)
        ol_tcr_nt, ol_umi_nt = calculate_percentage_overlap(ol_dataframe_nt, f1_total_umi, f2_total_umi)
        overlap_tcr_df_aa = ol_dataframe_aa[(ol_dataframe_aa['TCR Count 1'] > 0) &
                                            (ol_dataframe_aa['TCR Count 2'] > 0)]
        overlap_tcr_df_nt = ol_dataframe_nt[(ol_dataframe_nt['TCR Count 1'] > 0) &
                                            (ol_dataframe_nt['TCR Count 2'] > 0)]

        ############# Morisita-Horn Index ####################
        mh_aa = morisita_horn(overlap_tcr_df_aa,f1_total_umi,f2_total_umi)
        mh_nt = morisita_horn(overlap_tcr_df_nt, f1_total_umi, f2_total_umi)

        ############## Bhattacharyya Index ###################
        bh_aa = bhattacharyya(overlap_tcr_df_aa, f1_total_umi, f2_total_umi)
        bh_nt = bhattacharyya(overlap_tcr_df_nt, f1_total_umi, f2_total_umi)
        sample_1_name = '{}{}{}{}-{}'.format(donor_1, timepoint_1, cell_type_1, subset_1, chain_1)
        sample_2_name = '{}{}{}{}-{}'.format(donor_2, timepoint_2, cell_type_2, subset_2, chain_2)
        ol_results_row = [[sample_1_name,donor_1,timepoint_1,cell_type_1,subset_1,chain_1,
                           sample_2_name,donor_2, timepoint_2, cell_type_2, subset_2, chain_2,
                           file_comparison,len(f1_tcr_cdr3_aa.keys()),len(f1_tcr_cdr3_nt.keys()),f1_total_umi,
                           len(f2_tcr_cdr3_aa.keys()), len(f2_tcr_cdr3_nt.keys()), f2_total_umi,
                           ol_tcr_aa,ol_umi_aa,ol_tcr_nt,ol_umi_nt,mh_aa,mh_nt,bh_aa,bh_nt]]

        ol_results_df = pd.DataFrame(ol_results_row,columns = comparison_columns)
        all_data_df = pd.concat([all_data_df,ol_results_df],ignore_index= True)

    output_filename = 'Donor_{}_OL_Summary.csv'.format(donor_1)
    output_filename = os.path.join(args.output_directory,output_filename)

    all_data_df.to_csv(output_filename)


if __name__ == '__main__':
    Main()
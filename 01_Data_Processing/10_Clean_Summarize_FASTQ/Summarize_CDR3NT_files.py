"""
FILE NAME:      generate_umi_grouping_part1.py

COMPLEMENTARY
SCRIPT(S)/
FILE(S):        umi_grouping_part1.py; submit_umi_regrouping1.sh

DESCRIPTION:    This script will generate a submission script for each combined file recognized containing Unix
                commands to schedule a job for the file.

INPUT
FILES(S):       1) Combined Files

OUTPUT
FILES(S):       1) A submission script (no extension) for each combined file

INPUT
ARGUMENT(S):    1) Input Directory: The directory in which the DivE Summary Files are located in.
                2) Output Directory: The directory where the files
                3) Input Compression: Specify whether the combined files are compressed or not
                4) Output Compression: Specify whether the output file will be compressed or not
                5) Verbosity (optional): Used primarily for debugging purpses (will list TCR & UMI sequence, cluster
                number, etc.)

CREATED:        07Mar18

MODIFICATIONS
LOG:
04Jun18         1) Added code to check if submission script already exists (prevents duplicate
                submission script generation)
25Jun18         2) Added code to check if FASTQ file from input directory had already been converted
                to a summary file(prevents duplicate submission script generation).

LAST MODIFIED
BY:             Thomas Nguyen

PYTHON
VERSION USED
TO WRITE
SCRIPT:         2.7.14

VERSION:        2.1

AUTHOR(S):      Thomas Nguyen

STATUS:         Working

TO DO LIST:     None at this moment.

DISCLAIMER(S):  The gzip python module is poorly optimized to create large compressed files so run times for
                writing compressed files (especially large ones) will take longer than usual; if need be, write
                the output file as uncompressed but compress the file after processing is done.
"""

import glob, argparse, re, sys, os
from collections import namedtuple

tcr_sample = namedtuple("tcr_sample",['barcode','donor','chain','tcr_cdr3_aa','tcr_cdr3_nt','umi','reads'])

def make_dir(directory):
    """
    This function creates the directory if it doesn't already exist.
    """
    if not os.path.exists(directory):
        try:
            os.makedirs(directory)
        except OSError as exc:
            if exc.errno != os.errno.EEXIST:
                raise
            pass


def parse_file(file):
    tcrs_by_aa = set()
    tcrs_by_nt = set()
    umi_counts_tot = 0
    reads_count_tot = 0
    match = re.search(r"P7M1S-(.*)-(.*)-(alpha|beta)combinedfixed_cut.txt",os.path.basename(file))
    if match:
        barcode,donor,chain = match.group(1),match.group(2),match.group(3)
    else:
        print("Error! File {} does not match the naming conventions defined by the script please change it.".format(file))
        sys.exit(1)
    with open(file) as input:
        for line in input:
            info = line.rstrip().split('\t')
            v_gene,cdr3_aa,j_gene,cdr3_nt = info[:4]
            umi,reads = int(info[4]),int(info[5])
            tcr_aa = '{}|{}|{}'.format(v_gene,cdr3_aa,j_gene)
            tcr_nt = '{}|{}|{}|{}'.format(v_gene,cdr3_aa,j_gene,cdr3_nt)
            tcrs_by_aa.add(tcr_aa)
            tcrs_by_nt.add(tcr_nt)
            umi_counts_tot += umi
            reads_count_tot += reads
    tcr_summary = tcr_sample(barcode,donor,chain,len(tcrs_by_aa),len(tcrs_by_nt),umi_counts_tot,reads_count_tot)
    return tcr_summary



def Main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input_directory", help="specify the directory where the CDR3 NT Fixed Cut Files are located.")
    parser.add_argument("output_directory", help='Specify the location where Summary file file will reside in.')
    args = parser.parse_args()

    make_dir(args.output_directory)
    all_cut_files_cdr3_nt = sorted(glob.glob(args.input_directory + '*fixed_cut.txt'))

    output_filename = 'CDR3_NT_Summary_Counts.txt'
    output_filename = os.path.join(args.output_directory,output_filename)
    outfile = open(output_filename,'w')
    header = (('{}\t')*6 + '{}\n').format('Bardcode',"Donor","Chain","TCRs by CDR3 AA","TCRs by CDR3 NT","UMI Counts",
                                          "Reads Count")
    outfile.write(header)
    for file in all_cut_files_cdr3_nt:
        tcr_summary = parse_file(file)
        output_line = (('{}\t')*6 + '{}\n').format(tcr_summary.barcode,tcr_summary.donor,tcr_summary.chain,
                                                   tcr_summary.tcr_cdr3_aa,tcr_summary.tcr_cdr3_nt,tcr_summary.umi,
                                                   tcr_summary.reads)
        outfile.write(output_line)
    outfile.close()

if __name__ == '__main__':
    Main()

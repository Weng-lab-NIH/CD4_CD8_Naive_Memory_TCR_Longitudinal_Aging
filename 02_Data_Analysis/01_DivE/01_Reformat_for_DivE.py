"""
FILE NAME:      01_Reformat_for_DivE.py

COMPLEMENTARY
SCRIPT(S)/
FILE(S):        None

DESCRIPTION:    This script will convert the fixed cut files to a format that DivE can recognize:
                a tab-delineated file with unique species in the first column and frequency of occurrence (UMI)
                in the second column.

INPUT
FILES(S):       1) Fixed Cut Files

OUTPUT
FILES(S):       1) A DivE-reformatted file.

INPUT
ARGUMENT(S):    1) Input Directory: The directory in which the DivE Summary Files are located in.
                2) Output Directory: The directory where the files

CREATED:        19Apr18

MODIFICATIONS
LOG:            None at this moment.

PYTHON
VERSION USED
TO WRITE
SCRIPT:         2.7.14 --> 3.7

VERSION:        1.0

AUTHOR(S):      Thomas Nguyen

STATUS:         Working

TO DO LIST:     None at this moment
"""

import argparse,re,os,time,glob
from collections import defaultdict


def reformat_file(file):
    """
    This function reads in the fixed cut files and concatenates the TCR (V Gene, CDR3 AA, J Gene)
    and keeps track of its associated cumulative UMI Count (if the duplicate TCR's haven't been removed from the fixed
    cut file yet).
    """
    tcr_umi_count = defaultdict(int)
    with open(file) as input:
        for line in input:
            info = line.rstrip().split('\t')
            v_gene,cdr3_aa,j_gene,umi_count = info[0],info[1],info[2],int(info[3])
            reformatted_tcr = '{}|{}|{}'.format(v_gene,cdr3_aa,j_gene)
            tcr_umi_count[reformatted_tcr] += umi_count
    return tcr_umi_count


parser = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("input_directory",
                    help = "The Directory that contains the Cut & Fix Files")
parser.add_argument("output_directory",
                    help = "The Directory that contains the Fixed Cut Files")
args = parser.parse_args()

if not args.input_directory.endswith("/"):
    print("Error! You must end your path arguments with a \"/\"!")
    exit(1)

if not os.path.exists(args.output_directory):   #Create the directory if it doesn't exist already
    try:
        os.makedirs(args.output_directory)
    except OSError as exc:
        if exc.errno != os.errno.EEXIST:
            raise
    time.sleep(1)

input_files = sorted(glob.glob(args.input_directory + "P7M1S-*combinedfixed_cut.txt"))

if not input_files:
    print ("Error! No Fixed Cut Files have been detected! Please change the name of the fixed cut files" \
          "to the following format: P7M1S-<unique identifier>combinedfixed_cut.txt or change the " \
          "glob function in the script!")
    exit(1)

for file in input_files:
    try:
        match = re.search("P7M1S-(.*)combinedfixed_cut.txt",file)
        unique_file_subset = match.group(1)
    except AttributeError:
        print ("Error! Your Fixed Cut Files are not in the correct format! Please fix the name of the following " \
              "fixed cut file: {} OR change the regular expression search in the script!".format(file))
    new_file_name = 'g-P7M1S-{}_reformatted.txt'.format(unique_file_subset)
    outfilename = os.path.join(args.output_directory,new_file_name)
    outfile = open(outfilename,'w')

    tcr_umi_counts = reformat_file(file)

    for tcr,umi_freq in tcr_umi_counts.items(): #Access all of the Unique TCR's Found
        output_line = '{}\t{}\n'.format(tcr,umi_freq)   #Format is a 2 tab-delineated file
        outfile.write(output_line)
    outfile.close()


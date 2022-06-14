"""
FILE NAME:      Cut_Fixed_Files.py

COMPLEMENTARY
SCRIPT(S)/
FILE(S):        None

COMPATIBLE
PIPELINE(S):    Human, Mouse

DESCRIPTION:    This script will remove TCRs (based on V-CDR3 AA-J) below a specified read threshold (normally 3 reads).

INPUT
FILES(S):       1) Fixed File

OUTPUT
FILES(S):       1) A Cut Fixed File

INPUT
ARGUMENT(S):    1) Input Directory: The directory that contains the fixed files.
                2) Read Threshold: Specify the minimum number of reads a TCR needs

CREATED:        May2016

MODIFICATIONS
LOG:
07Aug19			Modified script to be compatible with Python 3.7


LAST MODIFIED
BY:             Thomas Nguyen

PYTHON
VERSION USED
TO WRITE
SCRIPT:         2.7.14 --. 3.7

VERSION:        3.0

AUTHOR(S):      Annette Ko; Jay Sharma; Chen Chen; Thomas Nguyen

STATUS:         Working

TO DO LIST:     None at this moment.

DISCLAIMER(S):

NOTE(S):
"""

import glob
import os
import sys
import datetime

input_dir = sys.argv[1]  # Specify the input directory containing the files you wish to cut.
output_dir = input_dir + 'cut_files'  # The output directory can be changed if desired.
min_reads = int(sys.argv[2])  # Specify the minimum number of TCR reads you wish to cut. This is usually 3.
input_files = sorted(glob.glob(input_dir + 'P7M1S-??-*combinedfixed.txt'))  # Look for fixed files within the input directory.

if not input_dir[-1] == "/":  # If the use doesn't add the slash at the end of the argument, warn them.
    print ("Error! Make sure to end your Python argument path with a /!")
    sys.exit(1)

if not input_files:
    print ('Error! No fixed files found. Check the input directory and/or check the {} script to make sure the name search matches.'.format(
        os.path.basename(__file__)))
    sys.exit(1)

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

summary_filename = 'results_cut.txt'
summary_fullpath = os.path.join(output_dir, summary_filename)
summary_file = open(summary_fullpath, 'w')

current_date = datetime.date.today().strftime("%A, %B %d, %Y")
current_time = datetime.datetime.now().strftime("%I:%M:%S %p")
summary_file.write("Cutting procedure run on {} at {}.\nEach line is one UMI and one TCR.\n".format(current_date, current_time))

for in_file in input_files:

    input_linecount = 0  # Keep track of the number of lines for every file before cutting and after cutting.
    output_linecount = 0

    input_filename = os.path.split(in_file)[1]  # Read in the current filename and store it.
    input_shortname = os.path.splitext(input_filename)[0]
    input_extension = os.path.splitext(input_filename)[1]

    output_filename = input_shortname + '_cut' + input_extension  # Write the new file using the original filename + "_cut".
    output_fullpath = os.path.join(output_dir, output_filename)
    output_file = open(output_fullpath, 'w')

    with open(in_file) as f:
        print ('Cutting file {} for reads greater than or equal to {}.'.format(input_filename, min_reads))
        for line in f:
            input_linecount += 1
            line = line.rstrip('\n')
            columns = line.split('\t')  # The columns are separated by tabs, so split each line by the tab character.
            v = columns[-5]
            cdr3_aa = columns[-4]
            j = columns[-3]
            umi_reads = int(columns[-2])
            tcr_reads = int(columns[-1])
            if tcr_reads >= min_reads:  # If and only if the number of TCR reads is greater OR EQUAL TO than the minimum threshold, write that line.
                output_file.write(line + '\n')
                output_linecount += 1

    output_file.close()
    summary_file.write(
        '\nFilename: {}\nNumber of lines in original file: {}\nNumber of lines in new file: {}\nNumber of lines removed: {}\n'.format(
            input_filename, input_linecount, output_linecount, (input_linecount - output_linecount)))
summary_file.close()

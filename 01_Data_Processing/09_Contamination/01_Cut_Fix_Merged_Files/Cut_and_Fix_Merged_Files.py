"""
FILENAME:       Cut_and_Fix_Merged_Files.py
DESCRIPTION:    Takes merged-*.txt files and outputs a file in the same format but devoid
                of stop codons, bad sequences & sequences with total reads less
                than a user-specified value (usually 3).
INPUT:          A directory containing files with the name format "merged-P7M1S-??-*-*combined*.txt".
OUTPUT:         Files in a subdirectory called "cut_fixed_merged_files" with the same input name.
CREATED:        October 10, 2017.
LAST MODIFIED:  April 5, 2018.
MODIFICATIONS:  Fixed a bug with glob wildcards and regular expression search. Added version information.
VERSION:        5.0
AUTHOR(S):      Jay Sharma
STATUS:         Completed and working.
"""

from __future__ import print_function

import argparse
import os
import re
from collections import defaultdict
from datetime import datetime
from glob import glob
from sys import exit
from time import clock


def define_outputs(path_in, dir_out):
    """Function to read input filenames and create output files.
    The output file extension will change depending on what mode (--check)
    the program is run in.
    """
    base_in = os.path.split(path_in)[1]
    name_in, extension_in = os.path.splitext(base_in)[0], os.path.splitext(base_in)[1]
    name_out = "{}{}".format(name_in, extension_in)
    # name_out = "{0}_cut_fixed{1}".format(name_in, extension_in)
    path_out = os.path.join(dir_out, name_out)
    write_out = open(path_out, 'w')
    return write_out


def get_info(f):
    """Function that extracts the barcode, donor & chain from each sample."""
    match = re.search(r"merged-P7M1S-(?:\w{2,})-(?:\S*)-(alpha|beta)combined.txt", f)
    if not match:
        print("Error! The regular expression name search didn't find any files. Check the {} script.".format(
            os.path.basename(__file__)))
        exit(1)
    chain = match.group(1)
    return chain


def check_cdr3(cdr3_aa):
    """Function that checks for any "?" or "*" in the CDR3 AA sequence & returns false if found."""
    if "?" in cdr3_aa or "*" in cdr3_aa:
        return False
    return True


def clean(cdr3_aa, v, chain):
    """Function that checks if the CDR3 AA is a valid AA sequence & returns false if invalid."""
    check = check_cdr3(cdr3_aa)
    if check:
        if chain == 'alpha':
            if cdr3_aa[0] != 'C':
                return False
        elif chain == 'beta':
            if cdr3_aa[0] == 'C':
                if 'BV17' in v or 'BV26' in v:
                    return False
            elif cdr3_aa[0] == 'Y':
                if 'BV17' not in v or 'BV26' not in v:
                    return False
        return True
    return False


def parse(line):
    """Function that extracts the values from each column from each line in the merged file."""
    columns = line.rstrip('\n').split('\t')
    cdr3_nt, v, j, d, cdr3_aa, umi, tcr, reads = columns[0], columns[1], columns[2], columns[3], columns[4], columns[5], \
                                                 columns[6], int(columns[-1])
    return cdr3_nt, v, j, d, cdr3_aa, umi, tcr, reads


def build_dictionary(f_in, min_reads, chain):
    """Function to build a dictionary where V + CDR3 amino acid + J is the key
    and D + CDR3 nucleotide + UMI + TCR + reads are the values.
    """
    clonotype_dict = defaultdict(list)  # Create a defaultdict to hold info for every unique V, CDR3 AA & J (clonotype).
    line_count = 0  # Record the number of lines read.

    with open(f_in) as merged_file:
        print('Cutting and fixing file {} for total reads greater than or equal to {}.'.format(os.path.basename(f_in),
                                                                                               min_reads))
        # next(merged_file)  # Skip the first line of the file, which is the column headers.
        for line in merged_file:
            line_count += 1
            cdr3_nt, v, j, d, cdr3_aa, umi, tcr, reads = parse(line)  # Extract each column value.
            if clean(cdr3_aa, v, chain):  # If the CDR3 AA sequence is clean & correct...
                key = v, cdr3_aa, j  # ...store V, CDR3 AA & J as the key...
                value = [d, cdr3_nt, umi, tcr, str(reads)]  # ...& D, CDR3 NT, UMI, TCR & reads as the value.
                clonotype_dict[key].append(value)
    return clonotype_dict, line_count


def cut_merged(f_in, dir_out, min_reads, is_header):
    """Function to perform cutting of reads below the user-specified "min_reads" threshold. Lines with the
    same V + CDR3 amino acid + J but with two or more different D + CDR3 nucleotide + UMI + TCR + reads will
    remain untouched, since removing them would result in a loss of UMIs, TCRs and reads. But if there's just
    a single D + CDR3 nucleotide + UMI + TCR + reads, the reads can be checked against the threshold.
    Discard the line if it's a single value that doesn't meet the minimum reads requirement.
    """
    count_in, count_out = 0, 0  # Record the number of lines before cutting & after cutting.
    header = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("CDR3 NT", "V", "J", "D", "CDR3 AA", "UMI", "TCR",
                                                       "Reads")  # Define the column header.

    f_out = define_outputs(f_in, dir_out)  # Define the output file name & create it.
    chain = get_info(f_in)  # Get chain from the file name.
    tcr_dict, count_in = build_dictionary(f_in, min_reads, chain)  # Compile a dictionary for the input merged file.

    if is_header:  # Write the column headers in the output file if the user specifies the "--header" argument.
        f_out.write(header)
    for key, value in tcr_dict.iteritems():
        v, cdr3_aa, j = key[0], key[1], key[2]
        for i in value:  # i becomes a list of the values (which is CDR3 NT, UMI, TCR & reads).
            d, cdr3_nt, umi, tcr, reads = i[0], i[1], i[2], i[3], i[4]
            line_out = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(cdr3_nt, v, j, d, cdr3_aa, umi, tcr, reads)
            if len(value) == 1:  # If there's only 1 CDR3 NT, UMI, TCR for a given V, J, D & CDR3 AA...
                if int(reads) >= int(min_reads):  # ...verify that reads >= minimum threshold...
                    f_out.write(line_out)  # ...& write it only if it does.
                    count_out += 1
            else:  # If there's more than 1 CDR3 NT, UMI, TCR for a given V, J, D & CDR3 AA, write them all.
                f_out.write(line_out)
                count_out += 1
    count_removed = count_in - count_out
    print("Number of lines in: {}. Number of lines out: {}. Number of lines removed: {}.\n".format(count_in, count_out,
                                                                                                   count_removed))
    tcr_dict.clear()  # Empty the dictionary before moving to the next file.
    f_out.close()


def main():
    cmdline = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    cmdline.add_argument('input_directory', help='the directory containing the input merged files')
    cmdline.add_argument('output_directory', help='the directory that will contain the cut and fixed merged files')
    cmdline.add_argument('minimum_reads', help='the minimum value of TCR reads that will be kept')
    cmdline.add_argument('-v,' '--version', action='version', version='%(prog)s 5.0')
    cmdline.add_argument('-e', '--header', help='write the column headers in the output file', action='store_true',
                         default=False)
    args = cmdline.parse_args()

    # If the user doesn't put a "/" at the end of the "input_directory" argument or
    # the end of the "output_directory" argument, warn them & exit the program.
    if not args.input_directory.endswith("/") or not args.output_directory.endswith("/"):
        print("Error! You must end your path arguments with a \"/\"!")
        exit(1)

    # Record time & date of program start.
    start_time = clock()
    print("{} program started on {}.\n".format(os.path.basename(__file__),
                                               datetime.now().strftime("%A, %B %d at %I:%M:%S %p")))

    # Make a list of all matching input fastq files.
    input_files = sorted(glob(args.input_directory + "merged-P7M1S-*-*-*combined*.txt"))

    # Check if the list of input files is empty. If it is, warn the user & exit the program.
    if not input_files:
        print("Error! No input \"merged-P7M1S-??-*-*combined*.txt\" files found. Check the following:\n"
              "1. The input directory.\n"
              "2. The {} script.\n".format(os.path.basename(__file__)))
        exit(1)

    # If the output directory doesn't exist, create it.
    if not os.path.exists(args.output_directory):
        os.makedirs(args.output_directory)

    # Do computation for every merged file.
    for merged_file in input_files:
        cut_merged(merged_file, args.output_directory, args.minimum_reads, args.header)

    # Record time & date of program completion.
    print("{} program finished on {}.".format(os.path.basename(__file__),
                                              datetime.now().strftime("%A, %B %d at %I:%M:%S %p")))
    print("Time taken: {} seconds.".format(clock() - start_time))


if __name__ == "__main__":
    main()

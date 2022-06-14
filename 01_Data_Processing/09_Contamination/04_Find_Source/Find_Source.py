"""
FILENAME:       Find_Source.py
DESCRIPTION:    This script reads in a directory containing overlap files and first concatenates them.
                The UMI column and TCR column is then used as the key for a dictionary, and the CDR3
                amino acid, UMI reads, total reads, PCR time, and file name columns are used as the value.
                Next, the values (for each UMI-TCR key) are grouped by PCR time. Whichever value has the
                smallest PCR time is considered the "source" PCR time for that UMI-TCR pair across all
                overlap files, and is written to a temporary file-like object (an "intermediate file").
                Lastly, any duplicate lines are removed from the intermediate file, and the real file is written.
                The output is a file that contains all "source" PCR times (and associated values, like CDR3
                amino acid, total reads, etc.) for a given UMI-TCR pair.
INPUT:          A directory containing files with the name format "ol-compact-merged-*.txt".
OUTPUT:         A txt file that contains the "source" values for a UMI-TCR pair.
CREATED:        February 6, 2018.
LAST MODIFIED:  March 28, 2018.
MODIFICATIONS:  Created separate parse_overlap function. Minor edits to comments.
VERSION:        2.0
AUTHOR(S):      Jay Sharma
STATUS:         Completed and working.
"""

from __future__ import print_function

import argparse
import os
from collections import defaultdict
from datetime import datetime
from fileinput import input
from glob import glob
from hashlib import md5
from itertools import groupby
from operator import itemgetter
from sys import exit
from time import clock

try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO


def define_outputs(dir_out):
    """Function to read input filenames and create output files. Returns an opened file."""
    name_out = "source_files.txt"
    path_out = os.path.join(dir_out, name_out)
    write_out = open(path_out, 'w')
    return write_out


def compute_md5(my_string):
    """Function to convert a string to a MD5 hash. Returns the hexadecimal form of the original string."""
    m = md5()
    m.update(my_string.encode('utf-8'))  # The string must first be converted to unicode.
    return m.hexdigest()


def remove_duplicates(f_temp, f_out):
    """Function to remove duplicate lines from a file. Each line from the input file or file-like object
    is read, converted to a MD5 hash, and inserted into a set. If the line hasn't been seen before, it's
    written to the output file; otherwise, it's not written. Converting to a hash saves memory when looking
    up the string in the set, as opposed to looking up an entire line in the set.

    Refer to [https://stackoverflow.com/a/48673557] for information on iterating through a file-like object.
    """
    lines_seen = set()  # Build a set to hold lines already seen.
    f_temp.seek(0)  # "Open" the temporary file-like object for reading.
    for line in f_temp:  # Iterate through the file-like object to check for duplicate lines.
        md5_line = compute_md5(line)  # Convert the line to a MD5 hash.
        if md5_line not in lines_seen:  # Not a duplicate line.
            f_out.write(line)
            lines_seen.add(md5_line)
    f_out.close()  # Close the real output file.
    f_temp.close()  # Close the temporary file-like object.


def compute_sources(list_in):
    """Function that reads in every CDR3 AA, number of UMIs, number of reads, PCR time & filename
    (a list of lists called "list_in"; a sublist) for a TCR and determines how many different
    filenames have that sublist. It does so by first grouping each sublist by filename (the last
    element in the sublist). Then the len() function is called to determine the number of unique filenames.
    """
    list_in.sort(key=lambda x: (x[-1]))  # First sort the list by earliest PCR time.
    grouped_list = [(filename, list(group)) for filename, group in groupby(list_in, key=lambda x: (x[-1]))]
    filename_count = len(grouped_list)
    return filename_count


def get_earliest(list_in):
    """Function that reads in every CDR3 AA, number of UMIs, number of reads, PCR time & filename (a list of lists
    called "list_in"; a sublist) for a TCR and determines the earliest PCR time. It does so by first grouping each
    "list_in" sublist by PCR time (second-from-last element in the sublist). Then the min() function is called
    to obtain the "list_in" group with the smallest PCR time. If no "list_in" sublist has a smallest PCR
    time, they're all returned. The number of source lists is also computed and returned.

    Refer to [https://stackoverflow.com/a/48530763] for the methodology on groupby().
    """
    list_in.sort(key=lambda x: (x[-2]))  # First sort the list by earliest PCR time.
    grouped_list = [(pcr_time, list(group)) for pcr_time, group in groupby(list_in, key=lambda x: (x[-2]))]
    list_out = min(grouped_list, key=itemgetter(0))[1]  # itemgetter gets every "list_in" with the smallest PCR time.
    source_count = compute_sources(list_out)  # Count how many files contain "source" UMI-TCR pairs.
    # source_count = len(list_out)  # Count how many "source" UMI-TCR pairs are there across all files.
    return list_out, source_count


def parse_overlap(line):
    """Function that extracts information from every line in the temporary
    file-like object (i.e. concatenated overlap file).
    """
    columns = line.rstrip('\n').split('\t')
    tcr, umi, cdr3_aa, umi_reads, total_reads, pcr_time, file_name = columns[0], columns[1], columns[2], columns[3], \
                                                                     columns[4], columns[5], columns[6]
    return tcr, umi, cdr3_aa, umi_reads, total_reads, pcr_time, file_name


def build_dictionary(concatenated_file):
    """Function that builds a dictionary with UMI, TCR as the key and CDR3 AA, UMI counts, total reads,
    PCR time and filename as the value. This is done on the concatenated file, so it contains every
    UMI-TCR pair found across all overlap files. It returns a dictionary of lists.
    """
    umi_tcr_dict = defaultdict(list)
    line_count = 0
    for line in concatenated_file:
        line_count += 1
        tcr, umi, cdr3_aa, umi_reads, total_reads, pcr_time, file_name = parse_overlap(line)
        key = umi, tcr
        value = (cdr3_aa, umi_reads, total_reads, pcr_time, file_name)
        umi_tcr_dict[key].append(value)
    return umi_tcr_dict


def concatenate_lines(fl_in):
    """Function to join all lines together from each file in the file list."""
    overlap_files = input(fl_in)
    return overlap_files


def write_temp(f_cat):
    """Function that writes the line from the input file (the overlap file) to a temporary file-like object
    (defined as "intermediate_file"). Duplicate lines are then removed from this temporary file-like object
    before being written to the real file. The lines that are written are those lines that contain the
    earliest PCR time for that UMI-TCR pair. Column headers can also be optionally written.
    """
    intermediate_output = StringIO()  # Create a temporary file-like object.
    umi_tcr_dict = build_dictionary(f_cat)  # Build a dictionary of UMI-TCR pairs from the concatenated overlap file.
    for k, v in umi_tcr_dict.iteritems():  # Process the dictionary (with UMI-TCR pair as the key) for writing.
        umi, tcr = k[0], k[1]  # Extract each UMI and TCR from the key in the UMI-TCR dictionary.
        smallest_list, number_of_sources = get_earliest(v)  # Perform grouping and determine earliest PCR time.
        for i in smallest_list:  # i is a tuple containing: CDR3 AA, UMI count, reads, earliest PCR time & filename.
            cdr3_aa, number_umis, reads, pcr_time, file_name = i[0], i[1], i[2], i[3], i[4]
            line_out = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(tcr, umi, cdr3_aa, number_umis, reads, pcr_time,
                                                                 file_name, number_of_sources)
            intermediate_output.write(line_out)  # Write to the temporary file-like object.
    return intermediate_output


def main():
    cmdline = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    cmdline.add_argument('input_directory', help='the directory containing overlap files')
    cmdline.add_argument('output_directory', help='the directory that will contain the output of this script')
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
    print("{} program started on {} for directory {}.".format(os.path.basename(__file__),
                                                              datetime.now().strftime("%A, %B %d at %I:%M:%S %p"),
                                                              args.input_directory))

    # If the output directory doesn't exist, create it.
    if not os.path.exists(args.output_directory):
        os.makedirs(args.output_directory)

    # Make a list of all matching input fastq files.
    input_files = sorted(glob(args.input_directory + "ol-compact-merged-*.txt"))

    # Create the output file.
    output_file = define_outputs(args.output_directory)

    # Define the column headers & write it to the output file if the user specifies the "--header" argument.
    header = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("100 Nucleotide TCR", "12 Nucleotide UMI", "CDR3 Amino Acid",
                                                       "UMI Count", "Number of Reads", "Earliest PCR Time",
                                                       "File Name Containing Earliest PCR Time",
                                                       "Number of Files with Earliest PCR Time")
    if args.header:
        output_file.write(header)

    # Check if the list of input files is empty. If it is, warn the user & exit the program.
    if not input_files:
        print("Error! No input txt files found. Check the following:\n"
              "1. The input directory.\n"
              "2. The {} script.".format(os.path.basename(__file__)))
        exit(1)

    # Do main computation.
    concatenated_file = concatenate_lines(input_files)  # First, combine all overlap files in the input directory.
    intermediate_file = write_temp(concatenated_file)  # Second, find earliest time & write to a file-like object.
    remove_duplicates(intermediate_file, output_file)  # Finally, remove duplicate sources & write to the real file.

    # Record time & date of program completion.
    print("{} program finished on {} for directory {}.".format(os.path.basename(__file__),
                                                               datetime.now().strftime("%A, %B %d at %I:%M:%S %p"),
                                                               args.input_directory))
    print("Time taken: {} seconds.".format(clock() - start_time))


if __name__ == "__main__":
    main()

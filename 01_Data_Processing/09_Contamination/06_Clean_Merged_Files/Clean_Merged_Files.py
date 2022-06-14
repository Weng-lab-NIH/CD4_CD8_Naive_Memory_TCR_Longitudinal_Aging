"""
FILENAME:       Clean_Merged_Files.py
DESCRIPTION:    This script reads in an overlap file and a merged file for a given sample, and outputs
                a "cleaned" merged file that's devoid of any UMI-TCR pairs found in the overlap file.
INPUT:          A txt file with the name format "*ol-compact-merged-*.txt", also known as an overlap file.
                A txt file with the name format "merged-*.txt", also known as a merged file.
OUTPUT:         A txt with the name format "cl-merged-*.txt", also known as a cleaned merged file.
CREATED:        February 21, 2018.
LAST MODIFIED:  March 28, 2018.
MODIFICATIONS:  Minor edits to comments. Added percentage calculation and version information.
VERSION:        1.0
AUTHOR(S):      Jay Sharma.
STATUS:         Working for removal of UMI-TCR pairs. UMI-only removal code not implemented (yet, if ever).
"""

from __future__ import print_function

import argparse
import os
from datetime import datetime
from sys import exit
from time import clock, sleep


def define_outputs(path_in, dir_out):
    """Function to read input filename and create the output file."""
    base_in = os.path.split(path_in)[1]
    name_out = "cl-{}".format(base_in)
    # name_out = base_in  # Keep the exact name as the input file name.
    path_out = os.path.join(dir_out, name_out)
    write_out = open(path_out, 'w')
    return write_out


def parse_file(line, file_type):
    """Function to read a line from the overlap file or the merged file and extract the UMI and TCR.
    "file_type" determines whether the program is parsing through the overlap file or the merged file;
    this distinction needs to be made because the UMI and TCR are found on different columns depending
    on the file type. Thus, this single function can parse two types of files.
    """
    columns = line.rstrip('\n').split('\t')
    if file_type == "overlap_file":
        tcr, umi = columns[0], columns[1]
    elif file_type == "merged_file":
        tcr, umi = columns[6][:100], columns[5]  # Return only the first 100 nucleotides of the TCR.
    else:  # This is error-checking. It shouldn't ever trigger, but I added it here anyway.
        tcr, umi = None, None
        print("Error! A UMI and/or TCR wasn't found! Check both the overlap file and the merged file for this sample.")
        exit(1)
    return umi, tcr


def collapse_umi(umi):
    """Function to collapse the UMI to twelve nucleotides so it can be compared to the overlap list."""
    umi = umi[0:16]  # Take only the first 16 nucleotides as the UMI (in case the UMI is 19 nucleotides).
    collapsed_umi = umi[1:5] + umi[6:10] + umi[11:15]  # Collapse the UMI to 12 nucleotides.
    return collapsed_umi


def build_overlap_list(f_overlap):
    """Function to parse the overlap file (the file containing the contaminated UMI-TCR
    pairs for a sample) and store every UMI-TCR in a dictionary for a given sample.
    It returns a dictionary with the tuple of (barcode, donor, chain) as the key
    and a list of tuples of (UMI, TCR) as value.
    """
    overlaps = []
    with open(f_overlap) as fh:
        for line in fh:
            umi, tcr = parse_file(line, "overlap_file")  # Extract the UMI & TCR from the overlap file
            overlaps.append((umi, tcr))  # Put the UMI & TCR into a list as a tuple.
    return overlaps


def write_clean_file(f_merged, o_list, f_out):
    """Function to parse the merged file and create a new, cleaned merged file that has UMI-TCRs
    that are not found in the overlap list (list of UMI-TCRs to remove). First, the UMI and TCRs
    are extracted from every line in the merged file, and collapsed to twelve nucleotides and
    one hundred nucleotides, respectively. They're then put into a tuple, which allows them to
    be checked against the list of tuples of contaminated UMI-TCRs. If the tuples match, then
    that's a contaminated UMI-TCR pair, and it's not written to the output file.
    """
    f_merged_line_count, cl_merged_line_count = 0, 0  # Initialize counters for the merged file & clean merged file.
    with open(f_merged) as fh:
        for line in fh:
            f_merged_line_count += 1
            umi, tcr = parse_file(line, "merged_file")  # Extract the UMI & TCR from the merged file.
            umi = collapse_umi(umi)  # Collapse the UMI to 12 nucleotides.
            merged_umi_tcr = umi, tcr  # Create a tuple to check against the overlap list.
            if merged_umi_tcr not in o_list:  # The merged file UMI-TCR is not in the list of UMI-TCRs to be removed.
                cl_merged_line_count += 1
                f_out.write(line)
    lines_removed = f_merged_line_count - cl_merged_line_count
    percentage_removed = (float(lines_removed) / float(
        f_merged_line_count)) * 100  # Calculate the percentage of UMI-TCR pairs removed.
    truncated_percentage = "{0:.4f}".format(percentage_removed)  # Format the percentage to 4 decimals.
    print("Number of UMI-TCR pairs in merged file = {}\n"
          "Number of UMI-TCR pairs in clean merged file = {}\n"
          "Number of UMI-TCR pairs removed = {}\n"
          "Percentage of UMI-TCR pairs removed = {}\n".format(f_merged_line_count, cl_merged_line_count, lines_removed,
                                                              truncated_percentage))
    f_out.close()


def main():
    cmdline = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    cmdline.add_argument('merged_file', help='the merged file that contains sequences to be removed')
    cmdline.add_argument('overlap_file', help='the overlap file that contains sequences to remove')
    cmdline.add_argument('output_directory', help='the directory that will contain the clean merged files')
    cmdline.add_argument('-v,' '--version', action='version', version='%(prog)s 1.0')
    cmdline.add_argument('-e', '--header', help='write the column headers in the clean merged file',
                         action='store_true', default=False)
    cmdline.add_argument('-u', '--umi', help='remove only contaminated UMIs, not contaminated TCRs',
                         action='store_true', default=False)
    args = cmdline.parse_args()

    # If the user doesn't put a "/" at the end of the "output_directory" argument, warn them & exit the program.
    if not args.output_directory.endswith("/"):
        print("Error! You must end your path arguments with a \"/\"!")
        exit(1)

    # Record time & date of program start.
    start_time = clock()
    print("{} program started on {} for merged file {} and overlap file {}.\n".format(os.path.basename(__file__),
                                                                                      datetime.now().strftime(
                                                                                          "%A, %B %d at %I:%M:%S %p"),
                                                                                      os.path.basename(
                                                                                          args.merged_file),
                                                                                      os.path.basename(
                                                                                          args.overlap_file)))

    # If the output directory doesn't exist, create it. Add a 1-second wait period to account for race conditions.
    if not os.path.exists(args.output_directory):
        try:
            os.makedirs(args.output_directory)
        except OSError, e:
            if e.errno != os.errno.EEXIST:
                raise
            sleep(1)
            pass

    # Create the cleaned merged file from the input merged file.
    output_file = define_outputs(args.merged_file, args.output_directory)

    # Define the column headers & write it to the output file if the user specifies the "--header" argument.
    header = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("CDR3 Nucleotide", "V", "J", "D", "CDR3 Amino Acid", "UMI",
                                                       "TCR", "Number of Reads")
    if args.header:
        output_file.write(header)

    # Build a list of tuples of (UMI, TCR) from the overlap file. These are contaminated UMI-TCR pairs.
    contaminated_pairs = build_overlap_list(args.overlap_file)

    # Create a new, cleaned merged file that is devoid of any UMI-TCR pairs found in the overlap list of tuples.
    write_clean_file(args.merged_file, contaminated_pairs, output_file)

    # Record time & date of program completion.
    print("{} program finished on {} for merged file {} and overlap file {}.".format(os.path.basename(__file__),
                                                                                     datetime.now().strftime(
                                                                                         "%A, %B %d at %I:%M:%S %p"),
                                                                                     os.path.basename(args.merged_file),
                                                                                     os.path.basename(
                                                                                         args.overlap_file)))
    print("Time taken: {} seconds.".format(clock() - start_time))


if __name__ == "__main__":
    main()

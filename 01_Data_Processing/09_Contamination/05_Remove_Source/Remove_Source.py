"""
FILENAME:       Remove_Source.py
DESCRIPTION:    This script reads in a directory containing overlap files and a file containing
                "source" UMI-TCR pairs. ("Source" UMI-TCR pairs are those UMI-TCR pairs with the earliest
                PCR time.) It then writes a new, cleaned overlap file that is devoid of any lines
                containing a source UMI-TCR pair for that sample that contains the source UMI-TCR.
INPUT:          A directory containing files with the name format "ol-compact-merged-*.txt". A directory containing
                a file that has a list of UMI-TCR pairs to keep (i.e. the "source" UMI-TCRs).
OUTPUT:         A directory containing files that are devoid of the "source" UMI-TCR pairs (i.e. cleaned overlap files.)
CREATED:        February 8, 2018.
LAST MODIFIED:  April 6, 2018.
MODIFICATIONS:  The program now works on a single file, for use with a generate script ("Generate_Remove_Source.py").
                Added version information.
VERSION:        3.0
AUTHOR(S):      Jay Sharma.
STATUS:         Completed and working.
"""

from __future__ import print_function

import argparse
import os
import re
from collections import defaultdict
from datetime import datetime
from time import clock


def extract_info(f_in):
    """Function to extract the barcode, donor and chain from each sample file name."""
    info = re.search(r"merged-P7M1S-(\w{2,})-(\S*)-(alpha|beta)combined.txt", f_in)  # Create regex capture groups.
    try:
        barcode, donor, chain = info.group(1), info.group(2), info.group(3)
    except AttributeError:  # This should never trigger, but warn the user if there's an errant file name.
        barcode, donor, chain = None, None, None
        print("Error! Regular expression couldn't extract information from file {}.".format(f_in))
    return barcode, donor, chain


def define_outputs(path_in, dir_out):
    """Function to read input filenames and create output files for each sample."""
    base_in = os.path.split(path_in)[1]
    name_out = "cl-{}".format(base_in)
    # name_out = base_in  # Keep the exact name as the input file name.
    path_out = os.path.join(dir_out, name_out)
    write_out = open(path_out, 'w')
    return write_out


def parse_file(line):
    """Function that extracts the TCR, UMI, filename and PCR time from every line."""
    columns = line.rstrip('\n').split('\t')
    tcr, umi, pcr_time, file_name = columns[0], columns[1], columns[5], columns[6]  # Extract the UMI and TCR.
    return tcr, umi, pcr_time, file_name


def build_source_dictionary(f_source):
    """Function to parse the source file (the file containing the earliest PCR time
    for a UMI-TCR pair for a sample) and store every UMI-TCR in a dictionary for a given
    sample. It returns a dictionary with the tuple of (barcode, donor, chain) as the key
    and a list of tuples of (UMI, TCR) as value.
    """
    sources = defaultdict(list)
    with open(f_source) as fh:
        for line in fh:
            if "TCR" in line:  # Skip the line containing column headers, if it exists.
                continue
            tcr, umi, pcr_time, file_name = parse_file(line)  # Extract the TCR, UMI, PCR time & filename.
            key = extract_info(file_name)  # Extract the barcode, donor & chain for a sample.
            value = umi, tcr, pcr_time
            sources[key].append(value)  # Create a list of tuples (UMI, TCR, PCR time) for every sample.
    return sources


def remove_source(dict_source, f_in, f_out):
    """Function to remove the source (earliest) UMI-TCR pair from the overlap file of every sample that
    contained that source UMI-TCR pair. First the tuple of (barcode, donor, chain) is extracted from
    the current overlap file. Then this tuple is looked up in the dictionary of source values to find
    all (UMI, TCR, PCR time) tuples for the current overlap file. Lastly, a new, "cleaned" overlap file
    is written that's devoid of any lines that contain a source tuple (i.e. tuples to be removed).
    """
    with open(f_in) as overlap_file:
        overlap_info = extract_info(f_in)  # Extract the sample information from the overlap file...
        list_removal = dict_source[overlap_info]  # ...& use it to get all UMI, TCR & PCR time for that overlap file.
        for line in overlap_file:
            tcr, umi, pcr_time, _ = parse_file(line)  # "_" is the filename, but it's unused in this function.
            overlap_umi_tcr_pcr = umi, tcr, pcr_time  # Make a tuple of UMI, TCR & PCR time to check against the source tuple.
            if overlap_umi_tcr_pcr not in list_removal:  # The overlap tuple is not in the list of tuples to be removed.
                f_out.write(line)
    f_out.close()


def main():
    cmdline = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    cmdline.add_argument('input_file', help='the unlcean overlap file')
    cmdline.add_argument('output_directory', help='the directory that will contain clean overlap files')
    cmdline.add_argument('source_file', help='the file that contains source UMI-TCR pairs to remove')
    cmdline.add_argument('-v,' '--version', action='version', version='%(prog)s 3.0')
    cmdline.add_argument('-e', '--header', help='write the column headers in the clean overlap file',
                         action='store_true', default=False)
    args = cmdline.parse_args()

    # Record time & date of program start.
    start_time = clock()
    print("{} program started on {} for file {}.".format(os.path.basename(__file__),
                                                         datetime.now().strftime("%A, %B %d at %I:%M:%S %p"),
                                                         os.path.basename(args.input_file)))

    # If the output directory doesn't exist, create it.
    if not os.path.exists(args.output_directory):
        os.makedirs(args.output_directory)

    # Define the column headers.
    header = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("Self: 100 Nucleotide TCR",
                                                                               "Self: 12 Nucleotide UMI",
                                                                               "Self: CDR3 Amino Acid",
                                                                               "Self: UMI Count",
                                                                               "Self: Number of Reads",
                                                                               "Self: PCR Time", "Self: File Name",
                                                                               "Other: 100 Nucleotide TCR",
                                                                               "Other: 12 Nucleotide UMI",
                                                                               "Other: CDR3 Amino Acid",
                                                                               "Other: UMI Count",
                                                                               "Other: Number of Reads",
                                                                               "Other: PCR Time", "Other: File Name")

    # Do computation.
    source_dict = build_source_dictionary(args.source_file)  # Build a dictionary of UMI-TCR pairs.
    output_file = define_outputs(args.input_file, args.output_directory)  # Create the output file for each input file.
    if args.header:  # Write the column headers to the output file if the user specifies the "--header" argument.
        output_file.write(header)
    remove_source(source_dict, args.input_file, output_file)  # Write data to the output file.

    # Record time & date of program completion.
    print("{} program finished on {} for file {}.".format(os.path.basename(__file__),
                                                          datetime.now().strftime("%A, %B %d at %I:%M:%S %p"),
                                                          os.path.basename(args.input_file)))
    print("Time taken: {} seconds.".format(clock() - start_time))


if __name__ == "__main__":
    main()

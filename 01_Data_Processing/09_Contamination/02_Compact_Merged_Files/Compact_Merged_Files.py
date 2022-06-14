"""
FILENAME:       Compact_Merged_Files.py
DESCRIPTION:    This script reads in merged-*.txt files and compacts the UMI to 12 nucleotides and the TCR to 100
                nucleotides. If the UMI is 16 nucleotides, it just collapses. If the UMI is 19 nucleotides,
                it strips the last 3 nucleotides -- ostensibly, the translated PCR tag -- and translates it
                back to the original PCR time tag. The 16 nucleotide UMI is then collapsed. Additionally, each
                CDR3 amino acid sequence (one per line) is checked for any stop codons or trash sequences.
                Finally, the UMI sequence, TCR sequence, CDR3 amino acid sequence, number of UMIs, number of reads
                and earliest PCR time are written on every line for a given TCR.
INPUT:          A merged-*.txt file. A txt file containing the translation tag (e.g. "tag_conversion.txt").
OUTPUT:         A compact-merged-*.txt file.
CREATED:        January 29, 2018.
LAST MODIFIED:  April 6, 2018.
MODIFICATIONS:  Disabled check for proper CDR3 amino acid sequences, which is done by Cut_and_Fix_Merged_Files.py.
                Enabled summation of duplicate lines.
VERSION:        5.0
AUTHOR(S):      Jay Sharma
STATUS:         Completed and working.
"""

from __future__ import print_function

import argparse
import os
from collections import defaultdict, Counter
from datetime import datetime
from itertools import groupby
from operator import itemgetter
from time import clock, sleep


def define_outputs(path_in, dir_out):
    """Function to read input filenames and create output files."""
    base_in = os.path.splitext(os.path.split(path_in)[1])[0]
    name_out = "compact-{0}.txt".format(base_in)
    path_out = os.path.join(dir_out, name_out)
    write_out = open(path_out, 'w')
    return write_out


def check_cdr3(aa):
    """Function to check if the CDR3 amino acid sequence contains any stop codons
    or other trash sequences. It looks for "?" or "*" characters in the string.
    """
    if "?" in aa or "*" in aa:
        return False
    return True


def reverse_translate(trn_tag, trn_dict):
    """Function to convert any translated tags (three nucleotides)
    back to their original PCR time based on the translation dictionary.
    """
    try:
        org_tag = trn_dict[trn_tag]  # Look up trn_tag in the translation dictionary & get the corresponding org_tag.
    except KeyError:
        print(
            "Error! The tag {} wasn't found in the translation table. Make sure the translation table is up to date. "
            "The program will continue with the rest of the file, but this tag will be skipped!".format(trn_tag))
        return None
    return org_tag


def change_umi(umi):
    """Function to change the UMI depending on its length and extract the PCR nucleotide tag if it has one.
    The length should be between sixteen and twenty nucleotides. Lengths of seventeen or eighteen nucleotides
    indicate a fake UMI for an unmodified (tagged) UMI, and twenty nucleotides indicates a fake UMI for
    an unmodified (untagged) UMI. No matter the length, the UMI is ultimately collapsed to twelve nucleotides.
    """
    umi_length = len(umi)
    if 16 == umi_length:  # Unmodified (untagged) UMI; use it as is, with no tag.
        original_umi, translated_tag = umi, None
    else:  # Modified (tagged) UMI; separate it from the tag.
        original_umi, translated_tag = umi[0:16], umi[16:19]
    collapsed_umi = original_umi[1:5] + original_umi[6:10] + original_umi[11:15]  # Collapse the UMI to 12 nucleotides.
    return original_umi, collapsed_umi, translated_tag


def parse_line(line_in, trn_dict):
    """Function to process each line of the merged file. For the UMI, any nucleotides past the 19th will be stripped.
    For the TCR, any nucleotides past the 100th will be stripped. The CDR3 amino acid will also be checked for
    stop codons or garbage sequences.
    """
    split = line_in.rstrip('\n').split('\t')
    try:
        aa, modified_umi, tcr, reads = split[-4], split[-3][:19], split[-2][:100], int(split[-1])  # Extract & collapse.
        original_umi, collapsed_umi, translated_tag = change_umi(modified_umi)
        if translated_tag is not None:  # The UMI is 19 nucleotides long.
            pcr_time = reverse_translate(translated_tag, trn_dict)  # Get the original PCR tag back.
            if pcr_time is not None:  # The translation succeeded.
                return tcr, collapsed_umi, aa, reads, pcr_time
            else:  # The translation failed; the last 3 nucleotides aren't in the translation table.
                return tcr, collapsed_umi, aa, reads, "Incorrect tag"
        else:  # The UMI is 16 nucleotides long.
            return tcr, collapsed_umi, aa, reads, "No tag"
    except IndexError:
        return None, None, None, None, None


def count_umis(list_in):
    """Function that reads in every UMI, CDR3 AA, number of reads & PCR time (a list of lists
    called "list_in"; a sublist) for a TCR and calculates how many UMIs exist for that TCR. This is
    done by first using collections.Counter(), which creates a new dictionary of UMI as the key and
    count as the value. The last step is to iterate through through each sublist and insert
    the respective UMI count in that list.

    Note that this function can determine the number of _individual_ UMIs, or the _total_ number of UMIs
    -- the latter via the len(umi_counter) function -- depending on what the user wants.

    For details on the methodology, refer to [https://stackoverflow.com/a/48629163].
    """
    umi_counter = Counter(sublist[0] for sublist in list_in)  # Count the number of UMIs; sublist[0] is the UMI.
    # umi_count = len(umi_counter)  # Get the _total_ number of unique UMIs (if desired).
    for sublist in list_in:  # Iterate through each sublist, look up the UMI count, & insert it into the sublist.
        sublist.insert(2, umi_counter[sublist[0]])
    return list_in  # "list_in" now has the UMI count inserted (before the reads) in it.


def sum_reads(list_in):
    """Function that reads in every UMI, CDR3 AA, number of UMIs, number of reads & PCR time
    (a list of lists called "list_in"; a sublist) for a TCR and sums the number of reads that have
    the same UMI, CDR3 AA, UMI count & PCR time. This means they're duplicates from the merged file
    and need to be collapsed. First each sublist is grouped by UMI, CDR3 AA, UMI count & PCR time.
    Then the reads of these grouped sublists are summed together; if there are no
    duplicates, nothing will be summed, and the output list will be the same as the input list.

    For details on the methodology, refer to [https://stackoverflow.com/a/48530763].
    """
    list_in.sort(key=lambda x: (x[0], x[1], x[2], x[-1]))  # First sort the list by UMI, CDR3 AA, UMI count & PCR time.
    grouped_list = [(umi_cdr3_pcr, list(group)) for umi_cdr3_pcr, group in
                    groupby(list_in, key=lambda x: (x[0], x[1], x[2], x[-1]))]
    list_out = [[collapsed_umi, cdr3_aa, umi_count, sum(i[3] for i in values_list), pcr_time] for
                (collapsed_umi, cdr3_aa, umi_count, pcr_time), values_list in grouped_list]
    return list_out


def get_earliest(list_in):
    """Function that reads in every UMI, CDR3 AA, number of UMIs, number of reads & PCR time (a list
    of lists called "list_in"; a sublist) for a TCR and determines the earliest PCR time. It does so by
    first grouping each "list_in" sublist by PCR time. Then the min() function is called to obtain the
    sublist with the smallest PCR time. If no sublist has a smallest PCR time, they're all returned.
    """
    list_in.sort(key=lambda x: (x[-1]))  # First sort the list by earliest PCR time.
    grouped_list = [(pcr_time, list(group)) for pcr_time, group in groupby(list_in, key=lambda x: (x[-1]))]
    list_out = min(grouped_list, key=itemgetter(0))[1]  # itemgetter gets every "list_in" with the smallest PCR time.
    return list_out


def finalize_values(v, keep_earliest):
    """Function that reads in every UMI, CDR3 AA, number of reads & PCR time (the "value",
    which is a list) for a given TCR (the "key") and determines the unique UMI count and
    smallest PCR time. This function handles all three cases:
    1. Different PCR times for every TCR --> Return the value with the earliest PCR time.
    2. Same PCR times but different UMIs and CDR3 AAs for every TCR --> Return all values.
    3. Same PCR times, UMIs and CDR3 AAs for every TCR --> Sum the reads of these same values and return all values.

    The "keep_earliest" argument is a boolean that's specified by the user as "--earliest". If it's false
    (or the argument isn't specified), then ALL sequences will be compacted and output. If it's true, then
    just the sequences that have the earliest PCR time will be output.

    If there's no single earliest PCR time, the values are checked for duplicate UMIs and CDR3 AAs.
    If any duplicates are found, their reads are summed. If no duplicates are found, and no PCR time is determined
    to be the earliest, all values are returned. In any case, the value returned is a list.
    """
    umi_list = count_umis(v)  # Obtain the number of UMIs.
    if keep_earliest:
        smallest_list = get_earliest(umi_list)  # Get only the sequence(s) with the earliest PCR time.
    else:
        smallest_list = umi_list  # Don't get the sequence with the earliest PCR time; output ALL sequences.
    if len(smallest_list) > 1:  # If this condition is met, there may be duplicates...
        smallest_list = sum_reads(smallest_list)  # ...so sum their reads.'''
    return smallest_list  # This list has the UMI count inserted, any duplicates summed & only the earliest PCR time.


def convert_merged(m_dict, f_out, f_in, is_earliest):
    """Main function to perform the compacting of the merged file. The merged file dictionary is
    iterated over, and each value -- which is a list containing UMI, CDR3 AA, reads & PCR time --
    moves on to final processing, which can include determining earliest PCR time, UMI count,
    and summation of reads for duplicate values.
    """
    for k, v in m_dict.iteritems():  # Process the dictionary (with TCR as the key) for writing.
        final_list = finalize_values(v, is_earliest)  # Perform grouping & final processing (e.g. count UMIs).
        for i in final_list:  # i is a list containing: UMI, CDR3 AA, UMI count, reads & earliest PCR time.
            collapsed_umi, cdr3_aa, number_umis, reads, pcr_time = i[0], i[1], i[2], i[3], i[4]
            line_out = '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(k, collapsed_umi, cdr3_aa, number_umis, reads, pcr_time,
                                                             os.path.basename(f_in))
            f_out.write(line_out)  # Write the compact-merged-*.txt file.
    f_out.close()


def build_merged_dictionary(f_in, trn_dict):
    """Function to read in the merged file and create a dictionary of lists. The TCR sequence is the
    key and a list containing 12-nucleotide UMI, CDR3 AA, number of reads and PCR time is the value.
    Note that the value must be a list in this case, since later on the UMI count will be inserted
    into it, so a tuple won't work (since tuples are immutable).
    """
    tcr_dict = defaultdict(list)  # Create a defaultdict to hold info for every TCR pair.
    with open(f_in) as merged_file:  # Read in the merged file & process it.
        for line_in in merged_file:
            tcr, collapsed_umi, cdr3_aa, reads, pcr_time = parse_line(line_in, trn_dict)
            if cdr3_aa is not None:  # The input line passes the CDR3 AA check.
                value = [collapsed_umi, cdr3_aa, reads, pcr_time]
                tcr_dict[tcr].append(value)  # Build the dictionary.
    return tcr_dict


def build_translation_dictionary(f_in):
    """Function to read a tab-delimited tag conversion file (e.g. "tag_conversion.txt") and
    create a dictionary with three-nucleotide tag as the key and original PCR time as the value.
    """
    trn_dict = {}  # Make a new dictionary to create a mapping from old tag to new tag.
    with open(f_in) as trn_tab:
        for line in trn_tab:
            if "Tag" in line or "tag" in line:
                continue  # Skip the line containing the column headers, if it exists.
            else:
                columns = line.rstrip('\n').split('\t')
                old_tag, new_tag = columns[0], columns[1]  # Extract the PCR time & nucleotide tag.
                trn_dict[new_tag] = old_tag  # key = nucleotide tag; value = PCR time
    return trn_dict


def main():
    cmdline = argparse.ArgumentParser(prog='Compact_Merged_Files.py', description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    cmdline.add_argument('input_file', help='the merged-*.txt file to be compacted')
    cmdline.add_argument('output_directory', help='the directory that will contain the compact-merged-*.txt files')
    cmdline.add_argument('path_to_table',
                         help='the full path to the conversion table which is used to translate old tags to new tags')
    cmdline.add_argument('-v,' '--version', action='version', version='%(prog)s 5.0')
    cmdline.add_argument('-a', '--earliest', help='keep only the sequences with the earliest PCR time',
                         action='store_true', default=False)
    cmdline.add_argument('-e', '--header', help='write the column headers in the output file', action='store_true',
                         default=False)
    args = cmdline.parse_args()

    # Record time & date of program start.
    start_time = clock()
    print("{} program started on {} for file {}.".format(os.path.basename(__file__),
                                                         datetime.now().strftime("%A, %B %d at %I:%M:%S %p"),
                                                         os.path.basename(args.input_file)))

    # If the output directory doesn't exist, create it. Add a 1-second wait period to account for race conditions.
    if not os.path.exists(args.output_directory):
        try:
            os.makedirs(args.output_directory)
        except OSError, e:
            if e.errno != os.errno.EEXIST:
                raise
            sleep(1)
            pass

    # Create the output file from the input file.
    output_file = define_outputs(args.input_file, args.output_directory)

    # Define the column headers & write it to the output file if the user specifies the "--header" argument.
    header = "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("100 Nucleotide TCR", "12 Nucleotide UMI", "CDR3 Amino Acid",
                                                   "UMI Count", "Number of Reads", "Earliest PCR Time", "File Name")
    if args.header:
        output_file.write(header)

    # Read in the translation table file & convert it to a dictionary for translation.
    translation_table = build_translation_dictionary(args.path_to_table)

    # Read in the merged file & build a dictionary (with TCR as the key) out of it.
    merged_dictionary = build_merged_dictionary(args.input_file, translation_table)

    # Perform the compacting of the merged file.
    convert_merged(merged_dictionary, output_file, args.input_file,
                   args.earliest)  # "args.input_file" is only used to extract the filename.

    # Record time & date of program completion.
    print("{} program finished on {} for file {}.".format(os.path.basename(__file__),
                                                          datetime.now().strftime("%A, %B %d at %I:%M:%S %p"),
                                                          os.path.basename(args.input_file)))
    print("Time taken: {} seconds.".format(clock() - start_time))


if __name__ == "__main__":
    main()

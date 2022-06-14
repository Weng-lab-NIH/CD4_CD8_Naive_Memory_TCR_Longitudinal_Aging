"""
FILE NAME:      generate_umi_grouping_part2.py

COMPLEMENTARY
SCRIPT(S)/
FILE(S):        umi_grouping_part2.py; submit_umi_regrouping2.sh

COMPATIBLE
PIPELINE(S):    Human

DESCRIPTION:    FASTQ files have the following stucture:

                1) @D00787:87:H7HJVBCXY:1:1107:7107:2097 UMI:TTCCGTGATGTAAAATAAA:H@CF<<<1<C<1<1DG
                2) GAGCANTGTCAGCCGGGTGCCTGGGCCAAAATACTGCGTATCTCGATACCCTGGCAAGCTGCTGGCACAGAGATACAT...
                +
                3) @B@@D#<<<<DGH@H?DHHHHIFHIHHICHHIHIIIIEHCGIIGIHHEECCHFHHH@CCHIHIHHIEHHIIIIHFHGH...

                Where lines 1, 2, and 3 represent the Read header, T-cell receptor (TCR)
                sequence, and Quality score respectively. In the Read Header, there exists a 16
                base-pair Unique Molecular Identifier (UMI) that is added onto the TCR sequence
                prior to sequencing. For contamination analysis, we need to be able to link
                each TCR sequence to when it was PCR amplified. By doing so, if the same
                UMI-TCR sequence pair are found across multiple samples, we can identify when
                the UMI-TCR sequence pair was first sequenced (the "source"). Once we isolate
                the source, any other instances of the same UMI-TCR sequence pair are as a
                result of contamination and we can eliminate them from other samples.In order
                to identify the source, each UMI is affixed an additional 3 base-pair
                nucleotide sequence to indicate when that sequence was amplified (a "PCR
                Time/Tag") where a smaller PCR time (e.g.: a PCR time of 5) indicates that
                specific sequence was sequenced earlier than another sequence with a larger PCR
                time (e.g.: a PCR time of 8).

                Once we cluster similar TCR sequences together, we need to
                identify the earliest PCR associated with the cluster. After doing so, we
                generate a new UMI for that specific sequence by using the old UMI and affixing
                it with the 3bp sequence corresponding to the earliest PCR time for PCR source
                identification purposes. If this new UMI has already been used, then we add a
                random base-pair onto the UMI (it now becomes a 20bp or more nucleotide
                sequence or 19bp nucleotide sequence). Each sequence and its corresponding
                quality score gets re-written to a new FASTQ file with its newly associated
                UMI.

INPUT
FILES(S):       1) Summary File for Combined Files
                2) PCR conversion table

OUTPUT
FILES(S):       1) A FASTQ file with new UMI's based on clustered sequences

INPUT
ARGUMENT(S):    1) Input Directory: The directory in which the DivE Summary Files are located in.
                2) Output Directory: The directory where the files
                3) PCR Conversion Table: The direct path to the PCR Conversion Table that has the PCR times
                    associated with each 3bp nucleotide
                4) Input Compression: Specify whether the combined files are compressed or not
                5) Output Compression: Specify whether the output file will be compressed or not
                6) Verbosity (optional): Used primarily for debugging purpses (will list TCR & UMI sequence, cluster
                number, etc.)

CREATED:        07Mar18

MODIFICATIONS
LOG:            09Mar18     Changed structure & reformatted code for faster processing.
                22Mar18     Split original script into 2 parts

LAST MODIFIED
BY:             Thomas Nguyen

PYTHON
VERSION USED
TO WRITE
SCRIPT:         2.7.14

VERSION:        1.1

AUTHOR(S):      Thomas Nguyen

STATUS:         Working

TO DO LIST:     None at this moment

DISCLAIMER(S):  The gzip python module is poorly optimized to create large compressed files so run times for
                writing compressed files (especially large ones) will take longer than usual; if need be, write
                the output file as uncompressed but compress the file after processing is done.
"""

from __future__ import print_function

import  os, time, argparse, logging,sys,gzip
from itertools import chain, product
from collections import defaultdict, namedtuple

Read = namedtuple("Read", ['original_umi', 'time_tag', 'id', 'umi_qual', 'seq', 'qual'])


def read_grouping_part1(file,compression):
    """
    This function reads in the summary file and extracts the information contained in each tab-delineated column
    and stores it into a defaultdict.
    """
    cluster_info = defaultdict(list)
    reads_count = 0
    if compression == 'N':
        with open(file) as input:
            for line in input:
                if "Original UMI" in line:
                    continue
                else:
                    line_info = line.rstrip().split('\t')
                    original_umi, time_tag, id = line_info[0], line_info[1], line_info[2]
                    umi_qual, seq, qual, cluster_n = line_info[3],line_info[4], line_info[5], line_info[6]

                    r = Read(original_umi, time_tag, id, umi_qual, seq, qual)
                    cluster_info[cluster_n].append(r)
                    reads_count += 1
    elif compression == 'Y':
        with gzip.open(file,'rb') as input:
            for line in input:
                if "Original UMI" in line:
                    continue
                else:
                    line_info = line.rstrip().split('\t')
                    original_umi, time_tag, id = line_info[0], line_info[1], line_info[2]
                    umi_qual, seq, qual, cluster_n = line_info[3],line_info[4], line_info[5], line_info[6]

                    r = Read(original_umi, time_tag, id, umi_qual, seq, qual)
                    cluster_info[cluster_n].append(r)
                    reads_count += 1

    logging.info("Finished Reading in %d clusters and %d reads",
                 len(cluster_info),
                 reads_count)
    return cluster_info


def build_translation_dictionary(f_in):
    """
    Function to read in a tab-delimited tag conversion file (e.g.
    "tag_conversion.txt") and create a dictionary with three-nucleotide PCR
    time as the key and original PCR time as the value.
    """
    trn_dict = {}
    with open(f_in) as trn_tab:
        for line in trn_tab:
            if "Tag" in line or "tag" in line:
                continue
            else:
                fields = line.strip().split('\t')
                # old tag: PCR time
                # new tag: nucleotide sequence representing PCR time
                old_tag, new_tag = int(fields[0][1:]), fields[1]
                trn_dict[new_tag] = old_tag
    return trn_dict


def find_earliest_pcr_time(translation_dictionary, new_tag_list):
    """
    This function reads in the translation tag-PCR dictionary and a list of the
    translation tags associated with the TCR sequences in each cluster. From
    this, it calculates the earliest PCR time and uses that.
    """

    pcr_time_list = []
    for new_tag in new_tag_list:
        # WR: idiomatic python doesn't need the .keys()
        # WR: if new_tag in translation_dictionary.keys():
        if new_tag in translation_dictionary:
            # WR: add nucleotide tag along with integer time
            pcr_time_list.append((translation_dictionary[new_tag], new_tag))
    # WR: lowest_pcr_time = min(int(new_tag) for new_tag in pcr_time_list)
    # WR: pcr times are already integers and are called 'old_tag' in previous function
    pcr_time_list.sort(key=lambda a: a[0])
    lowest_pcr_time_tag = pcr_time_list[0][1]
    return lowest_pcr_time_tag


#WR: I think build_new_umi should be deterministic and should prefer single
#    nt extension. Let's build a fixed extension list of 1, 2, and 3 nt
#    combinations. That's 84 possible extensions. That should be enough
umi_extension_tags = [''.join(a) for a in
        chain(product('ACGT', repeat=1), product('ACGT', repeat=2), product('ACGT', repeat=3),
              product('ACGT', repeat=4),product('ACGT', repeat=5),product('ACGT', repeat=6),
              product('ACGT', repeat=7),product('ACGT', repeat=8))]

umi_extension_tags = sorted(umi_extension_tags,key=len)


def build_new_umi(current_umi, used_umi_set,umi_extension_tags):
    """
    This function builds a new UMI if the current UMI already exists in a set.
    It builds the UMI by adding a random nucleotide to the existing UMI. The
    number of nucleotides added depends on how many variations of the UMI are
    used (UMI's can be greater than 20 bp in some cases).
    """
    for extension in umi_extension_tags:
        possible_umi = current_umi + extension
        if possible_umi not in used_umi_set:
            return possible_umi
    # this should never happen. If it does, may need to expand
    # umi_extension_tags list
    raise ValueError("could not create a unique extended UMI")


def write_regrouped_file(output_file,read_info,new_umi):
    """
    This function uses the old and new FASTQ information with (FASTQ ID, UMI Sequence, UMI Quality score, TCR
    sequence, Quality Score) to rewrite the FASTQ file.
    """
    output_header = "@{} UMI:{}:{}".format(read_info.id,new_umi,read_info.umi_qual)
    output_seq = read_info.seq
    output_qual = read_info.qual
    output_read = "{}\n{}\n+\n{}\n".format(output_header,output_seq,output_qual)
    output_file.write(output_read)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("input_summary",
                        help="The input summary file")
    parser.add_argument("output_directory",
                        help="The directory where the regrouped FASTQ files will be located in")
    parser.add_argument("PCR_Conversion_Table",
                        help="The absolute file path to the PCR conversion table")
    parser.add_argument("input_compression",
                        choices = {'Y','N'},
                        help = "Specify whether the input FASTQ files are compressed or not")
    parser.add_argument("output_compression",
                        choices = {'Y','N'},
                        help = "Specify whether the output summary files should be compressed or not")
    parser.add_argument('--verbose', '-v', action='store_true', default=False)

    args = parser.parse_args()

    translation_tags = build_translation_dictionary(args.PCR_Conversion_Table)

    # create output directory if it doesn't exist. account for concurrent processes
    if not os.path.exists(args.output_directory):
        try:
            os.makedirs(args.output_directory)
        except OSError, e:
            if e.errno != os.errno.EEXIST:
                raise
        time.sleep(1)

    outfile_prefix = os.path.splitext(os.path.basename(args.input_summary))[0]
    outfile_prefix = outfile_prefix.replace('summary_file-','')

    if args.output_compression == 'N':
        outfile_name = os.path.join(args.output_directory, outfile_prefix + '_regrouped.fastq')
        outfile = open(outfile_name,'w')
    elif args.output_compression == 'Y':
        outfile_name = os.path.join(args.output_directory, outfile_prefix + '_regrouped.fastq.gz')
        outfile = gzip.open(outfile_name,'wb')

    # WR: logging setup
    if args.verbose:
        log_level = logging.DEBUG
    else:
        log_level = logging.INFO
    logging.basicConfig(level=log_level,
                        format='%(asctime)s|%(filename)s|%(levelname)s|%(message)s',
                        datefmt="%A, %B %d at %I:%M:%S %p")

    start_time = time.clock()
    logging.info("started for file %s", os.path.basename(args.input_summary))

    clusters_dict = read_grouping_part1(args.input_summary,args.input_compression)

    if args.verbose:
        print ("Cluster Number\tUMI\tTime Tag\tSequence",file=sys.stderr)
        for cluster_n,cluster in clusters_dict.iteritems():
            for read in cluster:
                print ("%3s\t%s\t%3s\t%s" % (cluster_n,read.original_umi,
                                               read.time_tag,read.seq),file=sys.stderr)



    used_umi_set = set()

    for cluster_n,cluster_reads in clusters_dict.iteritems():
        umi_list = [r.original_umi for r in cluster_reads]
        umi_list = set(umi_list)
        if len(umi_list) == 1:  #check to see if each cluster has 1 UMI associated with it
            umi = umi_list.pop()
        else:
            raise ValueError("There are multiple original UMI's in this cluster")

        all_pcr_tags = [r.time_tag for r in cluster_reads]
        earliest_pcr_time = find_earliest_pcr_time(translation_tags,all_pcr_tags)

        new_translated_umi = umi + earliest_pcr_time
        if new_translated_umi in used_umi_set:
            new_translated_umi = build_new_umi(new_translated_umi,used_umi_set,umi_extension_tags)
        used_umi_set.add(new_translated_umi)

        for indiv_read in cluster_reads:
            write_regrouped_file(outfile,indiv_read,new_translated_umi)


    outfile.close()

    finish_time = time.clock()
    logging.info("Creating New UMI's and rewriting FASTQ file is finished.")
    logging.info(" Processing time for Rewriting file was %d sec",
                 finish_time - start_time)


if __name__ == "__main__":
    main()

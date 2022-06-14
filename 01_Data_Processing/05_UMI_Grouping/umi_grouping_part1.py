"""
FILE NAME:      umi_grouping.py

COMPLEMENTARY
SCRIPT(S)/
FILE(S):        generate_umi_grouping_part1.py; submit_umi_regrouping.sh


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

                With our current FASTQ files, there are several UMI's that have over 300,000
                associated TCR sequences with them. We would like to identify how many
                "clusters" of similar TCR sequences there are within a 16 bp UMI because each
                "cluster" may encode for a different TCR. Within 1 16bp UMI, we identify
                similar TCR sequences by calculating the Hamming Distance between each sequence
                and cluster them based on similarity (the current cutoff criteria for Hamming
                Distance is 30). Once the sequences have been clustered together based on
                similiarity, the remaining FASTQ file information (FASTQ ID, UMI, PCR Time Tag,
                UMI Quality Score, Sequence, Sequence Quality Score) is written to a tab-delineated
                .txt file along with the "cluster number".

INPUT
FILES(S):       1) Combined Files

OUTPUT
FILES(S):       1) Tab Delinated .txt file containing TCR sequence, UMI sequence, PCR Time tag, FASTQ ID header,
                UMI Quality Score, Sequence Quality Score, and Cluster Number

INPUT
ARGUMENT(S):    1) Input FASTQ File: The combined FASTQ file that will be processed by this script
                2) Output Directory: The directory where the files
                3) Input Compression: Specify whether the combined files are compressed or not
                4) Output Compression: Specify whether the output file will be compressed or not
                5) Verbosity (optional): Used primarily for debugging purpses (will list TCR & UMI sequence, cluster
                number, etc.)

CREATED:        07Mar18

MODIFICATIONS
LOG:            09Mar18     Changed structure & reformatted code for faster processing.
                22Mar18     Split original script into 2 parts

LAST MODIFIED
BY:             Wolfgang Resch; Thomas Nguyen

PYTHON
VERSION USED
TO WRITE
SCRIPT:         2.7.14

VERSION:        1.2

AUTHOR(S):      Thomas Nguyen

STATUS:         Working

TO DO LIST:     Look into multithreading/multiprocessing for faster analysis.

DISCLAIMER(S):  The gzip python module is poorly optimized to create large compressed files so run times for
                writing compressed files (especially large ones) will take longer than usual; if need be, write
                the output file as uncompressed but compress the file after processing is done.
"""

from __future__ import print_function

import sys, os, time, argparse, logging,Levenshtein,gzip
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from collections import defaultdict, namedtuple, Counter
import operator


Read = namedtuple("Read", ['original_umi', 'time_tag', 'id', 'umi_qual', 'seq', 'qual'])


def read_fastq(filename,compression):
    """
    This function reads the FASTQ and creates 3 dictionaries: 1 dictionary
    that contains the original 16bp UMI as the key, and the translation tag and
    TCR sequence as the values inside a tuple; another dictionary that contains
    the original 16bp UMI as the key and the TCR sequences associated with it
    as the value; a dictionary with the sequence as the key and the FASTQ Read
    header/quality score as the values.
    """
    #WR: doc string does not match code. in the code all three dics are keyed by
    #    original_umi.
    #WR: for a large file this will be a lot of memory - by creating three dictionaries
    #    you are nearly tripling the amount of memory required. I replaced that with
    #    a single original_umi keyed dictionary holding namedtuples with all information
    #    necessary

    reads = defaultdict(list)

    nreads = 0
    if compression == 'N':
        with open(filename) as input:
            for title, seq, qual in FastqGeneralIterator(input):
                nreads += 1
                #WR: splitting the UMI string and header here since
                #    it needs to be reassembled again for the output with the new UMI
                seqid, umistr = title.strip().split(' ')

                _, translated_umi, umi_qual = umistr[:3],umistr[4:23],umistr[24:]
                original_umi = translated_umi[:16]
                time_tag = translated_umi[16:]
                if len(seq) == 150:
                    r = Read(original_umi, time_tag, seqid, umi_qual, seq, qual)
                    reads[original_umi].append(r)
    else:
        with gzip.open(filename,'rb') as input:
            for title, seq, qual in FastqGeneralIterator(input):
                nreads += 1
                #WR: splitting the UMI string and header here since
                #    it needs to be reassembled again for the output with the new UMI
                seqid, umistr = title.strip().split(' ')

                _, translated_umi, umi_qual = umistr[:3],umistr[4:23],umistr[24:]
                original_umi = translated_umi[:16]
                time_tag = translated_umi[16:]
                if len(seq) == 150:
                    r = Read(original_umi, time_tag, seqid, umi_qual, seq, qual)
                    reads[original_umi].append(r)
    logging.info("read %d reads in %d UMI groups", nreads, len(reads))
    return reads


def cluster_sequences(sequence_list):
    """
    This function clusters sequneces in sequence_list such that all pairwise
    distances within a cluster are below a given cutoff.
	"""

    # WR: not really sure if this is the most efficient way of doing this
    # WR: using the pairwise identity function from biopython instead
    #    of a home made hamming function

    # there are likely to be identical sequences. Each of those should only have
    # to be clustered once. So let's collapse all the identical sequences
    seqd = defaultdict(list)
    for s in sequence_list:
        seqd[s.seq].append(s)

    sequences_to_cluster = seqd.keys()
    # the first cluster
    clusters = [set([sequences_to_cluster.pop()])]
    while sequences_to_cluster:
        new_seq = sequences_to_cluster.pop()
        in_any_cluster = False
        for cluster in clusters:
            in_this_cluster = True
            for cluster_seq in cluster:
                if Levenshtein.hamming(new_seq, cluster_seq) > 29:
                    in_this_cluster = False
                    break
            if in_this_cluster:
                cluster.add(new_seq)
                in_any_cluster = True
                break
        if not in_any_cluster:
            clusters.append(set([new_seq]))

    # expand the identical sequences again
    expanded_clusters = []
    for cluster in clusters:
        cl = []
        for seq in cluster:
            cl.extend(seqd[seq])
        expanded_clusters.append(cl)
    return expanded_clusters


def main():
    parser = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--verbose', '-v',
                        help = "specify whether you want the Logging Function to display extra information ( \
                               TCR Sequence, UMI Sequence, PCR Time Tag, Cluster Number)",
                        action='store_true', default=False)
    parser.add_argument("Input_FASTQ",
            help="The Input FASTQ file")
    parser.add_argument("output_directory",
            help="The directory where the summary files will be located in")
    parser.add_argument("input_compression",
                        choices = {'Y','N'},
                        help = "Specify whether the input FASTQ files are compressed or not")
    parser.add_argument("output_compression",
                        choices = {'Y','N'},
                        help = "Specify whether the output summary files should be compressed or not")

    args = parser.parse_args()

    # create output directory if it doesn't exist. account for concurrent processes
    if not os.path.exists(args.output_directory):
        try:
            os.makedirs(args.output_directory)
        except OSError, e:
            if e.errno != os.errno.EEXIST:
                raise
        time.sleep(1)

    #WR: logging setup
    if args.verbose:
        log_level = logging.DEBUG
    else:
        log_level = logging.INFO
    logging.basicConfig(level=log_level,
            format='%(asctime)s|%(filename)s|%(levelname)s|%(message)s',
            datefmt="%A, %B %d at %I:%M:%S %p")

    start_time = time.clock()
    logging.info("started for file %s", os.path.basename(args.Input_FASTQ))

    summary_file_suffix = os.path.basename(args.Input_FASTQ)
    summary_file_suffix = os.path.splitext(summary_file_suffix)[0]

    if args.output_compression == 'N':
        summary_file_name = 'summary_file-{}.txt'.format(summary_file_suffix)
        summary_file_name = os.path.join(args.output_directory, summary_file_name)
        summary_file = open(summary_file_name, 'w')
    else:
        summary_file_name = 'summary_file-{}.txt.gz'.format(summary_file_suffix)
        summary_file_name = os.path.join(args.output_directory, summary_file_name)
        summary_file = gzip.open(summary_file_name, 'wb')

    header = 'Original UMI\tTime Tag\tFASTQ ID\tUMI Quality Score\tTCR Sequence\tSequence Quality Score\tCluster Number\n'
    summary_file.write(header)


    readdict = read_fastq(args.Input_FASTQ,args.input_compression)

    build_time = time.clock()
    logging.info("finished reading input data [%d sec]",
            build_time - start_time)


    cluster_n = 0   #number of clusters that remain after filtering
    cluster_tot = 0 #total number of clusters (filtering not accounted for)

    unique_umis = len(readdict.keys())
    read_count = 0  #Reads at the end of filtering
    recovered_umi = 0   #UMI's that don't split as a result of clustering

    umi_cluster_count = defaultdict(int)
    umi_tcr_count = defaultdict(int)


    for umi, reads in readdict.iteritems():
        umi_tcr_count[umi] = len(reads)
        clusters = cluster_sequences(reads)

        logging.debug("UMI:%s -- %d reads in %d clusters", umi, len(reads), len(clusters))
        if args.verbose:
            for cluster_num, cluster_reads in enumerate(clusters):
                for r in cluster_reads:
                    print("%3i %5s %s" % (cluster_num, r.time_tag, r.seq), file=sys.stderr)

        if len(clusters) == 1:
            recovered_umi +=1
        for cluster_reads in clusters:
            cluster_tot += 1
            if len(cluster_reads) > 1:
                cluster_n += 1
                umi_cluster_count[umi]=len(cluster_reads)
                for read in cluster_reads:
                    read_count += 1
                    summary_line = '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(read.original_umi,
                                                                 read.time_tag,
                                                                 read.id,
                                                                 read.umi_qual,
                                                                 read.seq,
                                                                 read.qual,
                                                                 cluster_n)
                    summary_file.write(summary_line)
            else:
                if len(clusters) ==1:
                    recovered_umi -=1


    umi_freq_count = Counter([value for key,value in umi_tcr_count.iteritems() if value == 1])
    single_read_umi = umi_freq_count[1]
    max_cluster_umi_num = max(umi_cluster_count.iteritems(),key=operator.itemgetter(1))

    summary_file.close()
    finish_time = time.clock()

    logging.info("Finished Clustering %s UMI's into %s different clusters.",
                 unique_umis,
                 cluster_tot)
    logging.info("Clusters with reads less than 2 were removed from the file.")
    logging.info("%d reads were clustered into %d different clusters.",
                 read_count,
                 cluster_n)
    logging.info("%d of the original UMI's were not split and clustered",
                 recovered_umi)
    logging.info("The number of original UMI's with 1 read is: %d",
                 single_read_umi)
    logging.info("Highest clustering UMI: %s with %d different clusters",
                 max_cluster_umi_num[0],
                 max_cluster_umi_num[1])
    logging.info(" Processing time for clustering was %d sec",
                 finish_time - build_time)




if __name__ == "__main__":
    main()
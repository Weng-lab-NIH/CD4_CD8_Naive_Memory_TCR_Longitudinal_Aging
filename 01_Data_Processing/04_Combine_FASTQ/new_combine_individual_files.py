"""
FILE NAME:      new_combine_individual_files.py

COMPLEMENTARY
SCRIPT(S)/
FILE(S):        generate_new_combine_individual_files.py; submit_new_combine.sh

DESCRIPTION:    Donors are sequenced multiple times and have multiple files and thus we would like to compile
                all of their files into 1 file. Additionally, a particular combined file for a donor may not
                be the most updatedsince that donor may have been sequenced again. This file can combine
                individidual files together into 1 combined file or can append an individual file onto
                a combined file.v

COMPATIBLE
PIPELINES:		Human, Mouse

INPUT(S):       1) Input Directory that contains all of the FASTQ files that needs to be combined.

OUTPUT(S):      1) A FASTQ file with all of the reads compiled for each file of the same donor.

INPUT
ARGUMENT(S):    1) Input Directory (-ip): The directory where the individual/combined FASTQ files reside in.
                2) Output Directory (-op): The directory where the combined FASTQ files will reside in.
                3) Summary File Directory (-sum): The directory where the summary file (number of reads written
                to each combined file) will reside in.
                4) File Type (-ft): The type of FASTQ file (individual, combined, both[individual and combined])
                that needs to be combined.
                5) File Type 2(-ft2): Specify (naive, memory, reproducibility, mouse, neither) what type of files
                are being combined (primarily used for the summary file name purposes).
                6) Input Compression (-ic; optional): Specify whether ALL of the input FASTQ files are
                compressed or not.
                7) Second Input Directory (-ip2; optional): The second input directory that contains FASTQ files
                that needs to be combined with the first input directory FASTQ files. Specify this option if the
                File Type (-ft) argument is "both". If so, the first Input Directory argument (-ip) should contain
                individual files whereas specify the combined FASTQ file directory with this argument.


CREATED:        11Jan18

MODIFICATIONS
LOG: 
           
02Feb18			1) Formatted the script to combine individual files with combined files
03Mar18			2) Added code to handle compressed output files.
03Mar18			3) Changed wording of code to handle compressed input files.
03Mar18			4) Removed ability to compress output files due to slow running time.
11Apr18         5) Added list comprehension syntax to choose barcodes for combined files that don't have
                "a" or "b" as a suffix (e.g.: 24a or 11b).
11Apr18         6) Due to the many arguments required for this script, the optional argument tags
                (-arg,--argument,etc.) have been added to each of the required arguments. Bringing
                up the -h option with python will indicate which of these arguments are required
                and which are.
25Jun18         7) Added additional code to choose most common barcode for combining individual files
                together. Added additional code to choose barcode of combined file when combining
                combined file with individual file.
02Aug19         8) Changed the args.file_type_2 argument to a user-specified argument instead of a
                preset argument.
02Aug19         9) Added extra details to regular expression search of files to handle tagged PCR files.

VERSION:        7.0

AUTHOR(S):      Thomas Nguyen

STATUS:         Working. 

DISCLAIMER(S):	Creating uncompressed combined FASTQ files is quicker than compressed combined FASTQ files.
				The gzip module is poorly optimized for writing large amounts of data.

"""

import argparse, re, os, glob, gzip, time, datetime
from collections import defaultdict, namedtuple, Counter
from Bio.SeqIO.QualityIO import FastqGeneralIterator

File_Info = namedtuple("file_info_list",['file','barcode','batch','file_type','compression_type'])

def getinfo(file,file_list,input_compression):
    '''
    This function obtains the barcode,batch number, and donor of the file and determines what type of files
    (compressed/not compressed; individual/combined). This information is then stored in a tuple and subsequently
    in a defaultdict.
    '''
    if input_compression is False:
        match_1 = re.search('g-P7M1S-(.*)_(.*)_(.*).fastq', file)
        match_2 = re.search('g-P7M1S-(.*)_(.*)_(.*)_tag.fastq', file)
        match_3 = re.search('g-P7M1S-(.*)-(.*)-combined.fastq', file)
        if match_1:
            barcode, batch, donor = match_1.group(1), match_1.group(2), match_1.group(3)
            file_type = "individual"
            compression_type = False
        if match_2:
            barcode, batch, donor = match_2.group(1), match_2.group(2), match_2.group(3)
            file_type = "individual"
            compression_type = False
        elif match_3:
            barcode, batch, donor = match_3.group(1), "combined", match_3.group(2)
            file_type = "combined"
            compression_type = False
    elif input_compression is True:
        match_3 = re.search('g-P7M1S-(.*)_(.*)_(.*).fastq.gz', file)
        match_4 = re.search('g-P7M1S-(.*)_(.*)_(.*)_tag.fastq.gz', file)
        match_5 = re.search('g-P7M1S-(.*)-(.*)-combined.fastq.gz', file)
        if match_3:
            barcode, batch, donor = match_3.group(1), match_3.group(2), match_3.group(3)
            file_type = "individual"
            compression_type = True
        if match_4:
            barcode, batch, donor = match_4.group(1), match_4.group(2), match_4.group(3)
            file_type = "individual"
            compression_type = True
        elif match_5:
            barcode, batch, donor = match_5.group(1), "combined", match_5.group(2)
            file_type = "combined"
            compression_type = True
    info = File_Info(file,barcode,batch,file_type,compression_type)
    file_list[donor].append(info) #Append the file along with the collected information to a dictionary
    return file_list

def combine_and_write(file_list,summary_file,combine_file_directory,files_type):
    '''
    This function iterates through each file under 1 donor and writes it to the combined file. Subsequently, it writes
    to the summary file telling which files were written as well as how many reads were being written.
    '''
    for donor,file_info in file_list.iteritems():
        if files_type != 'both':
            barcode_list = [info.barcode for info in file_info]
            #We don't want to use barcodes for files that have been split by Decombinator
            barcode_list_filtered = [barcode for barcode in barcode_list if 'a' not in barcode and 'b' not in barcode]
            # We want to use the barcode that has been used the most for the combined file
            bc_count = Counter(barcode_list_filtered)
            sorted_bc_list = bc_count.most_common()
            barcode_to_use = sorted_bc_list[0][0]
        elif files_type == 'both':
            #We want to use the barcode of the combined file
            barcode_to_use = [info.barcode for info in file_info if info.file_type == 'combined'][0]
        new_combine_file = 'g-P7M1S-{}-{}-combined.fastq'.format(barcode_to_use,donor)
        new_combine_file = os.path.join(combine_file_directory,new_combine_file)
        summary_header = "Writing {}...\n".format(new_combine_file)
        summary_header_2 = "Donor\tBarcode\tBatch\tNumber of Reads written to new file\n"
        summary_file.write(summary_header)
        summary_file.write(summary_header_2)
        output_file = open(new_combine_file,'w')

        total_reads =0
        for const_file in file_info:
            read_counter = 0
            if const_file.compression_type is True:
                with gzip.open(const_file.file,'rb') as input:
                    # Write the constituent file to the combined file
                    for title, seq, qual in FastqGeneralIterator(input):
                        output_read =  '@%s\n%s\n+\n%s\n' % (title, seq, qual)
                        output_file.write(output_read)
                        read_counter +=1
                        total_reads +=1
            elif const_file.compression_type is False:
                with open(const_file.file,) as input:
                    for title, seq, qual in FastqGeneralIterator(input):
                        # Write the constituent file to the combined file
                        output_read =  '@%s\n%s\n+\n%s\n' % (title, seq, qual)
                        output_file.write(output_read)
                        read_counter +=1
                        total_reads +=1
            output_summary = "{}\t{}\t{}\t{}\n".format(donor,const_file.barcode,const_file.batch,read_counter)
            summary_file.write(output_summary)
        final_summary_line = "Total Number of reads written to the combined file:{}\n\n".format(total_reads)
        summary_file.write(final_summary_line)
        output_file.close()


def Main():

    parser = argparse.ArgumentParser(description = __doc__,
                                     formatter_class= argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-ip","--input_directory",
                        help = "The directory where the FASTQ files reside in (REQUIRED)",
                        required = True)
    parser.add_argument("-op","--combined_file_output_directory",
                        help = "The directory where the newly combined FASTQ files will reside in (REQUIRED)",
                        required = True)
    parser.add_argument("-sum","--summary_file_output_directory",
                        help="The directory where the summary file will reside in (REQUIRED)",
                        required = True)
    parser.add_argument("-ic","--input_compression",
                        help = "Specify this argument if the input FASTQ files are compressed (OPTIONAL)",
                        action = 'store_true',
                        default = False)
    parser.add_argument("-ft","--file_type",
                        help = "The type of FASTQ file that needs to be combined with (REQUIRED)",
                        choices = ['individual','combined','both'],
                        required = True)
    parser.add_argument('-ft2','--file_type_2',
                        help='Specify what types of files need to be combined; you can input any parameter'
                         'such as naive, memory, reproducibility, etc.',type = 'str',
                        required = True)
    parser.add_argument('-ip2','--second_input_directory',
                        help = "Note: Use this option if you're combining individual files with combined files."
                               "If so, make sure the individual files and combined files are in separate directories.",
                        type = str)

    args = parser.parse_args()

    if not os.path.exists(args.combined_file_output_directory):  # Checks if the output directory exists, and creates it if it doesn't exist.
        try:
            os.makedirs(args.combined_file_output_directory)
        except OSError, e:
            if e.errno != os.errno.EEXIST:
                raise
        time.sleep(1)
        pass

    if not os.path.exists(args.summary_file_output_directory):  # Checks if the output directory exists, and creates it if it doesn't exist.
        try:
            os.makedirs(args.summary_file_output_directory)
        except OSError, e:
            if e.errno != os.errno.EEXIST:
                raise
        time.sleep(1)
        pass

    file_lst= defaultdict(list)
    current_date = datetime.datetime.now().strftime("%y-%m-%d")
    summary_file_name = os.path.join(args.summary_file_output_directory + 'Combined_Files_Summary-{}-{}.txt'.format(args.file_type_2,current_date))
    summary_file = open(summary_file_name,'w')  # The summary file will show what are the constituents of each combined file are

    if args.file_type != "both":
        files_to_be_combined = sorted(glob.glob(args.input_directory + '*.fastq*'))

        for fastq_file in files_to_be_combined:
            file_lst = getinfo(fastq_file,file_lst,args.input_compression)

    elif args.file_type =="both":
        indiv_files = sorted(glob.glob(args.input_directory + "*.fastq*"))
        combined_files = sorted(glob.glob(args.second_input_directory + '*.fastq*'))

        for sub_file in indiv_files:
            file_lst = getinfo(sub_file,file_lst,args.input_compression)

        for big_file in combined_files:
            file_lst = getinfo(big_file,file_lst,args.input_compression)

    combine_and_write(file_lst,summary_file, args.combined_file_output_directory, args.file_type)


    summary_file.close()


if __name__ == "__main__":
    Main()

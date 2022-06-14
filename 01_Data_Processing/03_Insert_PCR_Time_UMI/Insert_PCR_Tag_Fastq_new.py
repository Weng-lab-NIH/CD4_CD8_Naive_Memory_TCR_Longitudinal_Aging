"""
FILE NAME:      Insert_PCR_Tag_FASTQ_new.py

COMPLEMENTARY
SCRIPT(S)/
FILE(S):        submit_PCR_tag.sh; batch_info_pcr_tags.txt

DESCRIPTION:    Inserts PCR Tag (converted from PCR Time) into the end of each UMI.

COMPATIBLE
PIPELINES:		Human, Mouse

INPUT
FILE(S):        1) Individual FASTQ file (reformatted FASTQ file with barcode, donor, and batch in the file name).
                2) batch_info_pcr_tags.txt: File that contains the PCR tag to be embedded into the UMI.

OUTPUT
FILES(S):       1) Tagged FASTQ file that will reside in the <Input_Directory>/PCR_tagged directory.

INPUT
ARGUMENTS:      1) Input Directory: Directory where the untagged FASTQ files reside in.
                2) PCR Master list: Tab-delineated .txt file that lists the PCR Tag for each individual file.

CREATED:        01Nov17

MODIFICATIONS
LOG:

15Nov17         Prints out file name that correctly got PCR tagged.
17Jan18         Changed PCR Tagging from end of each FASTQ header to end of UMI
02Jul18         Reformatted code for better language & incorporation of argparse.

PYTHON VERSION
USED TO
WRITE
SCRIPT:         2.7.14

VERSION:        1.0

AUTHOR(S):      Thomas Nguyen

STATUS:         Working

DISCLAIMER(S):  The mouse pipeline does not currently utilize PCR Times.
"""

import glob, re, os, sys, time,argparse,gzip
from Bio.SeqIO.QualityIO import FastqGeneralIterator


def batchwriter(input_file,output_directory,pcr_tag,input_compression,output_compression):
    """
    This function appends the PCR Tag (3bp nucleotide) onto the end of the UMI in each read.
    """
    file_name = os.path.basename(input_file)
    if input_compression == 'N':
        output_filename = file_name.replace('.fastq','_tag.fastq')
    else:
        output_filename = file_name.replace('.fastq.gz', '_tag.fastq.gz')
    output_filename = os.path.join(output_directory,output_filename)
    read_count = 0
    if output_compression == 'N':
        outfile = open(output_filename,'w')
    else:
        outfile = gzip.open(output_filename, 'wb')

    if output_compression == 'N':
        with open(input_file) as input:
            try:
                for title,seq,qual in FastqGeneralIterator(input):
                    read_count +=1
                    info = title.rstrip().split(' ')
                    umi_info = info[1]
                    match = re.search(r"UMI:(\w{16}):(.*)",umi_info)
                    umi,umi_qual = match.group(1),match.group(2)
                    if not match:
                        print "Error! The UMI can not be found in the FASTQ read ID. Please check Read {}!".format(
                            read_count + 1)
                        sys.exit(1)
                    new_umi = umi + pcr_tag
                    new_umi_info = 'UMI:{}:{}'.format(new_umi,umi_qual)
                    new_title = ' '.join([info[0],new_umi_info])
                    read = '@{}\n{}\n+\n{}\n'.format(new_title,seq,qual)
                    outfile.write(read)
            except ValueError:
                print "Error! Read {] has issues! Please check this read!".format(read_count + 1)
    else:
        with gzip.open(input_file,'rb') as input:
            try:
                for title,seq,qual in FastqGeneralIterator(input):
                    read_count +=1
                    info = title.rstrip().split(' ')
                    umi_info = info[1]
                    match = re.search(r"UMI:(\w{16}):(.*)",umi_info)
                    umi,umi_qual = match.group(1),match.group(2)
                    if not match:
                        print "Error! The UMI can not be found in the FASTQ read ID. Please check Read {}!".format(
                            read_count + 1)
                        sys.exit(1)
                    new_umi = umi + pcr_tag
                    new_umi_info = 'UMI:{}:{}\n'.format(new_umi,umi_qual)
                    new_title = ' '.join(info[0],new_umi_info)
                    read = '@{}\n{}\n+\n{}\n'.format(new_title,seq,qual)
                    outfile.write(read)
            except ValueError:
                print "Error! Read {] has issues! Please check this read!".format(read_count + 1)
    outfile.close()


def Main():
    command_line = argparse.ArgumentParser(description=__doc__,
                                           formatter_class=argparse.RawDescriptionHelpFormatter)
    command_line.add_argument("input_directory",
                              help="The completed filler match FASTQ file that needs to be renamed")
    command_line.add_argument("pcr_master_list",
                              help="The tab-delineated .txt file that contains the PCR tags for each individual file")
    command_line.add_argument("input_compression",
                              help = "Specify whether the input files are compressed or not",
                              choices = ['N',"Y"])
    command_line.add_argument("output_compression",
                              help = "Specify whether the output files need to be compressed or not",
                              choices = ['N',"Y"])
    args = command_line.parse_args()

    if not args.input_directory[-1] == "/":  # If the use doesn't add the slash at the end of the argument, warn them.
        print "Error! Make sure to end your Python argument path with a /!"
        sys.exit(1)

    output_directory = os.path.join(args.input_directory + 'PCR_tagged/')

    if not os.path.exists(output_directory):  # Create the output directory if it doesn't already exist.
        try:
            os.makedirs(output_directory)
        except OSError, e:
            if e.errno != os.errno.EEXIST:
                raise
            time.sleep(1)
            pass

    input_fastq_files = glob.glob(args.input_directory + '*.fastq')

    #Keep track of which files did
    files_tagged = 0
    files_to_be_tagged = len(input_fastq_files)
    master_list_length = sum(1 for line in open(args.pcr_master_list) if "PCR TIME" not in line)

    for file in input_fastq_files:

        if args.input_compression == 'N':
            match = re.search(r"g-P7M1S-(.*)_(.*)_(.*).fastq",file)
        else:
            match = re.search(r"g-P7M1S-(.*)_(.*)_(.*).fastq.gz", file)
        f_barcode,f_batch,f_donor = match.group(1),match.group(2),match.group(3)

        if "_" in f_barcode:    #During old runs of PCR Tagging, Barcode included an underscore in them
            f_barcode = f_barcode.replace('_','')


        for line in open(args.pcr_master_list,'r'):
            tags = line.rstrip().split('\t')
            master_donor = tags[0]
            barcode = tags[1]
            if "P7M1S-" in barcode:						# some barcodes have P7M1S in the barcode name
                master_barcode = barcode[6:]
            else:
                master_barcode = barcode
            master_batch = tags[2]
            pcr_time = tags[3]


            #Check to see if the file barcode, batch, and donor matches with the info from the master list
            if master_barcode == f_barcode and master_batch == f_batch and master_donor == f_donor:
                batchwriter(file,output_directory,pcr_time,args.input_compression,args.output_compression)
                files_tagged +=1

    print "Number of Files in the Input Directory: {}".format(files_to_be_tagged)
    print "Number of Files that ended up getting PCR tagged: {}".format(files_tagged)
    print "Number of Files that needed to be tagged from the Master List: {}".format(master_list_length)



if __name__ == "__main__":
    Main()

"""
FILE NAME:      gc_reformat-filler.py

COMPLEMENTARY
SCRIPT(S)/
FILE(S):        generate-gc-ref.py; submit-gc_ref-filler.sh

DESCRIPTION:    Take raw sequences generated from Illumina sequencing machine and convert the FASTQ file
				into a MiGEC-compatible file.

COMPATIBLE
PIPELINES:		Human, Mouse

INPUT
FILE(S):        1) Individual FASTQ file (unformatted FASTQ file generated from end of the filler match step
				or rename step.

OUTPUT
FILES(S):       1) MiGEC-compatible FASTQ file.

INPUT
ARGUMENTS:      1) Input File Name: Specify the full path of the file name that will be read by this script.
                2) Output File Name: Specify the full path of the file name that will be generated
                by this script.

CREATED:        May16

MODIFICATIONS
LOG:			None

PYTHON VERSION
USED TO
WRITE
SCRIPT:         2.6

VERSION:        1.0

AUTHOR(S):      Annette Ko, Jay Sharma

STATUS:         Working

NOTES:			None
"""


import argparse
from Bio import SeqIO
import time 

def generate_fastq(record):
	"""
	This function converts the raw sequence files to a standardized format
	accepted by MiGEC.
	"""
	dna_seq, dna_qual, umi_seq, umi_qual, seq_id = parse(record)
	record.description = '%s UMI:%s:%s'%(seq_id, umi_seq, umi_qual)
	record.id, record.name = seq_id, seq_id
	return record

def parse(record):
	"""
	This function sifts through the raw sequence file to obtain the DNA sequence,
	its associated quality score, UMI, UMI quality score, and Illumina header.
	"""
	split1, split2 = record.description.split(' '), record.description.split('\t')
	dna_seq, dna_qual, umi_seq, umi_qual, seq_id = str(record.seq), record.letter_annotations['phred_quality'], split2[2][0:16], split1[-1][0:16], split1[0]
	return dna_seq, dna_qual, umi_seq, umi_qual, seq_id

def main():
	parser = argparse.ArgumentParser(description = 'Reformatting tool for custom decombinator outputs back into fastq')
	parser.add_argument('-in', '--input', type = str, required=True, help = 'Input file path')
	parser.add_argument('-out', '--output', type = str, required=True, help = 'Output file path')
	args = parser.parse_args()

	in_file = args.input
	output = open(args.output, 'w')
 
	for record in SeqIO.parse(in_file, 'fastq'):
		fastq_obj = generate_fastq(record)
		SeqIO.write(fastq_obj, output, "fastq")	#Write standardized MiGEC-compatible object to new file


if __name__ == '__main__':
	start_time = time.time()
	main()
	print("--- %s seconds ---" % (time.time() - start_time)) #Take note of how long processing takes


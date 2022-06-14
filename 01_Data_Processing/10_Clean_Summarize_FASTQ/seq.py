"""
FILE NAME:      seq.py

COMPLEMENTARY
SCRIPT(S)/
FILE(S):        None

COMPATIBLE
PIPELINE(S):    Human, Mouse

DESCRIPTION:    This script will summarize the number of unique TCRs (V-CDR3 AA-J) each sample has.

INPUT
FILES(S):       1) Fixed Cut File

OUTPUT
FILES(S):       1) A summary file

INPUT
ARGUMENT(S):    1) Input Directory: The directory that contains the fixed cut files.


CREATED:        May2016

MODIFICATIONS
LOG:
07Aug19			Modified script to be compatible with Python 3.7


LAST MODIFIED
BY:             Thomas Nguyen

PYTHON
VERSION USED
TO WRITE
SCRIPT:         2.7.14 --. 3.7

VERSION:        3.0

AUTHOR(S):      Annette Ko; Jay Sharma; Chen Chen; Thomas Nguyen

STATUS:         Working

TO DO LIST:     None at this moment.

DISCLAIMER(S):

NOTE(S):
"""


import re, sys,  glob, time, ntpath

input_dir = sys.argv[1]
files_dir = sorted(glob.glob(input_dir + "*fixed*.txt")) 

summary_header = "%s\t%s\t%s\t%s\t%s\t%s"%("barcode","donor","chain","total TCR","total UMI","total reads")
output_file = open('summary_counts.txt','w')
output_file.write(summary_header+"\n")

start_time = time.time()
for _file in files_dir:
	head, filename = ntpath.split(_file)
	match = re.search('(.*)P7M1S-(.*)-(.*)-(.*)combinedfixed(.*).txt', filename)
	barcode = match.group(2)
	donor = match.group(3)
	chain = match.group(4)
	tot_umi_count = 0
	tot_reads = 0
	tot_tcr = 0
	f = open(_file, 'r')
	for line in f:
		line = line.rstrip()
		components = line.split()
		tot_umi_count += int(components[-2]) #Running count of Total UMI Count
		tot_reads += int(components[-1])	#Running count of Total Reads
		tot_tcr += 1	#each line is a unique TCR
	f.close()

	summary_line = "%s\t%s\t%s\t%s\t%s\t%s"%(barcode,donor,chain, int(tot_tcr), int(tot_umi_count), int(tot_reads))
	print (summary_line)
	output_file.write(summary_line+"\n")
output_file.close()
print ("---", time.time() - start_time, "seconds---")
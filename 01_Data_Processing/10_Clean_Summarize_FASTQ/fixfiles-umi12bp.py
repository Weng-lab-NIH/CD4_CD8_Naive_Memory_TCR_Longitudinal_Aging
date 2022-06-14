"""
FILE NAME:      fixfiles-umi12bp.py

COMPLEMENTARY
SCRIPT(S)/
FILE(S):        None

COMPATIBLE
PIPELINE(S):    Human, Mouse

DESCRIPTION:    This script will eliminate TCRs with stop codons, frame shifts, and summarize the UMI/Reads count
				for each TCR.

INPUT
FILES(S):       1) Merged File

OUTPUT
FILES(S):       1) A Fixed File

INPUT
ARGUMENT(S):    1) Input Directory: The directory that contains the merged files.

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

############################################################################################################
# DESCRIPTION: ELIMINATE TCRS WITH WRONG / BAD CDR3 AA (STOP CODONS, FRAMESHIFT, ETC) & SUMMARIZE
# NOW EACH LINE BECOMES: TCR unique by V-CDR3 AA-J or V-CDR3 NT-J with UMI counts, read counts
# INPUT: CURRENT DIRECTORY (DIR MUST CONTAIN 'merged-*txt' or 'c-*txt' (c for cleaned)
# OUTPUT: 'P7M1S-*fixed.txt' files in 'fixed_uniqueAA', 'fixed_uniqueNT' dirs, 'results-fixfiles.txt' in ./
# LAST MODIFIED: 11/3/16
############################################################################################################

import os, sys, re, glob
from collections import defaultdict

dir = sys.argv[1]
files = sorted(glob.glob(dir+'merged-*.txt'))

def getinfo(file):
	match = re.search("(.*)P7M1S-(.*)-(.*)-(.*)combined.txt", file)
	barcode, donor, chain = match.group(2), match.group(3), match.group(4)
	return barcode, donor, chain
def checkcdr3(cdr3aa):
	if re.search(r'\*|\?', cdr3aa):
		return False
	return True
def clean(cdr3aa, v, chain):
	check = checkcdr3(cdr3aa)
	if check:
		if chain=='alpha':
			if cdr3aa[0]!='C':
				return False
		elif chain=='beta':
			if cdr3aa[0]=='C':
				if 'BV17' in v or 'BV26' in v:
					return False	
			elif cdr3aa[0]=='Y':
				if 'BV17' not in v or 'BV26' not in v:
					return False	
		return True
	return False
def parse(line):
	line = line.rstrip()
	split = line.split('\t')
	cdr3nt, v, j, d, cdr3aa, umi, tcr150, rds = split[0], split[1], split[2], split[3], split[4], split[5], split[6], int(split[-1])
	if len(umi) == 16:
		umi = umi[1:5]+umi[6:10]+umi[11:15]
	if len(umi) > 16:
		pass
	return cdr3nt, v, j, d, cdr3aa, umi, tcr150, rds
def makedicts(file):
	barcode, donor, chain = getinfo(file)
	dict_aa, dict_nt = defaultdict(list), defaultdict(list)
	dict_aa_rds, dict_nt_rds = defaultdict(int), defaultdict(int)
	for line in open(file, 'r'):
		cdr3nt, v, j, d, cdr3aa, umi, tcr150, rds = parse(line)
		if clean(cdr3aa, v, chain):
			tcr_aa = '%s\t%s\t%s'%(v, cdr3aa, j)
			tcr_nt = '%s\t%s\t%s\t%s'%(v, cdr3aa, j, cdr3nt)
			dict_aa[tcr_aa] += [umi]
			dict_nt[tcr_nt] += [umi]
			dict_aa_rds[tcr_aa] += int(rds)
			dict_nt_rds[tcr_nt] += int(rds)
	return dict_aa, dict_nt, dict_aa_rds, dict_nt_rds
def write(newfilename, dirname, dict, dict_rds):
	filename = dirname+'/'+newfilename
	outfile = open(filename, 'w')
	tot_umi, tot_rds = 0, 0
	for tcr, umilist in dict.iteritems():
		umi, rds = len(set(umilist)), dict_rds[tcr]
		tot_umi += umi
		tot_rds += rds
		line = '%s\t%s\t%s'%(tcr, umi, rds)
		outfile.write(line+'\n')
	outfile.close()
	return filename, len(dict), tot_umi, tot_rds

dirlist = ['fixed_uniqueNT', 'fixed_uniqueAA']

for dir in dirlist:
	if not os.path.exists(dir):
		os.makedirs(dir)
	
outfile = open('results-fixfiles.txt','w')
header = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%('barcode','donor','chain','tot_tcr_aa','tot_umi_aa','tot_rds_aa','tot_tcr_nt','tot_umi_nt','tot_rds_nt')
outfile.write(header+'\n')
print (header)
for file in files:
	dict_aa, dict_nt, dict_aa_rds, dict_nt_rds = makedicts(file)	
	barcode, donor, chain = getinfo(file)
	if len(dict_aa)>0:
		newfilename = 'P7M1S-'+barcode+'-'+donor+'-'+chain+'combinedfixed.txt'
		filename_aa, tot_tcr_aa, tot_umi_aa, tot_rds_aa = write(newfilename, dirlist[1], dict_aa, dict_aa_rds)
		filename_nt, tot_tcr_nt, tot_umi_nt, tot_rds_nt = write(newfilename, dirlist[0], dict_nt, dict_nt_rds)
		summary = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%(barcode,donor,chain,tot_tcr_aa,tot_umi_aa,tot_rds_aa,tot_tcr_nt,tot_umi_nt,tot_rds_nt)
		print (summary)
		outfile.write(summary+'\n')


"""
FILE NAME:      contam-V3-parallel.py

COMPLEMENTARY
SCRIPT(S)/
FILE(S):        generate-contam-parallel.py; submit-contam-para.sh

COMPATIBLE
PIPELINE(S):    Human, Mouse

DESCRIPTION:    For a collapsed consensus 100bp seqeunce, collapsed 12bp UMI, and CDR3 AA in a particular sample,
				if the UMIs and CDR3 are found in another sample and if the consensus sequence is either a perfect match
				or at most 3bp mismatch away, the "contaminant" sequence will be written to this file.

INPUT
FILES(S):       1) Compact Merged File

OUTPUT
FILES(S):       1) Overlap File that lists which sequences were found in other samples.

INPUT
ARGUMENT(S):    1) Input Directory: The directory that contains the compact merged files.
				2) Input File: The compact merged file that will serve as a reference to compare all other files with.


CREATED:        Jan2016

MODIFICATIONS
LOG:
06Apr18			Changed Column Order to extract correct information


LAST MODIFIED
BY:             Jay Sharma

PYTHON
VERSION USED
TO WRITE
SCRIPT:         2.6

VERSION:        2.0

AUTHOR(S):      Annette Ko; Jay Sharma; Chen Chen; Thomas Nguyen

STATUS:         Working

TO DO LIST:     Clean up script, optimization of functions. Could potentially rewrite

DISCLAIMER(S):

NOTE(S):
"""

import os, sys, ntpath, glob, time, Levenshtein as lev
from collections import defaultdict

dir = sys.argv[1]
my_file = [sys.argv[2]]
files = sorted(glob.glob(dir+'compact-merged-*'))
chainlist = ['alpha', 'beta']
mis = 2


def parse(line):
	"""
	This function reads through the compact merged file
	and extracts the corresponding information.
	"""
	line = line.rstrip()
	split = line.split('\t')
	aa, umi, seq = split[2], split[1], split[0][:100]  # Adjusted for use with new compact files.
	return aa, umi, seq


def makedict(file):
	"""
	This function returns a dictionary of the CDR3 AA and UMI (TCR information that requires a perfect match) along
	with the corresponding merged file information.
	"""
	dict = defaultdict(list)
	cdr3list = set()
	for line in open(file, 'r'):
		aa, umi, seq = parse(line)
		if len(umi) > 12:
			umi = umi[1:5]+umi[6:10]+umi[11:15]
		cdr3list.add(aa)
		key = '%s\t%s'%(aa, umi)
		if line not in dict[key]:
			dict[key] += [line]
	return dict, len(cdr3list)


def compareseq(list1, list2):
	"""
	This function checks the hamming distance between 2 sequences. If it is at most 3bp, it will append it to
	a list for further analysis.
	"""
	linelist, count = [], []
	for line1 in list1:
		for line2 in list2:
			aa1, umi1, seq1 = parse(line1)
			aa2, umi2, seq2 = parse(line2)
			if aa1==aa2:
				hamming = lev.hamming(seq1[0:100], seq2[0:100])
				if hamming < mis:
					count.append(line1)
					line1 = line1.rstrip()
					line2 = line2.rstrip()
					newline = '%s\t%s\t%s'%(line1, line2, hamming)
					linelist.append(newline)
	return list(set(linelist)), list(set(count))


def sumdict(dict):
	"""
	This function returns the total number of sequences stored in the dictionary.
	"""
	sum = 0
	for key, vallist in dict.iteritems():
		sum += len(set(vallist))
	return sum


def getolcdr3(keylist):
	"""
	This function returns all overlapping CDR3 AA's.
	"""
	newlist = list()
	for key in keylist:
		split = key.split('\t')
		newlist.append(split[0])
	return list(set(newlist))


def getolumi(keylist):
	"""
	This function returns all overlapping UMIs.
	"""
        newlist = list()
        for key in keylist:
                split = key.split('\t')
                newlist.append(split[1])
        return list(set(newlist))


start = time.time()
header = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%('file1', 'file2','cdr31','umi1','rds1','cdr32','umi2','rds2','ol-cdr3','ol-umi','ol-dna')
print header
for file in my_file:
	head, filename = ntpath.split(file)
	resultf = 'results-ol-'+filename
	resultf2 = 'results-ol-summary-'+filename
	resultfile = open(resultf, 'w')
	resultfile.write(header+'\n')
	resultfile2 = open(resultf2, 'w')
	contam, cdr3c, umic, dnac = [], [], [], []
	chain = None
	dict1, totaa1 = makedict(file)
	if 'alpha' in file:
		chain = 'alpha'
	elif 'beta' in file:
		chain = 'beta'	
	for j, file2 in enumerate(files):
		if chain in file2 and file not in file2:
				start1 =time.time()
				head, filename2 = ntpath.split(file2)
				dict2, totaa2  = makedict(file2)
				info = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%(filename, filename2, totaa1, len(dict1),
														 sumdict(dict1), totaa2, len(dict2), sumdict(dict2))
				ol = set(dict1.keys()).intersection(dict2.keys())	#Find overlapping CDR3's and UMI's
				cdr3ol = getolcdr3(ol)
				umiol = getolumi(ol)
				ollinelist = []
				count_dna_contam = []
				if len(ol)>0:
					for key in ol:
						# Check the hamming distance between the consensus sequences
						ollinelist, count_dna_contam = compareseq(dict1[key], dict2[key])
						if len(ollinelist)>0:
							for line in ollinelist:
								newline = '%s'%(line)
								contam.append(newline)
				result = '%s\t%s\t%s\t%s'%(info, len(set(cdr3ol)), len(umiol),
										   len(count_dna_contam)) # cdr3 level, umi level, dna level
				cdr3c += (cdr3ol)
				umic += (umiol)
				if len(count_dna_contam)>0:
					dnac += (count_dna_contam)	#Count how many true contaminated DNA sequences there are
				resultfile.write(result+'\n')
				print '%s\t%0.2f  SEC =='%(result,time.time()-start1)
	result2 = '%s\t%s\t%s\t%s\t%s\t%s\t%s'%(filename, totaa1, len(dict1), sumdict(dict1), len(set(cdr3c)), len(set(umic)), len(set(dnac)))
	print result2
	resultfile2.write(result2+'\n')
	if len(contam)>0:
		outfilename = 'ol-'+filename
		outfile = open(outfilename, 'w')
		for contamline in contam:
			outfile.write(contamline+'\n')
		outfile.close()
	resultfile.close()
	resultfile2.close()

os.system("cat results-ol-compact-merged-* > results-ol.txt && cat results-ol-summary-compact-merged-* > results-ol-summary.txt")
		
print '== TOTAL TIME TAKEN: %0.2f SEC =='%(time.time()-start)
			
	

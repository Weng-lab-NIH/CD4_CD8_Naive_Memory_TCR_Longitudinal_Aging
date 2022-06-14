"""
FILE NAME:     filler_match_v6.py

COMPLEMENTARY
SCRIPT(S)/
FILE(S):        matchbarcodes.py; submitfiller6.sh; fillerscript6.sh

DESCRIPTION:    Identify files that need to be converted into a MiGEC-compatible format.

COMPATIBLE
PIPELINES:		Human, Mouse

INPUT
FILE(S):        1) None

OUTPUT
FILES(S):       1) Submission script for each file recognized

INPUT
ARGUMENTS:      1) Input Directory: Specify where the FASTQ files generated from the filler match
				or the rename step reside in.

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

NOTES:
"""


#take decombinator output and then analyze UMI
import os, sys, glob, time
import Levenshtein as lev
from itertools import izip
from Bio import SeqIO
from acora import AcoraBuilder


def new_record(record_seq, UMI_seq, barcode_seq, mismatch, UMI_q, barcode_q, umi_phred):
	return_record = record_seq
	return_record.description = barcode_seq + "\t" + mismatch + "\t" + UMI_seq + "\t" + str(UMI_q) + "\t" + str(barcode_q) + "\t" + str(record.seq) + " " + str(umi_phred)
	return return_record

def match(filler_seq, find_seq):
	"""
	This function returns boolean and hamming distance. True if hamming < 8, False if otherwise.
	"""
	#print len(filler_seq)
	#print len(find_seq)
	hamming = lev.hamming(filler_seq, find_seq)
	if hamming < 8:
		return True, hamming
	else:
		return False, hamming

def find_filler(record_seq):
	"""
	This function returns a boolean, where the barcode is in the sequence,
	and the hamming distnace bt barcode (19b "filler_seq") and the first matching sequence.
	"""
	global barcode_acora, filler_seq
	seq = record_seq.seq
	search_filler = barcode_acora.findall(str(record.seq))
	found = False
	if search_filler:
		start_loc = search_filler[0][1]
		return True, start_loc, 0
	else: 
		for idx in xrange(8, 12):
			to_match_seq = seq[idx:idx+19]
			to_match_seq = ''.join(to_match_seq)
			_bool, hamming = match(filler_seq, to_match_seq)
			if _bool and hamming < 8:
				return True, idx, hamming
	return False, 0, 999999

def convert(list):
	return ''.join([chr(int(x)+33) for x in list])

def chop_UMI_barcode(record, record2, filler_start_loc):
	# reads returned properties from file and converts to unicode when desired
	# q indicates quality of nucleotide (phred score)
	seq = record.seq
	UMI_q = record.letter_annotations["phred_quality"][filler_start_loc+19:filler_start_loc+35]
	UMI = seq[filler_start_loc+19:filler_start_loc+35]
	UMI = ''.join(UMI)
	barcode = seq[filler_start_loc-8:filler_start_loc]
	barcode = ''.join(barcode)
	barcode_q = record.letter_annotations["phred_quality"][filler_start_loc-8:filler_start_loc]
	TCR_q = record2.letter_annotations['phred_quality']
	str_UMI_q = convert(UMI_q)
	str_barcode_q = convert(barcode_q)
	str_TCR_q = convert(TCR_q)
	return barcode, UMI, barcode_q, UMI_q, TCR_q, str_barcode_q, str_UMI_q, str_TCR_q

def count_N(UMI):
	count = 0
	for char in list(UMI):
		if char == 'N':
			count += 1
	return count

def score_barcode(barcode, errors):
	"""
	This function goes through barcodes from new_barcodes.txt,
	returns if sequenced barcode is in list, the best hamming distance,...
	a "last local hamming" distance defined as len(_barcode[i])-8,
	and the index of the expected barcode in _barcodes.
	"""
	global _barcodes
	found = 0 
	local_minimum_hamming = 9999999
	last_local_hamming = 9999999
	last_barcode_id = 99999999
	for idx, val in enumerate(_barcodes):
		if len(val) == 8:
			_comp = val
			_comp = ''.join(_comp)	#turn barcode into string
			# hamming distance between mean read barcode and (expected) barcode
			# from new_barcodes.txt
			hamming_seg1 = lev.hamming(_comp, barcode)
			
			# seq sequenced barcode with min.
			# Hamming from expected to be last_barcode_id.
			if hamming_seg1 < errors:
				found = 1
				if hamming_seg1 < local_minimum_hamming:
					local_minimum_hamming = hamming_seg1
					last_local_hamming = 0
					last_barcode_id = idx 
		if len(val) == 9:
			_comp = val[1:]
			_comp = ''.join(_comp)
			hamming_seg1 = lev.hamming(_comp, barcode)
			if hamming_seg1 < errors:
				found = 1
				if hamming_seg1 < local_minimum_hamming:
					local_minimum_hamming = hamming_seg1
					last_local_hamming = 1
					last_barcode_id = idx 
		elif len(val) == 10:
			_comp = val[2:]
			_comp = ''.join(_comp)
			hamming_seg2 = lev.hamming(_comp, barcode)
			if hamming_seg2 < errors:
				found = 2
				if hamming_seg2 < local_minimum_hamming:
					local_minimum_hamming = hamming_seg2
					last_local_hamming = 2
					last_barcode_id = idx
	return found, local_minimum_hamming, last_local_hamming, last_barcode_id

def score_UMI_T(UMI):
	score_1, score_2, score_3, score_4 = 0,0,0,0
	if UMI[0] == 'T':
		score_1 += 1
	if UMI[5] == 'T':
		score_2 += 1
	if UMI[10] == 'T':
		score_3 += 1
	if UMI[15] == 'T':
		score_4 += 1
	return score_1, score_2, score_3, score_4

def score_UMI_N(UMI):
	score_1, score_2, score_3, score_4 = 0,0,0,0
	if UMI[0] == 'N':
		score_1 += 1
	if UMI[5] == 'N':
		score_2 += 1
	if UMI[10] == 'N':
		score_3 += 1
	if UMI[15] == 'N':
		score_4 += 1
	return score_1, score_2, score_3, score_4

#################################################################################################

start = time.time()

barcode_builder = AcoraBuilder()
barcode_builder.add('CAGTGGTATCAACGCAGAG')			#Conserved sequence in each barcode
barcode_acora = barcode_builder.build()
filler_seq = 'CAGTGGTATCAACGCAGAG'

dir = sys.argv[1]
errors = int(sys.argv[2])	#Maximum of 2 errors

_input = glob.glob(dir+'*_R2*fastq')[0]
input_file = open(_input, 'r')

_input2 = glob.glob(dir+'*_R1*fastq')[0]
input_file_2 = open(_input2, 'r')

os.system("python matchbarcodes.py ./")

barcodefile = glob.glob(dir+'new_barcodes.txt')[0]
tags = open(barcodefile, 'r')			#This file contains the barcodes tags

_barcodes = []
_barcodes_index = []
_barcodes_names = []
for line in tags:
	# Reads file titled "new_barcodes.txt" in specified directory and adds right column to _barcodes,...
	# ...the index of the barcode in the second column to _barcodes_index, and left column of the file...
	# ...to _barcodes_names.
	line = line.rstrip()
	components = line.split()
#	print(components)
	if len(components) > 0:
		if components[1][1] == 'n':
			_barcodes.append(components[1][:10])	#use only tags associated with each barcode
			_barcodes_index.append(2)
			_barcodes_names.append(components[0])
		#	print "Appending %s to barcodes!\n"%(components[1][:10])
		elif components[1][2] == 'n':
			_barcodes.append(components[1][:9])
			_barcodes_index.append(1)
			_barcodes_names.append(components[0])
		#	print "Appending %s to barcodes!\n"%(components[1][:9])	
		elif len(components[1]) > 2:
			_barcodes.append(components[1])
			_barcodes_index.append(0)
			_barcodes_names.append(components[0])
		#	print "Appending %s to barcodes!\n"%(components[1][:9])					

print "Size of barcodes array:\t%i\n"%(len(_barcodes))

### calculate minimum hamming distance between barcodes ... ###

minimum_hamming = 999999999
for idx, val in enumerate(_barcodes):
	for idx2, val2 in enumerate(_barcodes):
		if idx<idx2: 
			current_hamming = lev.hamming(val, val2)
			if current_hamming < minimum_hamming:
				minimum_hamming = current_hamming	
print "Minimum Hamming Between Barcodes:\t%i\n"%(minimum_hamming)

### initalize variables for hamming distance calculations for barcodes ###

num_records = 0
num_found = 0
num_ambiguous = 0
num_unfound = 0
num_T0 = 0
num_T5 = 0
num_T10 = 0
num_T15 = 0
limit = 8
barcode_counter_array = [0 for x in xrange(len(_barcodes))]	#how many sequences best match this barcode
barcode_found_index_array = [0 for x in xrange(11)]
N_distibution = [0 for x in xrange(16)]
filler_hamming_distribution = [0 for x in xrange(11)]
barcode_hamming_array = [0 for x in xrange(errors)]
insufficient_length = 0
_sample_print_array = []
for idx, val in enumerate(_barcodes):
	_sample_print_array.append([])	#each element matching a _barcode will be filled with (SeqIO) record objects

for record, record2 in izip(SeqIO.parse(input_file, "fastq"), SeqIO.parse(input_file_2, "fastq")):
	num_records += 1
	_filler_bool, _filler_loc, filler_hamming = find_filler(record)
	if _filler_bool and _filler_loc >= 8:
		# if the first record has a match to filler_seq and this takes place past 9b in the sequence...
		num_ambiguous += 1
		filler_hamming_distribution[filler_hamming] += 1
		barcode, UMI, barcode_q, UMI_q, TCR_q, str_barcode_q, str_UMI_q, str_TCR_q = chop_UMI_barcode(record, record2, _filler_loc)
		found, local_minimum_hamming, last_local_hamming, last_barcode_id = score_barcode(barcode, errors)
		if found and len(UMI)>=16 and len(barcode) >= 8:
			outputname1 = 'd-'+ _barcodes_names[last_barcode_id] + '_R1.fastq'
			outputname2 = 'd-'+_barcodes_names[last_barcode_id] + '_R2.fastq'
	#		outputname3 = 'QS-'+ _barcodes_names[last_barcode_id] + '.txt'
			output1 = open(outputname1, 'a')
			output2 = open(outputname2, 'a')
	#		output3 = open(outputname3, 'a')
			SeqIO.write(record2, output1, 'fastq')
			SeqIO.write(record, output2, 'fastq')
			newline = '%s\t%s\t%s\t%s\t%s\t%s\t%s'%(record.name, barcode, str_barcode_q, UMI, str_UMI_q, str(record2.seq), str_TCR_q)
	#		output3.write(newline+'\n')

			_Ns = count_N(UMI)
			if _Ns >= 0 and _Ns < 16:
				N_distibution[_Ns] += 1
			UMI_q_average = sum(UMI_q)/float(len(UMI_q))
			barcode_q_average = sum(barcode_q)/float(len(barcode_q))
			umi_phred = convert(UMI_q)
			num_found += 1
			_1, _2, _3, _4 = score_UMI_T(UMI)
			_a, _b, _c, _d = score_UMI_N(UMI)
			UMI = ''.join(UMI)
			_record = new_record(record2, UMI, barcode, str(local_minimum_hamming),UMI_q_average, barcode_q_average, umi_phred)
			_sample_print_array[last_barcode_id].append(_record) # append to appropriate barcode index the record object
			barcode_counter_array[last_barcode_id] += 1
			barcode_found_index_array[last_local_hamming] += 1
			barcode_hamming_array[local_minimum_hamming] += 1
			num_T0 += _1
			num_T5 += _2
			num_T10 += _3
			num_T15 += _4
		else:
			num_unfound += 1

	if num_records % 500000 == 0:
		print '%s\t%s\t%s\t%s\t%s\t%s\t%s'%('# records','# found','# ambiguous','# unfound','% found', '% ambiguous','% unfound')
		print "%i\t%i\t%i\t%i\t%0.03f\t%0.03f\t%0.03f"%(num_records, num_found, num_ambiguous, num_unfound, 100*num_found/float(num_records), 100*num_ambiguous/float(num_records), 100*num_unfound/float(num_ambiguous))
		#	print "\nHamming Array Stats"
		#	for idx, val in enumerate(barcode_hamming_array):
		#		print "%i\t%i\t%0.03f"%(idx, val, val/float(num_found))

		#	print "\nBarcode Found Index Array Stats"
		#	for idx, val in enumerate(barcode_found_index_array):
		#		print "%i\t%i\t%0.03f"%(idx, val, val/float(num_found))

		#	print "\nBarcode Counter Array Stats"
		#	for idx, val in enumerate(barcode_counter_array):
		#		print "%i\t%i\t%0.03f"%(idx, val, val/float(num_found))

		#	print "\nUMI N distribution"
		#	for idx, val in enumerate(N_distibution):
		#		print "%i\t%i\t%0.03f"%(idx, val, val/float(sum(N_distibution)))

		#	print "\nFiller Hamming"
		#	for idx, val in enumerate(filler_hamming_distribution):
		#		print "%i\t%i\t%0.03f"%(idx, val, val/float(sum(filler_hamming_distribution)))

		#	print "\nUMI error checking"
		#	print '%s\t%s\t%s\t%s\t%s'%('% T0/found','% T5/found', '% T10/found','% T15/found','% insufficient length/found')
		#	print "%0.03f\t%0.03f\t%0.03f\t%0.03f\t%0.05f"%(num_T0/float(num_found), num_T5/float(num_found), num_T10/float(num_found), num_T15/float(num_found), insufficient_length/float(num_found))
		#	print '\n'
		for idx, val in enumerate(_sample_print_array):
		#	print "Outputting %s with %i number of reads!"%(_barcodes_names[idx], len(val))
			_file_name = _barcodes_names[idx] + "_completed.fastq"
			_output_handle = open(_file_name, 'a')
			SeqIO.write(val, _output_handle, "fastq")
		_sample_print_array = []
		for idx, val in enumerate(_barcodes):
			_sample_print_array.append([])

		print '=== TIME TAKEN FOR FILLER MATCH (every 500000 records): %s sec ==='%(time.time()-start)

print '%s\t%s\t%s\t%s\t%s\t%s\t%s'%('# records','# found','# ambiguous','# unfound','% found', '% ambiguous','% unfound')
print "%i\t%i\t%i\t%i\t%0.03f\t%0.03f\t%0.03f"%(num_records, num_found, num_ambiguous, num_unfound, 100*num_found/float(num_records), 100*num_ambiguous/float(num_records), 100*num_unfound/float(num_ambiguous))
print "\nBarcode Hamming Array Stats"
for idx, val in enumerate(barcode_hamming_array):
	try:
		print "%i\t%i\t%0.03f"%(idx, val, val/float(num_found))
	except ZeroDivisionError:
		#pass
		print "%i\t%i\t%0.03f"%(idx, val, 0)

print "\nBarcode Found Index Array Stats"
for idx, val in enumerate(barcode_found_index_array):
	try:
		print "%i\t%i\t%0.03f"%(idx, val, val/float(num_found))
	except ZeroDivisionError:
		#pass
		print "%i\t%i\t%0.03f"%(idx, val, 0)

print "\nBarcode Counter Array Stats"
for idx, val in enumerate(barcode_counter_array):
	try:
		print "%i\t%i\t%0.03f"%(idx, val, val/float(num_found))
	except ZeroDivisionError:
		#pass
		print "%i\t%i\t%0.03f"%(idx, val, 0)

print "\nUMI N distribution"
for idx, val in enumerate(N_distibution):
	try:
		print "%i\t%i\t%0.03f"%(idx, val, val/float(sum(N_distibution)))
	except ZeroDivisionError:
		#pass
		print "%i\t%i\t%0.03f"%(idx, val, 0)

print "\nFiller Hamming"
for idx, val in enumerate(filler_hamming_distribution):
	try:
		print "%i\t%i\t%0.03f"%(idx, val, val/float(sum(filler_hamming_distribution)))
	except ZeroDivisionError:
		#pass
		print "%i\t%i\t%0.03f"%(idx, val, 0)

print "\nUMI error checking"
print '%s\t%s\t%s\t%s\t%s'%('% T0/found','% T5/found', '% T10/found','% T15/found','% insufficient length/found')
try:
	print "%0.03f\t%0.03f\t%0.03f\t%0.03f\t%0.05f"%(num_T0/float(num_found), num_T5/float(num_found), num_T10/float(num_found), num_T15/float(num_found), insufficient_length/float(num_found))
except ZeroDivisionError:
	print "%0.03f\t%0.03f\t%0.03f\t%0.03f\t%0.05f"%(0, 0, 0, 0, 0)
print '\n'

for idx, val in enumerate(_sample_print_array):
	print "Outputting %s with %i number of reads!"%(_barcodes_names[idx], len(val))
	_file_name = _barcodes_names[idx] + "_completed.fastq"
	_output_handle = open(_file_name, 'a')
	SeqIO.write(val, _output_handle, "fastq")

print '\n=== TIME TAKEN FOR FILLER MATCH: %s sec ==='%(time.time()-start)

try:
    sys.stdout.close()
except:
    pass
try:
    sys.stderr.close()
except:
    pass
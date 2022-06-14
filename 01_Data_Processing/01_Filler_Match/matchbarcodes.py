# match barcode ID with corresponding NT seq
import os, sys, re, time, glob
import argparse 

common_seq = 'CAGTGGTATCAACGCAGAG'

def main():
#    parser = argparse.ArgumentParser(description = 'Queries complete barcode list for a set of barcodes we want to extract and then print them out.')
#    parser.add_argument('-in', '--input', type = str, required=True, help = 'Set of query barcodes')
#    parser.add_argument('-in2', '--input2', type = str, required=True, help = 'Set of all barcodes')
#    parser.add_argument('-out', '--output', type = str, required=True, help = 'Returns 8 nucleotide barcode sequences for query barcodes in .txt file')
#    args = parser.parse_args()
    dir = sys.argv[1]
    targetbarcodes = glob.glob(dir+'target*txt')[0]
    completebarcodes = glob.glob(dir+'complete*txt')[0]
    barcode_set = set()
    global common_seq
    
    for line in open(targetbarcodes, 'r'):
        line.rstrip()
        line = line.upper()
        split_line = list(line)
        split_line[3] = '1'
        line = ''.join(split_line)
        line = line.rstrip()
 #       print line
        barcode_set.add(line)

    outfile = open('new_barcodes.txt', 'w')
    
    matched = 0
    for line in open(completebarcodes, 'r'):
        line.rstrip()
        split_line = line.split("\t")
        barcode = split_line[0]
        barcode = barcode.rstrip()
        dna = split_line[1]
        dna = dna.rstrip()
        #print barcode
        if barcode in barcode_set:
            matched += 1
            common_index = dna.index(common_seq)
            barcode_index = common_index-8
            new_barcodes = dna[barcode_index:common_index]
            outfile.write("%s\t%s\n"%(barcode, new_barcodes))
    outfile.close()

    print '# target barcodes: %s'%(sum(1 for line in open(targetbarcodes, 'r')))
    print '# all barcodes: %s'%(sum(1 for line in open(completebarcodes, 'r')))
    print '# matched barcodes: %s'%(matched)
    os.system("cat %s"%('new_barcodes.txt'))

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))

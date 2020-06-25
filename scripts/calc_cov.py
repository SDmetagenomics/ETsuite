import pysam
import sys
import pandas as pd
import argparse
from Bio import SeqIO
def read_bam(scaf2bin, bam, fasta, min_ani, min_mapq):
	f = open(scaf2bin)
	scaf2bin = {}
	lengths = {}
	for record in SeqIO.parse(fasta, 'fasta'):
		lengths[record.id] = len(record.seq)
	for line in f.readlines():
		scaf2bin[line.split()[0]] = line.split()[1]
	coverage = {}
	samfile = pysam.AlignmentFile(bam)
	print("SCAFFOLD\tGENOME\tCOVERAGE\tREADS")
	for contig in scaf2bin:
		read_lengths = []
		count = 0
		for read in samfile.fetch(contig):
			if read.get_reference_positions() != []:
				ani = 1-float(read.get_tag("NM")) / float(len(read.get_reference_positions())) 
				if read.mapping_quality >= min_mapq and ani >= min_ani:
					count +=1
					read_lengths.append(len(read.get_reference_positions()))
		if len(read_lengths) > 0:
			print(contig + "\t" + scaf2bin[contig] + "\t" + str(float(count)*sum(read_lengths) / len(read_lengths) / lengths[contig]) + "\t" + str(count))
		else:
			print(contig + "\t" + scaf2bin[contig] + "\t" + "0")
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calculate filtered coverages for scaffolds",
             formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # Required positional arguments
    parser.add_argument("scaf2bin", help="scaffold to bin file, tab separated.")    
    parser.add_argument("bam", help="Sorted, indexed .bam file.")
    parser.add_argument("fasta", help="FASTA file.")
    parser.add_argument("-c", "--min_ani", action="store", default=0.95, type=float, \
        help='Minimum required percent identity of read pairs to consensus. (>=)')
    parser.add_argument("--min_mapq", action="store", default=2, type=int,\
        help='Minimum mapq score of EITHER read in a pair to use that pair. (>=) ')
    args = parser.parse_args()
    run = read_bam(args.scaf2bin, args.bam, args.fasta, args.min_ani, args.min_mapq)

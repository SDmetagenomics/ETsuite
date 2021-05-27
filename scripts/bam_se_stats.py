# bam_se_stats v0.2

import pysam
import sys
import pandas as pd
import argparse
from collections import defaultdict
def read_bam(scaf2bin, bam):
    f = open(scaf2bin)
    scaf2bin = {}
    for line in f.readlines():
            scaf2bin[line.split()[0]] = line.split()[1].strip()
    f.close()
    samfile = pysam.AlignmentFile(bam)
    total_read_count = 0
    read_info = {}
    found_both = set()
    reads = defaultdict(list)
    for scaf in scaf2bin:
            for read in samfile.fetch(scaf):
                    total_read_count += 1
                    if read.get_reference_positions() != []:
                            if read.is_reverse:
                                strand = '-'
                            else:
                                strand = '+'
                            if read.query_name not in read_info:
                                    qseq = read.query_sequence
                                    refseq = read.get_reference_sequence()
                                    good_pos = 0
                                    true_nm = 0
                                    first_three = 0
                                    clipped = 'False'
                                    if len(qseq) == len(refseq):
                                        read_len_qual = len(read.query_alignment_qualities)
                                        for pos in range(0,read_len_qual):
                                            if read.query_alignment_qualities[pos] > 20:
                                                good_pos += 1
                                                if qseq[pos] != refseq[pos]:
                                                    true_nm += 1
                                            if strand == '+' and pos <= 2:
                                                if qseq[pos] != refseq[pos]:
                                                    first_three += 1
                                            elif strand == '-' and pos >= read_len_qual-3:
                                                if qseq[pos] != refseq[pos]:
                                                    first_three += 1
                                    else:
                                        clipped = 'True'
                                    reads['Read'].append(read.query_name)
                                    reads['SCAF1'].append(scaf)
                                    reads['GENOME1'].append(scaf2bin[scaf])
                                    reads['MAPQ1'].append(read.mapping_quality)
                                    reads['NM1'].append(read.get_tag('NM'))
                                    reads['LEN1'].append(len(read.get_reference_positions()))
                                    reads['START1'].append(read.get_reference_positions()[0])
                                    reads['END1'].append(read.get_reference_positions()[-1])
                                    reads['STRAND1'].append(strand)
                                    reads['QUAL1'].append(good_pos)
                                    reads['QUALNM1'].append(true_nm)
                                    reads['CLIPPED1'].append(clipped)
                                    reads['R1_START_NM'].append(first_three)
    pd.DataFrame(reads).to_csv(sys.stdout, index=False, sep="\t")
    return True
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Read Paired End BAM data",
             formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # Required positional arguments
    parser.add_argument("scaf2bin", help="scaffold to bin file, tab separated.")
    parser.add_argument("bam", help="Sorted, indexed .bam file.")
    args = parser.parse_args()
    run = read_bam(args.scaf2bin, args.bam)

import pysam
import sys
import pandas as pd
import argparse


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

    from collections import defaultdict

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
                                    read_info[read.query_name] = {"scaf": scaf, 'is_read1': read.is_read1, 'strand': strand, "mapq": read.mapping_quality, "NM": read.get_tag("NM"), 
                                    "len": len(read.get_reference_positions()), "end": read.get_reference_positions()[-1], "start": read.get_reference_positions()[0]}
                            else:
                                    found_both.add(read.query_name)
                                    if read.is_read1:
                                            reads['Read'].append(read.query_name)
                                            reads['SCAF1'].append(scaf)
                                            reads['GENOME1'].append(scaf2bin[scaf])
                                            reads['MAPQ1'].append(read.mapping_quality)
                                            reads['NM1'].append(read.get_tag('NM'))
                                            reads['LEN1'].append(len(read.get_reference_positions()))
                                            reads['START1'].append(read.get_reference_positions()[0])
                                            reads['END1'].append(read.get_reference_positions()[-1])
                                            reads['STRAND1'].append(strand)

                                            reads['SCAF2'].append(read_info[read.query_name]['scaf'])
                                            reads['GENOME2'].append(scaf2bin[read_info[read.query_name]['scaf']])
                                            reads['MAPQ2'].append(read_info[read.query_name]['mapq'])
                                            reads['NM2'].append(read_info[read.query_name]['NM'])
                                            reads['LEN2'].append(read_info[read.query_name]['len'])
                                            reads['START2'].append(read_info[read.query_name]['start'])
                                            reads['END2'].append(read_info[read.query_name]['end'])
                                            reads['STRAND2'].append(read_info[read.query_name]['strand'])

                                    elif read.is_read2:
                                            reads['Read'].append(read.query_name)
                                            reads['SCAF2'].append(scaf)
                                            reads['GENOME2'].append(scaf2bin[scaf])
                                            reads['MAPQ2'].append(read.mapping_quality)
                                            reads['NM2'].append(read.get_tag('NM'))
                                            reads['LEN2'].append(len(read.get_reference_positions()))
                                            reads['START2'].append(read.get_reference_positions()[0])
                                            reads['END2'].append(read.get_reference_positions()[-1])
                                            reads['STRAND2'].append(strand)

                                            reads['SCAF1'].append(read_info[read.query_name]['scaf'])
                                            reads['GENOME1'].append(scaf2bin[read_info[read.query_name]['scaf']])
                                            reads['MAPQ1'].append(read_info[read.query_name]['mapq'])
                                            reads['NM1'].append(read_info[read.query_name]['NM'])
                                            reads['LEN1'].append(read_info[read.query_name]['len'])
                                            reads['START1'].append(read_info[read.query_name]['start'])
                                            reads['END1'].append(read_info[read.query_name]['end'])
                                            reads['STRAND1'].append(read_info[read.query_name]['strand'])



    for read in read_info:
        if read not in found_both:
            r = read_info[read]
            if read_info[read]['is_read1']:
                reads['Read'].append(read)
                reads['SCAF1'].append(read_info[read]['scaf'])
                reads['GENOME1'].append(scaf2bin[read_info[read]['scaf']])
                reads['MAPQ1'].append(read_info[read]['mapq'])
                reads['NM1'].append(read_info[read]['NM'])
                reads['LEN1'].append(read_info[read]['len'])
                reads['START1'].append(read_info[read]['start'])
                reads['END1'].append(read_info[read]['end'])
                reads['STRAND1'].append(read_info[read]['strand'])

                reads['SCAF2'].append('NA')
                reads['GENOME2'].append('NA')
                reads['MAPQ2'].append('NA')
                reads['NM2'].append('NA')
                reads['LEN2'].append('NA')
                reads['START2'].append('NA')
                reads['END2'].append('NA')
                reads['STRAND2'].append('NA')
            else:
                reads['Read'].append(read)
                reads['SCAF2'].append(read_info[read]['scaf'])
                reads['GENOME2'].append(scaf2bin[read_info[read]['scaf']])
                reads['MAPQ2'].append(read_info[read]['mapq'])
                reads['NM2'].append(read_info[read]['NM'])
                reads['LEN2'].append(read_info[read]['len'])
                reads['START2'].append(read_info[read]['start'])
                reads['END2'].append(read_info[read]['end'])
                reads['STRAND2'].append(read_info[read]['strand'])

                reads['SCAF1'].append('NA')
                reads['GENOME1'].append('NA')
                reads['MAPQ1'].append('NA')
                reads['NM1'].append('NA')
                reads['LEN1'].append('NA')
                reads['START1'].append('NA')
                reads['END1'].append('NA')
                reads['STRAND1'].append('NA')

    print(pd.DataFrame(reads).to_csv(sys.stdout, index=False, sep="\t"))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Read Paired End BAM data",
             formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Required positional arguments
    parser.add_argument("scaf2bin", help="scaffold to bin file, tab separated.")
    
    parser.add_argument("bam", help="Sorted, indexed .bam file.")
    
    args = parser.parse_args()
    read_bam(args.scaf2bin, args.bam)
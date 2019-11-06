
import pysam
import sys

f = open(sys.argv[1])
scaf2bin = {}
for line in f.readlines():
	scaf2bin[line.split()[0]] = line.split()[1].strip()
f.close()


samfile = pysam.AlignmentFile(sys.argv[2])
total_read_count = 0

for scaf in scaf2bin:
	for read in samfile.fetch(scaf):
		total_read_count += 1
		if read.get_reference_positions() != []:
			if read.is_reverse:
				strand = "-"
			else:
				strand = "+"
			print(str(read.query_name) + "\t" + scaf2bin[scaf] + "\t" + str(read.mapping_quality) + "\t" + str(read.get_tag('NM')) + "\t" + str(len(read.get_reference_positions())) + "\t" + str(read.get_reference_positions()[0]) + "\t" + str(read.get_reference_positions()[-1]) + "\t" + str(strand))

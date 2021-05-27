# bam_se_stats v0.1
import pysam
import sys

f = open(sys.argv[1])
scaf2bin = {}
for line in f.readlines():
	scaf2bin[line.split()[0]] = line.split()[1].strip()
f.close()


samfile = pysam.AlignmentFile(sys.argv[2])
total_read_count = 0


print("Read\tGENOME\tSCAF\tMAPQ\tNM\tLEN\tSTART\tEND\tSTRAND\tQUAL\tQUALNM\tCLIPPED")
for scaf in scaf2bin:
	for read in samfile.fetch(scaf):
		total_read_count += 1
		if read.get_reference_positions() != []:

			qseq = read.query_sequence
			refseq = read.get_reference_sequence()
			good_pos = 0
			true_nm = 0
			clipped = 'False'
			if len(qseq) == len(refseq):
				for pos in range(0,len(read.query_alignment_qualities)):
					if read.query_alignment_qualities[pos] > 20:
						good_pos += 1
						if qseq[pos] != refseq[pos]:
							true_nm += 1
			else:
				clipped = 'True'

			if read.is_reverse:
				strand = "-"
			else:
				strand = "+"
			print(str(read.query_name) + "\t" + scaf2bin[scaf] + "\t" + scaf + "\t" + str(read.mapping_quality) + "\t" + str(read.get_tag('NM')) + "\t" + str(len(read.get_reference_positions())) + "\t" + str(read.get_reference_positions()[0]) + "\t" + str(read.get_reference_positions()[-1]) + "\t" + str(strand) + "\t" + str(good_pos) + "\t" + str(true_nm) + "\t" + str(clipped))

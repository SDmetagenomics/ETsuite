import pysam
import sys

f = open(sys.argv[1])
scaf2bin = {}
for line in f.readlines():
        scaf2bin[line.split()[0]] = line.split()[1].strip()
f.close()


samfile = pysam.AlignmentFile(sys.argv[2])
total_read_count = 0
read_info = {}
found_both = set()
print("Read\tGENOME1\tSCAF1\tGENOME2\tSCAF2\tMAPQ1\tNM1\tLEN1\tSTART1\tEND1\tSTRAND1\tMAPQ2\tNM2\tLEN2\tSTART2\tEND2\tSTRAND2")
for scaf in scaf2bin:
        for read in samfile.fetch(scaf):
                total_read_count += 1
                if read.get_reference_positions() != []:
                        if read.is_reverse:
                            strand = '-'
                        else:
                            strand = '+'

                        if read.query_name not in read_info:
                                read_info[read.query_name] = {"scaf": scaf, 'is_read1': read.is_read1, 'strand': strand, "mapq": read.mapping_quality, "NM": read.get_tag("NM"), "len1": len(read.get_reference_positions()), "end": read.get_reference_positions()[-1], "start": read.get_reference_positions()[0]}
                        else:
                                found_both.add(read.query_name)
                                if read.is_read1:
                                        print(str(read.query_name) + "\t" + scaf2bin[scaf] + "\t" + scaf + "\t" + scaf2bin[read_info[read.query_name]['scaf']] + "\t" + read_info[read.query_name]['scaf'] + "\t" + str(read.mapping_quality) + "\t" + str(read.get_tag('NM')) + 
                                                "\t" + str(len(read.get_reference_positions())) + "\t" + str(read.get_reference_positions()[0]) + "\t" +  str(read.get_reference_positions()[-1]) + "\t" + strand +
                                                "\t" + str(read_info[read.query_name]['mapq']) + "\t" + 
                                                str(read_info[read.query_name]['NM']) + "\t" + str(read_info[read.query_name]['len1']) + 
                                                "\t" + str(read_info[read.query_name]['start']) + "\t" + str(read_info[read.query_name]['end']) + "\t" + read_info[read.query_name]['strand'])
                                elif read.is_read2:
                                        print(str(read.query_name) + "\t" + scaf2bin[scaf] + "\t" + scaf + "\t" + scaf2bin[read_info[read.query_name]['scaf']] + "\t" + read_info[read.query_name]['scaf'] + "\t" + str(read_info[read.query_name]['mapq']) + "\t" + str(read_info[read.query_name]['NM']) + "\t" + 
                                                str(read_info[read.query_name]['len1']) + "\t" + str(read_info[read.query_name]['start']) + "\t" + str(read_info[read.query_name]['end']) + "\t" + read_info[read.query_name]['strand'] +
                                                 "\t" + str(read.mapping_quality) + "\t" + str(read.get_tag('NM')) + 
                                                "\t" + str(len(read.get_reference_positions())) + "\t" + str(read.get_reference_positions()[0]) + "\t" +  str(read.get_reference_positions()[-1]) + "\t" + strand)
for read in read_info:
    if read not in found_both:
        r = read_info[read]
        if read_info[read]['is_read1']:
            print(str(read) + "\t" + scaf2bin[r['scaf']] + "\t" + r['scaf'] + "\t" + "NA" + "\t" + "NA" + "\t" + str(r['mapq']) + "\t" + str(r['NM']) + 
                    "\t" + str(r['len1']) + "\t" + str(r['start']) + "\t" +  str(r['end']) + "\t" + r['strand'] +
                    "\t" + "NA" + "\t" + 
                    "NA" + "\t" + "NA" + 
                    "\t" + "NA" + "\t" + "NA" + "\t" + "NA")
        else:
            print(str(read) + "\t" + scaf2bin[r['scaf']] + "\t" + r['scaf'] + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + 
                    "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" +
                     "\t" + str(r['mapq']) + "\t" + str(r['mapq']) + 
                    "\t" + str(r['len1']) + "\t" + str(r['start']) + "\t" +  str(r['end']) + "\t" + r['strand'])

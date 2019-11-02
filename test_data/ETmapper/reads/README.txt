> Data Origin:
This test data is from ET-seq Experiment conducted on 8/26/2019. The Data contains 5
samples: 3 prepared with ETseq2 protocol and 2 prepared with ETseq 4 protocol, also
each uses a different 

> Data Overview

Sample	Total_Reads	Illumina_Adap	Total_Model
JD_ZZ_2ndcycle_1	139524	1585	126691
JD_ZZ_2ndcycle_2	223295	5110	199997
JD_ZZ_2ndcycle_3	172263	16248	154150
JD_ZZ_2ndcycle_4	166824	7898	151325
JD_ZZ_2ndcycle_5	243134	54893	200575

next thing look at X option



> Run read trimming on fwd and quality filter both
cutadapt -a file:/Users/Spencer/Dropbox/Banfield_Lab_Files/Projects/mCAFE/mCAFE_Project/ET_Mapper/db/ETseq_2_4_adap.fa -j 6 -O 5 -q 20 -o JD_ZZ_2ndcycle_5_S5_L001_R1_001.fastq.trim -p JD_ZZ_2ndcycle_5_S5_L001_R2_001.fastq.trim JD_ZZ_2ndcycle_5_S5_L001_R1_001.fastq JD_ZZ_2ndcycle_5_S5_L001_R2_001.fastq

=== Summary ===

Total read pairs processed:            243,134
  Read 1 with adapter:                  15,825 (6.5%)
  Read 2 with adapter:                       0 (0.0%)
Pairs written (passing filters):       243,134 (100.0%)

Total basepairs processed:    73,426,468 bp
  Read 1:    36,713,234 bp
  Read 2:    36,713,234 bp
Quality-trimmed:               1,480,070 bp (2.0%)
  Read 1:       329,171 bp
  Read 2:     1,150,899 bp
Total written (filtered):     71,719,182 bp (97.7%)
  Read 1:    36,156,847 bp
  Read 2:    35,562,335 bp





> Testing with Linked adapter 


 system(paste0("cutadapt -a file:",ad," -A file:",ad," -G file:",ad, # specify adapter types and file
                  " -j 4 -O ",al," -q ",qs, # specify trimming params
                  " -o ",batch_file[i,3],".trim", # fwd output
                  " -p ",batch_file[i,4],".trim", # rev output
                  " ",rd,batch_file[i,3], # fwd read
                  " ",rd,batch_file[i,4], # rev read
                  " > ",batch_file[i,1],".trim.log")) # log files

=== Summary ===

Total read pairs processed:            139,524
  Read 1 with adapter:                   2,137 (1.5%)
  Read 2 with adapter:                 128,004 (91.7%)
Pairs written (passing filters):       139,524 (100.0%)

Total basepairs processed:    42,136,248 bp
  Read 1:    21,068,124 bp
  Read 2:    21,068,124 bp
Quality-trimmed:                 355,729 bp (0.8%)
  Read 1:        50,773 bp
  Read 2:       304,956 bp
Total written (filtered):     38,861,807 bp (92.2%)
  Read 1:    20,943,884 bp
  Read 2:    17,917,923 bp

=== Second read: Adapter 3 ===

Sequence: TGCTGAACCGCTCTTCCGATCT...AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT; Type: linked; Length: 22+58; 5' trimmed: 127777 times; 3' trimmed: 999 times

No. of allowed errors:
0-9 bp: 0; 10-19 bp: 1; 20-22 bp: 2

No. of allowed errors:
0-9 bp: 0; 10-19 bp: 1; 20-29 bp: 2; 30-39 bp: 3; 40-49 bp: 4; 50-58 bp: 5

Overview of removed sequences at 5' end
length  count   expect  max.err error counts
8       1       2.1     0       1
13      2       0.0     1       0 2
15      1       0.0     1       0 1
20      23      0.0     2       0 0 23
21      964     0.0     2       59 514 391
22      126622  0.0     2       95991 26526 4105
23      129     0.0     2       31 90 8
24      35      0.0     2       7 1 27


#### NO -G this time
cutadapt -a file:/Users/Spencer/Dropbox/Banfield_Lab_Files/Projects/mCAFE/mCAFE_Project/ET_Mapper/db/ETseq_2_4_adap.fa -A file:/Users/Spencer/Dropbox/Banfield_Lab_Files/Projects/mCAFE/mCAFE_Project/ET_Mapper/db/ETseq_2_4_adap.fa -j 6 -O 5 -q 20 -o JD_ZZ_2ndcycle_5_S5_L001_R1_001.fastq.trim -p JD_ZZ_2ndcycle_5_S5_L001_R2_001.fastq.trim JD_ZZ_2ndcycle_5_S5_L001_R1_001.fastq JD_ZZ_2ndcycle_5_S5_L001_R2_001.fastq


system(paste0("cutadapt -a file:",ad," -A file:",ad, # specify adapter types and file
                  " -j 4 -O ",al," -q ",qs, # specify trimming params
                  " -o ",batch_file[i,3],".trim", # fwd output
                  " -p ",batch_file[i,4],".trim", # rev output
                  " ",rd,batch_file[i,3], # fwd read
                  " ",rd,batch_file[i,4], # rev read
                  " > ",batch_file[i,1],".trim.log")) # log files

=== Summary ===

Total read pairs processed:            139,524
  Read 1 with adapter:                   2,137 (1.5%)
  Read 2 with adapter:                 127,974 (91.7%)
Pairs written (passing filters):       139,524 (100.0%)

Total basepairs processed:    42,136,248 bp
  Read 1:    21,068,124 bp
  Read 2:    21,068,124 bp
Quality-trimmed:                 355,729 bp (0.8%)
  Read 1:        50,773 bp
  Read 2:       304,956 bp
Total written (filtered):     38,861,966 bp (92.2%)
  Read 1:    20,943,884 bp
  Read 2:    17,918,082 bp




### NO -G this time and Anchored linked adapter using X
system(paste0("cutadapt -a file:",ad," -A file:",ad, # specify adapter types and file
                  " -j 4 -O ",al," -q ",qs, # specify trimming params
                  " -o ",batch_file[i,3],".trim", # fwd output
                  " -p ",batch_file[i,4],".trim", # rev output
                  " ",rd,batch_file[i,3], # fwd read
                  " ",rd,batch_file[i,4], # rev read
                  " > ",batch_file[i,1],".trim.log")) # log files


=== Summary ===

Total read pairs processed:            139,524
  Read 1 with adapter:                   2,137 (1.5%)
  Read 2 with adapter:                 127,935 (91.7%)
Pairs written (passing filters):       139,524 (100.0%)

Total basepairs processed:    42,136,248 bp
  Read 1:    21,068,124 bp
  Read 2:    21,068,124 bp
Quality-trimmed:                 355,729 bp (0.8%)
  Read 1:        50,773 bp
  Read 2:       304,956 bp
Total written (filtered):     38,862,914 bp (92.2%)
  Read 1:    20,943,884 bp
  Read 2:    17,919,030 bp










### Testing model trimming from Fwd Read

cutadapt -g file:/Users/Spencer/Dropbox/Banfield_Lab_Files/Projects/mCAFE/mCAFE_Project/ET_Mapper/db/models.fa -e 0.02 -O 25 --discard-untrimmed --info-file tmp.txt -o JD_ZZ_2ndcycle_5_S5_L001_R1_001.fastq.clean -p JD_ZZ_2ndcycle_5_S5_L001_R2_001.fastq.clean JD_ZZ_2ndcycle_5_S5_L001_R1_001.fastq.trim JD_ZZ_2ndcycle_5_S5_L001_R2_001.fastq.trim

This is cutadapt 2.5 with Python 3.6.0
Command line parameters: -g file:/Users/Spencer/Dropbox/Banfield_Lab_Files/Projects/mCAFE/mCAFE_Project/ET_Mapper/db/models.fa -e 0.02 -O 25 --discard-untrimmed --info-file tmp.txt -o JD_ZZ_2ndcycle_5_S5_L001_R1_001.fastq.clean -p JD_ZZ_2ndcycle_5_S5_L001_R2_001.fastq.clean JD_ZZ_2ndcycle_5_S5_L001_R1_001.fastq.trim JD_ZZ_2ndcycle_5_S5_L001_R2_001.fastq.trim
Processing reads on 1 core in paired-end mode ...
[--------->8 ] 00:00:25       243,134 reads  @    106.1 µs/read;   0.57 M reads/minute
Finished in 25.83 s (106 us/read; 0.56 M reads/minute).

=== Summary ===

Total read pairs processed:            243,134
  Read 1 with adapter:                 193,455 (79.6%)
  Read 2 with adapter:                       0 (0.0%)
Pairs written (passing filters):       193,455 (79.6%)

Total basepairs processed:    68,918,195 bp
  Read 1:    33,355,860 bp
  Read 2:    35,562,335 bp
Total written (filtered):     38,017,440 bp (55.2%)
  Read 1:     9,275,627 bp
  Read 2:    28,741,813 bp

=== First read: Adapter model_magic_Tn5 ===

Sequence: NNNNNNGATGTCCACGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTCGACCAGCGGCCGGCCGGTTGAGATGTGTATAAGAGACAG; Type: regular 5'; Length: 98; Trimmed: 159455 times.

No. of allowed errors:
0-49 bp: 0; 50-72 bp: 1

Overview of removed sequences
length	count	expect	max.err	error counts
87	20	0.0	1	20
92	19	0.0	1	18 1
93	105	0.0	1	103 2
95	4	0.0	1	3 1
96	25	0.0	1	20 5
97	1355	0.0	1	253 1102
98	157464	0.0	1	153724 3740
99	456	0.0	1	168 288
100	4	0.0	1	4
102	2	0.0	1	2
103	1	0.0	1	1


=== First read: Adapter model_magic_mariner ===

Sequence: NNNNNNGATGTCCACGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTCGACCAGCGGCCGGCCAGACCGGGGACTTATCAGCCAACCTGT; Type: regular 5'; Length: 101; Trimmed: 34000 times.

No. of allowed errors:
0-49 bp: 0; 50-75 bp: 1

Overview of removed sequences
length	count	expect	max.err	error counts
90	4	0.0	1	4
95	5	0.0	1	5
96	29	0.0	1	29
97	1	0.0	1	0 1
98	1	0.0	1	0 1
99	5	0.0	1	4 1
100	256	0.0	2	53 203
101	33604	0.0	2	32705 899
102	94	0.0	2	28 66
105	1	0.0	2	1


=== First read: Adapter model_revcomp_mariner ===

Sequence: ACAGGTTGGCTGATAAGTCCCCGGTCTGGCCGGCCGCTGGTCGACCTGCAGCGTACGNNNNNNNNNNNNNNNNNNNNAGAGACCTCGTGGACATCNNNNNN; Type: regular 5'; Length: 101; Trimmed: 0 times.

=== First read: Adapter model_revcomp_Tn5 ===

Sequence: CTGTCTCTTATACACATCTCAACCGGCCGGCCGCTGGTCGACCTGCAGCGTACGNNNNNNNNNNNNNNNNNNNNAGAGACCTCGTGGACATCNNNNNN; Type: regular 5'; Length: 98; Trimmed: 0 times.



	- When the reverse complement of the model sequences were included in the model file this had no major effect on trimming the
	forward reads in this trial. NO REVERSE READS GET TRIMMED.Ideally we want to have both forward and reverse versions of these
	sequences in the same file so the file can be used for trimming both forward and reverse reads without the need for an additional file.
	



### Testing model trimming from Rev Read

cutadapt -A file:/Users/Spencer/Dropbox/Banfield_Lab_Files/Projects/mCAFE/mCAFE_Project/ET_Mapper/db/models.fa -e 0.02 -O 25 -o JD_ZZ_2ndcycle_5_S5_L001_R1_001.fastq.clean -p JD_ZZ_2ndcycle_5_S5_L001_R2_001.fastq.clean JD_ZZ_2ndcycle_5_S5_L001_R1_001.fastq.trim JD_ZZ_2ndcycle_5_S5_L001_R2_001.fastq.trim
This is cutadapt 2.5 with Python 3.6.0
Command line parameters: -A file:/Users/Spencer/Dropbox/Banfield_Lab_Files/Projects/mCAFE/mCAFE_Project/ET_Mapper/db/models.fa -e 0.02 -O 25 -o JD_ZZ_2ndcycle_5_S5_L001_R1_001.fastq.clean -p JD_ZZ_2ndcycle_5_S5_L001_R2_001.fastq.clean JD_ZZ_2ndcycle_5_S5_L001_R1_001.fastq.trim JD_ZZ_2ndcycle_5_S5_L001_R2_001.fastq.trim
Processing reads on 1 core in paired-end mode ...
[    8<------] 00:00:10       243,134 reads  @     45.2 µs/read;   1.33 M reads/minute
Finished in 11.02 s (45 us/read; 1.32 M reads/minute).

=== Summary ===

Total read pairs processed:            243,134
  Read 1 with adapter:                       0 (0.0%)
  Read 2 with adapter:                  98,738 (40.6%)
Pairs written (passing filters):       243,134 (100.0%)

Total basepairs processed:    68,918,195 bp
  Read 1:    33,355,860 bp
  Read 2:    35,562,335 bp
Total written (filtered):     61,641,031 bp (89.4%)
  Read 1:    33,355,860 bp
  Read 2:    28,285,171 bp

=== Second read: Adapter model_magic_Tn5 ===

Sequence: NNNNNNGATGTCCACGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTCGACCAGCGGCCGGCCGGTTGAGATGTGTATAAGAGACAG; Type: regular 3'; Length: 98; Trimmed: 0 times.

=== Second read: Adapter model_magic_mariner ===

Sequence: NNNNNNGATGTCCACGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTCGACCAGCGGCCGGCCAGACCGGGGACTTATCAGCCAACCTGT; Type: regular 3'; Length: 101; Trimmed: 0 times.

=== Second read: Adapter model_revcomp_mariner ===

Sequence: ACAGGTTGGCTGATAAGTCCCCGGTCTGGCCGGCCGCTGGTCGACCTGCAGCGTACGNNNNNNNNNNNNNNNNNNNNAGAGACCTCGTGGACATCNNNNNN; Type: regular 3'; Length: 101; Trimmed: 7282 times.

No. of allowed errors:
0-49 bp: 0; 50-75 bp: 1

Bases preceding removed adapters:
  A: 98.9%
  C: 0.0%
  G: 0.0%
  T: 0.0%
  none/other: 1.0%
WARNING:
    The adapter is preceded by "A" extremely often.
    The provided adapter sequence could be incomplete at its 3' end.

Overview of removed sequences
length	count	expect	max.err	error counts
25	4	0.0	0	4
26	6	0.0	0	6
27	13	0.0	0	13
28	6	0.0	0	6
29	3	0.0	0	3
30	11	0.0	0	11
31	24	0.0	0	24
32	12	0.0	0	12
33	10	0.0	0	10
34	6	0.0	0	6
35	6	0.0	0	6
36	15	0.0	0	15
37	15	0.0	0	15
38	10	0.0	0	10
39	9	0.0	0	9
40	5	0.0	0	5
41	10	0.0	0	10
42	13	0.0	0	13
43	11	0.0	0	11
44	22	0.0	0	22
45	15	0.0	0	15
46	10	0.0	0	10
47	13	0.0	0	13
48	11	0.0	0	11
49	1	0.0	0	1
50	14	0.0	1	14
51	5	0.0	1	5
52	8	0.0	1	7 1
53	4	0.0	1	4
54	13	0.0	1	13
55	4	0.0	1	4
56	6	0.0	1	6
57	5	0.0	1	5
58	17	0.0	1	16 1
59	9	0.0	1	8 1
60	8	0.0	1	8
61	2	0.0	1	2
62	7	0.0	1	7
63	8	0.0	1	8
64	3	0.0	1	3
65	4	0.0	1	4
66	3	0.0	1	3
67	3	0.0	1	3
68	3	0.0	1	3
69	7	0.0	1	6 1
70	9	0.0	1	9
71	7	0.0	1	7
72	3	0.0	1	3
73	10	0.0	1	10
74	35	0.0	1	34 1
75	7	0.0	1	7
76	12	0.0	1	9 3
77	17	0.0	1	16 1
78	25	0.0	1	7 18
79	9	0.0	1	7 2
80	3	0.0	1	3
81	5	0.0	1	5
82	5	0.0	1	4 1
83	5	0.0	1	5
84	5	0.0	1	5
85	11	0.0	1	9 2
86	103	0.0	1	100 3
87	80	0.0	1	79 1
88	71	0.0	1	67 4
89	226	0.0	1	220 6
90	58	0.0	1	55 3
91	81	0.0	1	79 2
92	122	0.0	1	118 4
93	215	0.0	1	210 5
94	177	0.0	1	174 3
95	84	0.0	1	78 6
96	108	0.0	1	106 2
97	43	0.0	1	41 2
98	60	0.0	1	58 2
99	101	0.0	1	100 1
100	318	0.0	2	304 14
101	78	0.0	2	73 5
102	101	0.0	2	99 2
103	458	0.0	2	446 12
104	126	0.0	2	121 5
105	181	0.0	2	175 6
106	52	0.0	2	48 4
107	98	0.0	2	94 4
108	47	0.0	2	45 2
109	41	0.0	2	40 1
110	49	0.0	2	47 2
111	110	0.0	2	107 3
112	40	0.0	2	34 6
113	310	0.0	2	296 14
114	71	0.0	2	67 4
115	81	0.0	2	78 3
116	158	0.0	2	155 3
117	79	0.0	2	77 2
118	81	0.0	2	79 2
119	48	0.0	2	48
120	153	0.0	2	150 3
121	91	0.0	2	89 2
122	141	0.0	2	136 5
123	98	0.0	2	98
124	91	0.0	2	86 5
125	95	0.0	2	92 3
126	77	0.0	2	77
127	23	0.0	2	23
128	37	0.0	2	35 2
129	89	0.0	2	86 3
130	65	0.0	2	64 1
131	129	0.0	2	122 7
132	24	0.0	2	21 3
133	99	0.0	2	97 2
134	30	0.0	2	30
135	37	0.0	2	36 1
136	268	0.0	2	265 3
137	35	0.0	2	33 2
138	43	0.0	2	41 2
139	38	0.0	2	37 1
140	104	0.0	2	101 3
141	56	0.0	2	56
142	58	0.0	2	56 2
143	142	0.0	2	139 3
144	47	0.0	2	45 2
145	49	0.0	2	44 5
146	77	0.0	2	76 1
147	117	0.0	2	112 5
148	187	0.0	2	177 10
149	52	0.0	2	49 3
150	40	0.0	2	39 1
151	72	0.0	2	39 33


=== Second read: Adapter model_revcomp_Tn5 ===

Sequence: CTGTCTCTTATACACATCTCAACCGGCCGGCCGCTGGTCGACCTGCAGCGTACGNNNNNNNNNNNNNNNNNNNNAGAGACCTCGTGGACATCNNNNNN; Type: regular 3'; Length: 98; Trimmed: 91456 times.

No. of allowed errors:
0-49 bp: 0; 50-72 bp: 1

Bases preceding removed adapters:
  A: 98.7%
  C: 0.1%
  G: 0.1%
  T: 0.0%
  none/other: 1.0%
WARNING:
    The adapter is preceded by "A" extremely often.
    The provided adapter sequence could be incomplete at its 3' end.

Overview of removed sequences
length	count	expect	max.err	error counts
25	1263	0.0	0	1263
26	2404	0.0	0	2404
27	1765	0.0	0	1765
28	1371	0.0	0	1371
29	1322	0.0	0	1322
30	1022	0.0	0	1022
31	1082	0.0	0	1082
32	978	0.0	0	978
33	654	0.0	0	654
34	1363	0.0	0	1363
35	432	0.0	0	432
36	1406	0.0	0	1406
37	759	0.0	0	759
38	653	0.0	0	653
39	1653	0.0	0	1653
40	2265	0.0	0	2265
41	844	0.0	0	844
42	670	0.0	0	670
43	1467	0.0	0	1467
44	939	0.0	0	939
45	381	0.0	0	381
46	1224	0.0	0	1224
47	485	0.0	0	485
48	542	0.0	0	542
49	1382	0.0	0	1382
50	2247	0.0	1	2160 87
51	1119	0.0	1	1067 52
52	696	0.0	1	661 35
53	803	0.0	1	762 41
54	537	0.0	1	505 32
55	586	0.0	1	561 25
56	2091	0.0	1	2031 60
57	1098	0.0	1	1067 31
58	852	0.0	1	815 37
59	2336	0.0	1	2246 90
60	847	0.0	1	820 27
61	1055	0.0	1	1002 53
62	1473	0.0	1	1424 49
63	932	0.0	1	907 25
64	1183	0.0	1	1147 36
65	672	0.0	1	652 20
66	571	0.0	1	546 25
67	589	0.0	1	560 29
68	492	0.0	1	468 24
69	384	0.0	1	362 22
70	1541	0.0	1	1490 51
71	633	0.0	1	605 28
72	849	0.0	1	811 38
73	1020	0.0	1	959 61
74	3840	0.0	1	3717 123
75	2653	0.0	1	587 2066
76	351	0.0	1	311 40
77	805	0.0	1	774 31
78	1400	0.0	1	1363 37
79	469	0.0	1	452 17
80	674	0.0	1	649 25
81	547	0.0	1	522 25
82	154	0.0	1	143 11
83	197	0.0	1	191 6
84	748	0.0	1	723 25
85	327	0.0	1	321 6
86	449	0.0	1	438 11
87	713	0.0	1	694 19
88	352	0.0	1	340 12
89	395	0.0	1	383 12
90	259	0.0	1	255 4
91	342	0.0	1	333 9
92	881	0.0	1	850 31
93	271	0.0	1	269 2
94	290	0.0	1	278 12
95	476	0.0	1	456 20
96	667	0.0	1	640 27
97	674	0.0	1	653 21
98	143	0.0	1	136 7
99	815	0.0	1	791 24
100	196	0.0	1	193 3
101	184	0.0	1	178 6
102	760	0.0	1	740 20
103	243	0.0	1	230 13
104	173	0.0	1	166 7
105	1450	0.0	1	1399 51
106	391	0.0	1	378 13
107	486	0.0	1	470 16
108	489	0.0	1	477 12
109	256	0.0	1	248 8
110	301	0.0	1	294 7
111	1011	0.0	1	980 31
112	196	0.0	1	190 6
113	367	0.0	1	353 14
114	344	0.0	1	327 17
115	287	0.0	1	280 7
116	472	0.0	1	458 14
117	298	0.0	1	286 12
118	301	0.0	1	286 15
119	194	0.0	1	190 4
120	667	0.0	1	645 22
121	582	0.0	1	563 19
122	551	0.0	1	521 30
123	131	0.0	1	129 2
124	666	0.0	1	641 25
125	345	0.0	1	337 8
126	329	0.0	1	323 6
127	231	0.0	1	221 10
128	494	0.0	1	482 12
129	300	0.0	1	292 8
130	368	0.0	1	361 7
131	426	0.0	1	412 14
132	794	0.0	1	764 30
133	530	0.0	1	514 16
134	416	0.0	1	407 9
135	277	0.0	1	273 4
136	170	0.0	1	165 5
137	219	0.0	1	213 6
138	311	0.0	1	308 3
139	203	0.0	1	193 10
140	237	0.0	1	230 7
141	99	0.0	1	93 6
142	100	0.0	1	96 4
143	212	0.0	1	208 4
144	83	0.0	1	78 5
145	266	0.0	1	259 7
146	956	0.0	1	933 23
147	77	0.0	1	75 2
148	118	0.0	1	111 7
149	259	0.0	1	252 7
150	159	0.0	1	150 9
151	657	0.0	1	409 248


WARNING:
    One or more of your adapter sequences may be incomplete.
    Please see the detailed output above.


	- When the reverse complement of the model sequences were included in the model file this had no major effect on trimming the
	reverse reads in this trial. We do not see any of the model representing the sequence in the 5' direction on the forward read
	get trimmed off given the -A option. NO FORWARD READS GET TRIMMED. Thus there is no harm in using this file for trimming the 
	reverse reads independently.
	

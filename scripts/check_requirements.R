#!/usr/bin/env Rscript

print("Checking python (>=3.4)")
system('python3 --version')

print("Checking samtools >= 1.9")
system("samtools --version")

print("Checking cutadapt (>=2.6)")
if(as.numeric(system('cutadapt --version', intern=TRUE)) >= 2.6){
	print("cutadapt >= 2.6 found")
}

print("Checking bowtie2 (>= 2.3)")
system('bowtie2 --version')
	print("bowtie2 found")

print("Checking prodigal")
system('prodigal -v')

print("Checking pysam (>= 0.15)")
system('pip3 list | grep pysam')

print("Checking pandas ( >= 0.2)")
system('pip3 list | grep pandas')

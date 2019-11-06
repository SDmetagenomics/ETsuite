### Check and Load Libraries
if("seqinr" %in% installed.packages() == F){
  print("Please install R package seqinr. Program quitting...")
  q(save="no")
}

library(seqinr)


### Set up path variables for associated scripts and databases

# Get relative path of ETmapper install directory
initial.options <- commandArgs(trailingOnly = F)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)



### Collect and Parse arguments
args <- commandArgs(trailingOnly = T)

# Display help if no args or -h
if("-h" %in% args | !("-w" %in% args) | !("-d" %in% args) | !("-b" %in% args) | length(args) == 0) {
  cat("
  
    Usage: ETdb.R -g [genomes] -o [output_dir] [Additonal_Options]

    Mandatory Arguments:
      
      -g: Genomes in fasta format (No Default)

      
    Additional Options:

      -ad: Adapter sequence file (Default: db/adap_small.fa)
      -am: Min length of adapter match (Default: 5)
      -qs: Min base quality (Default: 20)"
      
      
      )

      q(save="no")
}



### Notes:
# Name for genome can be recovered from file name
#


list.files(path = "test_data/ETdb/genomes/",pattern = "*.fa" )




### Concatenate 

### Calculate Genome Statistics

fuck <- read.fasta("test_data/ETdb/genomes/AMD830_Bacillus_Bacilli_39_7.contigs.fa", seqtype = "DNA")

genome_dat <- vector(length = length(fuck), mode = "numeric")

for (i in 1:length(fuck)){
gc_vec[i] <- GC(fuck[[i]])

}


mean(gc_vec)



### Build Scaffold2Bin File

scaff2bin <- data.frame()

for (i in 1:length(files)){
  
  # tmp_scaffs <- names(fuck)
  # tmp_genome <- files[i]
  
  # tmp_s2b <- data.frame(bin = tmp_genome, scaffold = tmp_scaffs)
  
  # scaff2bin <- rbind(scaff2bin, tmp_s2b)
  
}


### Build bowtie2 index 






### Predict Proteins (if desired) - This will gnerate coordinates for later plotting ect



### 














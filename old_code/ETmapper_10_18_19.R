#! /usr/local/bin/Rscript

## TO DO
## 1. Create ReadTheDocs
## 2. 
## 3. 
## 4. 
## 5. 


### Check and Load Libraries
if("data.table" %in% installed.packages() == F){
  print("Please install R package data.table. Program quitting...")
  q(save="no")
}

library(data.table)




### Set up path variables for associated scripts and databases

# Get relative path of ETmapper install directory
initial.options <- commandArgs(trailingOnly = F)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

# Set up default database paths (note: normalizePath creates a full path from relative path)
ad <- normalizePath(paste0(script.basename,"/../db/adap_small.fa")) #adapter sequences
md <- normalizePath(paste0(script.basename,"/../db/models.fa")) #model sequences


# Testing directory stuff above
print(normalizePath(script.basename))
print(ad)




### Collect and Parse arguments
args <- commandArgs(trailingOnly = T)

# Display help if no args or -h
if("-h" %in% args | !("-w" %in% args) | !("-d" %in% args) | !("-b" %in% args) | length(args) == 0) {
  cat("
  
    ####### #######                                           
    #          #    #    #   ##   #####  #####  ###### #####  
    #          #    ##  ##  #  #  #    # #    # #      #    # 
    #####      #    # ## # #    # #    # #    # #####  #    # 
    #          #    #    # ###### #####  #####  #      #####  
    #          #    #    # #    # #      #      #      #   #  
    #######    #    #    # #    # #      #      ###### #    # 
      
    Usage: ETmapper.R -w [workflow] -d [read_dir] -b [batch_file] [Additonal_Options]

    Mandatory arguments:
      
      -w: Workflow type (No Default)
          jm - Junction mapping
          lm - Lite metagenomics coverage
      -d: Directory containing read files (No Default)
      -b: Batch file with sample information (No Default)

    Adapter Filtering Options:

      -ad: Adapter sequence file (Default: db/adap_small.fa)
      -al: Min length of adapter match (Default: 5)
      -qs: Min base quality (Default: 20)

    Junction mapping options:
    
      -md: Junction model sequence file (Default: db/models.fa)
      -et: Model match error (Default: 0.02)
      -mm: Model min match length (Default: 25)
      
    Program Control:
    
      -o: Output directory (Will be created if not specified)
      -cpu: Number of cores (Default: 1)
      -h: Bring up this help menu\n\n")
  
  
  q(save="no")
}


### Arg Testing

# Work Flow Type
wf <- args[which(args == "-w") + 1]
print(paste0("workflow type is: ", wf))

# Read Directory
rd <- args[which(args == "-d") + 1]
print(paste0("read directory is: ", rd))

# Batch File
bf <- args[which(args == "-b") + 1]
print(paste0("batch file is: ", bf))
batch_file <- read.table(bf, sep = "\t", header = F)

# Adapter Database File
ad <- ad
if("-ad" %in% args){
  ad <- args[which(args == "-ad") + 1]
  ad <- normalizePath(ad)
}
print(paste0("Adapter Database is: ", ad))

# MMin Length of Adapter Match (Default: 5)
al <- 5
if("-m" %in% args){
  mm_cut <- as.numeric(args[which(args == "-m") + 1])
}

# Min base quality score (Default: 20)
qs <- 20
if("-q" %in% args){
  mq_cut <- as.numeric(args[which(args == "-q") + 1])
}

# Model Database File
md <- md
if("-md" %in% args){
  md <- args[which(args == "-md") + 1]
  md <- normalizePath(md)
}
print(paste0("Adapter Database is: ", md))



# ### TEST ARGS
# wf <- "jm"
# rd <- "~/Dropbox/Banfield_Lab_Files/Projects/mCAFE/mCAFE_Project/ET_Mapper/test_data/ETmapper_test_data/reads/"
# bf <- "~/Dropbox/Banfield_Lab_Files/Projects/mCAFE/mCAFE_Project/ET_Mapper/test_data/ETmapper_test_data/test_batch.txt"
#af <- "~/Dropbox/Banfield_Lab_Files/Projects/mCAFE/mCAFE_Project/ET_Mapper/db/adap_small.fa"


if (wf == "jm"){

  ### Create Log File - WORKING!
  cat(
  paste0("ETmapper Log    Created: ", date()),"\n\n",
  "Program Parameters:\n\n",
  paste0("Workflow type is: ", wf),"\n",
  paste0("Total Samples: ",nrow(batch_file)),"\n",
  paste0("Adapter Trim DB: ", ad,"\n"),
  file = "ETmapper.log")
  
  
  
  ### Trimming and filtering paired reads
  
  # Run trimming for loop
  
  # for (i in 1:nrow(batch_file)){
  # 
  #   system(paste0("cutadapt -a file:",ad," -j 4 -O ",al," -q ",qs,
  #                 " -o ",batch_file[i,3],".trim",
  #                 " -p ",batch_file[i,4],".trim",
  #                 " ",rd,batch_file[i,3],
  #                 " ",rd,batch_file[i,4],
  #                 " > ",batch_file[i,1],".trim.log"))
  # 
  # }
  
  # Put trim metadata into spreadsheet
  # for loop for pulling cutadapt.log files 
  
  
  
  
  ### Identifying and Trimming Models
  
  # Run Model Finding cutadapt on fwd reads but include reverse 
  # info file should just give fwd reads...double check this 
  # cutadapt -g file:/Users/Spencer/Dropbox/Banfield_Lab_Files/Projects/mCAFE/mCAFE_Project/ET_Mapper/db/models.fa -e 0.02 -O 25 --discard-untrimmed --info-file tmp.txt -o JD_ZZ_2ndcycle_5_S5_L001_R1_001.fastq.trim -p JD_ZZ_2ndcycle_5_S5_L001_R2_001.fastq.trim JD_ZZ_2ndcycle_5_S5_L001_R1_001.fastq JD_ZZ_2ndcycle_5_S5_L001_R2_001.fastq
  

  # Identify flanking region to use as barcode query
  bc_flank <- system("grep -o \"NN.*.NN\" /Users/Spencer/Dropbox/Banfield_Lab_Files/Projects/mCAFE/mCAFE_Project/ET_Mapper/db/models.fa | sed 's/N//g'", intern = T)
  bc_flank <- unique(bc_flank)
  
  # pull filtered model hit file
  ca_info <- fread("test_data/ETmapper_test_data/reads/tmp_filt.txt", header = F)
  
  # create barcode output file
  bc_out <- data.frame(read = ca_info$V1, model = ca_info$V8, mod_len = ca_info$V4)
  
  # create tmp barcode storage frame and loop finding barcodes  
  bc_tmp_store <- data.frame()
  for (i in 1:length(bc_flank)){

    # parse ca_info for those with exact primer match
    bc_flank_i_matches <- ca_info[grep(bc_flank[i], ca_info$V6),c(1,6)]

    # Find upstream 20bp of target sequence and sub out target for nothing
    matches <- regexpr(paste0(bc_flank[i],".{0,20}"), bc_flank_i_matches$V6, perl = T)
    bc_flank_i_matches$barcodes <- sub(bc_flank[i], "", regmatches(bc_flank_i_matches$V6, matches))
    
    # Create final data frame w/ column for flank and match and cbind to outupt
    bc_flank_i_matches <- data.frame(bc_flank_i_matches[,-2], flank_seq = bc_flank[i])
    bc_tmp_store <- rbind(bc_tmp_store, bc_flank_i_matches)
  
  }
    
  # merge barcodes to output file
  bc_out <- merge(bc_out, bc_tmp_store, by.x = "read", by.y = "V1", all.x = T)
  
  
  # step 3: output read barcodes to file as well as log number of barcodes per samples

} else {
  print("MetaG workflow not ready yet!")
}






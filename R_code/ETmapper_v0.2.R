#!/usr/bin/env Rscript

### Load Test Data
#foo <- read.table("test_data/ETmapper_test_data/test_batch.txt")



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
ad <- normalizePath(paste0(script.basename,"/../db/ETseq_2_4_adap.fa")) #adapter sequences
md <- normalizePath(paste0(script.basename,"/../db/models.fa")) #model sequences
scripts <- normalizePath(paste0(script.basename,"/../scripts")) #accesory scripts directory


# Testing directory stuff above
#print(normalizePath(script.basename))
#print(ad)

# NOT IN Operator for Arg Parseing
'%notin%' <- Negate('%in%')

### Collect and Parse arguments
args <- commandArgs(trailingOnly = T)

# Display help if no args or -h
if("-h" %in% args | !("-w" %in% args) | !("-d" %in% args) | !("-b" %in% args) | !("-g" %in% args) | length(args) == 0) {
  cat("
  
    ####### #######                                           
    #          #    #    #   ##   #####  #####  ###### #####  
    #          #    ##  ##  #  #  #    # #    # #      #    # 
    #####      #    # ## # #    # #    # #    # #####  #    # 
    #          #    #    # ###### #####  #####  #      #####  
    #          #    #    # #    # #      #      #      #   #  
    #######    #    #    # #    # #      #      ###### #    # 
      
    Usage: ETmapper.R -w [workflow] -d [read_dir] -b [batch_file] -g [genome_db] [Additonal_Options]

    Mandatory Arguments:
    
      -w: Workflow type (No Default)
          jm - Junction mapping
          lm - Lite metagenomics coverage
      -d: Directory containing read files (No Default)
      -b: Batch file with sample information (No Default)
      -g: Directory containing genome database (No Default)

    Adapter Filtering Options:
    
      -ad: Adapter sequence file (Default: db/adap_small.fa)
      -am: Min length of adapter match (Default: 5)
      -qs: Min base quality (Default: 20)

    Model Identification Options:
    
      -md: Junction model sequence file (Default: db/models.fa)
      -mm: Min length of model match (Default: 25)
      -et: Model match error (Default: 0.02)
      -rl: Min final read length (Default: 40)
      
    Read Mapping Options:
    
      -X: Maximum insert length (Default: 500)
      -F: Map forward reads only (Only for pe data)
      
    Program Control:
    
      -o: Output directory (Will be created if not specified)
      -cpu: Number of cores (Default: 1)
      -h: Bring up this help menu\n\n")
  
  
  q(save="no")
}


## Mandatory Arguments

# Work Flow Type
wf <- args[which(args == "-w") + 1]
if(wf %notin% c("jm","lm")){
  cat(paste0("\n",wf, " is not a kown workflow...exiting"))
  #q(save="no")
}

# Read Directory
rd <- args[which(args == "-d") + 1]

# Batch File
bf <- args[which(args == "-b") + 1]
batch_file <- read.table(bf, sep = "\t", header = F)

# Genome Database Directory
gd <- args[which(args == "-g") + 1]


## Adapter Filtering Arguments

# Adapter Database File
ad <- ad
if("-ad" %in% args){
  ad <- args[which(args == "-ad") + 1]
  ad <- normalizePath(ad)
}

# MMin Length of Adapter Match (Default: 5)
am <- 5
if("-m" %in% args){
  am <- as.numeric(args[which(args == "-am") + 1])
}

# Min base quality score (Default: 20)
qs <- 20
if("-q" %in% args){
  qs <- as.numeric(args[which(args == "-qs") + 1])
}


## Model Identification Options

# Model Database File
md <- md
if("-md" %in% args){
  md <- args[which(args == "-md") + 1]
  md <- normalizePath(md)
}

# Min length of model match  (Default: 25)
mm <- 25
if("-mm" %in% args){
  mm <- as.numeric(args[which(args == "-mm") + 1])
}

# Max model match error (Default: 0.02)
et <- 0.02
if("-et" %in% args){
  et <- as.numeric(args[which(args == "-et") + 1])
}

# Min final read length (Default: 40)
rl <- 40
if("-rl" %in% args){
  rl <- as.numeric(args[which(args == "-rl") + 1])
}


## Read Mapping Options

# Max insert length (Default: 500 bp)
isl <- 500
if("-X" %in% args){
  isl <- as.numeric(args[which(args == "-X") + 1])
}

# Map forward reads only (Default: FALSE)
fwd_only <- FALSE
if("-F" %in% args){
  fwd_only <- TRUE
}


## Program Control Options

# number of cores to use where applicable
cpu <- 1
if("-cpu" %in% args){
  cpu <- as.numeric(args[which(args == "-cpu") + 1])
}



###### HOW TO FIX NON-ARGUMENT ENTRY BUG i.e. ETmapper.R -w [workflow] -d -b [batch_file]
# for important arguments check that: 
# 
# which(argument) + 1 != c(all_other_arguments)
# 
# If somebody forgets to enter something then the thing +1 from the argument in the vector will
# be another argument 



#### BEGIN FUNCTION DEFINE ####



#### END FUNCTION DEFINE ####





#### BEGIN WORKFLOWS ####


###############
############### JUNCTION MAPPING WORKFLOW
###############


if (wf == "jm"){
  
  
  
  ### Create Log File - WORKING!
  cat(
  paste0("ETmapper Log    Created: ", date()),"\n\n",
  "Program Parameters:\n\n",
  paste0("Workflow type is: ", wf),"\n",
  paste0("Total Samples: ",nrow(batch_file)),"\n",
  paste0("Adapter Trim DB: ", ad,"\n"),
  file = "ETmapper.log")
  
  
  ### Determine if reads provided are se or pr end and initalize program
  if(ncol(batch_file) == 3){
    paired_end_data <- FALSE
    cat("Input Reads Identifed as Single End...Begining Analysis\n\n")
    cat("Paired End = FALSE", file = "ETmapper.log", append = T)
  }
  
  if(ncol(batch_file) == 4){
    paired_end_data <- TRUE
    cat("Input Data Identifed as Paired End...Begining Analysis\n\n")
    cat("Paired End = TRUE", file = "ETmapper.log", append = T)
  }
  
    
  
  ### Begin Trimming / Mapping Steps of Junction Mapping Workflow (PAIRED END)
  if(paired_end_data == TRUE){

    ### Run adapter/flanking sequence trimming AND Quality Score Filtering
    for (i in 1:nrow(batch_file)){
    
      # Indicate what program is doing
      cat(paste0("\nTrimming Adapters/Flanks and Qscore Filter for: ",batch_file[i,1],"\n"))
    
      # run cutadapt command
      system(paste0("cutadapt -a file:",ad," -A file:",ad, # specify adapter types and file
                    " -j ",cpu," -O ",am," -q ",qs, # specify trimming params
                    " -o ",batch_file[i,3],".trim", # fwd output
                    " -p ",batch_file[i,4],".trim", # rev output
                    " ",rd,batch_file[i,3], # fwd read
                    " ",rd,batch_file[i,4], # rev read
                    " > ",batch_file[i,1],".trim.log")) # log files

    }
  
  
    ### Identifying and Trimming Models (ONLY 1 THREAD POSSIBLE w/ --info-file)
  
    # Run model finding for loop on fwd reads only but include reverse 
    for (i in 1:nrow(batch_file)){
    
      # Indicate what program is doing
      cat(paste0("\nFinding and Trimming Models for: ",batch_file[i,1],"\n"))
    
      # Run cutadapt command 
      system(paste0("cutadapt -g file:",md, # specify model file
                    " -O ",mm," -e ",et, # specify trimming params
                    " --discard-untrimmed", # discard read pairs without model
                    " --info-file ",batch_file[i,1],".info", # print 
                    " -o ",batch_file[i,3],".clean", # fwd output
                    " -p ",batch_file[i,4],".clean", # rev output
                    " ",batch_file[i,3],".trim", # fwd read
                    " ",batch_file[i,4],".trim", # rev read
                    " > ",batch_file[i,1],".clean.log")) # log files
    
    }
  
  
    ### Identifying barcodes and creating barcode db file
  
    # Identify flanking region to use as barcode query
    bc_flank <- system(paste0("grep -o \"NN.*.NN\" ",md," | sed 's/N//g'"), intern = T)
    bc_flank <- unique(bc_flank)
  
    for (i in 1:nrow(batch_file)){
  
      # Indicate what program is doing
      cat(paste0("\nFinding Barcodes for: ",batch_file[i,1],"\n"))
    
      # Filter cutadapt info file for reads with model
      system(paste0("awk '{if ($3!=-1) print}' ",batch_file[i,1],".info > ",batch_file[i,1],".info.filt"))
  
      # pull filtered model hit file
      ca_info <- fread(paste0(batch_file[i,1],".info.filt"), header = F)

      # create barcode output file
      bc_out <- data.frame(Read = ca_info$V1, model = ca_info$V8, mod_len = ca_info$V4)
  
      # create tmp barcode storage data frame
      bc_tmp_store <- data.frame()
    
      # loop over possible barcode flanking sequences 
      for (j in 1:length(bc_flank)){
      
        # parse ca_info for those with exact primer match
        bc_flank_i_matches <- ca_info[grep(bc_flank[j], ca_info$V6),c(1,6)]
      
        # skip iteration if no flank sequence match
        if(nrow(bc_flank_i_matches) == 0){
          next
        }
      
        # Find upstream 20bp of target sequence and sub out target for nothing
        matches <- regexpr(paste0(bc_flank[j],".{0,20}"), bc_flank_i_matches$V6, perl = T)
        bc_flank_i_matches$barcodes <- sub(bc_flank[j], "", regmatches(bc_flank_i_matches$V6, matches))
      
        # Create final data frame w/ column for flank and match and cbind to outupt
        bc_flank_i_matches <- data.frame(bc_flank_i_matches[,-2], flank_seq = bc_flank[j])
        bc_tmp_store <- rbind(bc_tmp_store, bc_flank_i_matches)
      
      }
    
      # merge barcodes to output file and remove sequencing barcode line
      bc_out <- merge(bc_out, bc_tmp_store, by.x = "Read", by.y = "V1", all.x = T)
      bc_out$Read <- sub(pattern = " .*$", replacement = "", bc_out$Read, perl = T)
    
      # write barcodes out 
      write.table(bc_out, paste0(batch_file[i,1],".bc"), row.names = F)
  
    }
  
    
    ### Trimming Rev reads and filtering for length

    # Run model finding for loop on rev reads and include fwd reads
    for (i in 1:nrow(batch_file)){

      # Indicate what program is doing
      cat(paste0("\nTrimming Rev Read/Size Filter: ",batch_file[i,1],"\n"))
    
      # run cutadapt command
      system(paste0("cutadapt -A file:",md, # specify model file
                    " -O 10"," -j ",cpu," -e ",et, # specify trimming params
                    " --minimum-length ",rl, # discard read pairs if both reads are < rl bp
                    " --pair-filter=both", # At least 1 read must be >= rl
                    " -o ",batch_file[i,3],".clean2", # fwd output
                    " -p ",batch_file[i,4],".clean2", # rev output
                    " ",batch_file[i,3],".clean", # fwd read
                    " ",batch_file[i,4],".clean", # rev read
                    " > ",batch_file[i,1],".clean2.log")) # log files
    
    }
  
  
    ### Run bowtie mapping 
  
    # Map fwd and rev read pairs
    if(fwd_only == FALSE){

      # Run bowtie mapping for loop for PAIRED END mapping
      for (i in 1:nrow(batch_file)){

      # Indicate what program is doing
      cat(paste0("\nMapping: ",batch_file[i,1],"\n"))

      # run bowtie mapping
      system(paste0("bowtie2 -x ",gd,"/bt2/All_Genomes", # specify genome database
                    " -p ",cpu," -X ",isl, # specify bowtie options
                    " -1 ",batch_file[i,3],".clean2", # specify fwd reads
                    " -2 ",batch_file[i,4],".clean2", # specify rev reads
                    " -S ",batch_file[i,1],".sam", # specify sam file output
                    " 2> ",batch_file[i,1],".bowtie.log")) # specify log file output
    
      # Convert SAM to BAM format with sorting and indexing
      system(paste0("samtools view -S -b ",batch_file[i,1],".sam > ",batch_file[i,1],".bam; samtools sort ",batch_file[i,1],".bam -o ",batch_file[i,1],".bam.sorted; samtools index ",batch_file[i,1],".bam.sorted"))
                  
      # Run read hit stats script
      system(paste0("python3 ",scripts,"/bam_pe_stats.py ",gd,"/scaff2bin.txt ", batch_file[i,1],".bam.sorted > ",batch_file[i,1],".hits"))
    
      # Integrate hit reads with barcodes into combined output and write
      hit_dat <- fread(paste0(batch_file[i,1],".hits"), header = T, stringsAsFactors = F)
      bc_dat <- fread(paste0(batch_file[i,1],".bc"), header = T, stringsAsFactors = F)
      merge_dat <- merge(hit_dat, bc_dat, by = "Read", all.x = T)
    
      write.table(merge_dat, paste0(batch_file[i,1],".hits2"), row.names = F, quote = F, sep = "\t")
      cat(paste0("Finished reading a BAM, wrote output to ", paste0(batch_file[i,1],".hits2")))
      cat("\nSneak peak of unfiltered results:\n")
      print(table(merge_dat[which(merge_dat$GENOME1 == merge_dat$GENOME2),]$GENOME1))

    
      }
    
    }
  
  
    if(fwd_only == TRUE){
    
      # Run bowtie mapping for loop for SINGLE END mapping
      for (i in 1:nrow(batch_file)){
      
        # Indicate what program is doing
        cat(paste0("\nMapping: ",batch_file[i,1],"\n"))
      
        # run bowtie mapping
        system(paste0("bowtie2 -x ",gd,"/bt2/All_Genomes", # specify genome database
                      " -p ",cpu, # specify bowtie options
                      " -U ",batch_file[i,3],".clean2", # specify fwd reads
                      " -S ",batch_file[i,1],".sam", # specify sam file output
                      " 2> ",batch_file[i,1],".bowtie.log")) # specify log file output
      
        # Convert SAM to BAM format with sorting and indexing
        system(paste0("samtools view -S -b ",batch_file[i,1],".sam > ",batch_file[i,1],".bam; samtools sort ",batch_file[i,1],".bam -o ",batch_file[i,1],".bam.sorted; samtools index ",batch_file[i,1],".bam.sorted"))
      
        # Run read hit stats script
        system(paste0("python3 ",scripts,"/bam_se_stats.py ",gd,"/scaff2bin.txt ", batch_file[i,1],".bam.sorted > ",batch_file[i,1],".hits"))
      
        # Integrate hit reads with barcodes into combined output and write
        hit_dat <- fread(paste0(batch_file[i,1],".hits"), header = T, stringsAsFactors = F)
        bc_dat <- fread(paste0(batch_file[i,1],".bc"), header = T, stringsAsFactors = F)
        merge_dat <- merge(hit_dat, bc_dat, by= "Read", all.x = T)
      
        write.table(merge_dat, paste0(batch_file[i,1],".hits2"), row.names = F, quote = F, sep = "\t")
       
        cat(paste0("Finished reading a BAM, wrote output to ", paste0(batch_file[i,1],".hits2")))
        cat("\nSneak peak of unfiltered results:\n")
        print(table(merge_dat$GENOME1))

      }
  
    }

  } ### END Trimming / Mapping Steps of Junction Mapping Workflow (PAIRED END)
  
  
  
  
  ### Begin Trimming / Mapping Steps of Junction Mapping Workflow (SINGLE END)
  if(paired_end_data == FALSE){
    
    
    ### Run adapter/flanking sequence trimming AND Quality Score Filtering
    for (i in 1:nrow(batch_file)){
      
      # Indicate what program is doing
      cat(paste0("\nTrimming Adapters/Flanks and Qscore Filter for: ",batch_file[i,1],"\n"))
      
      # run cutadapt command
      system(paste0("cutadapt -a file:",ad, # specify adapter types and file
                    " -j ",cpu," -O ",am," -q ",qs, # specify trimming params
                    " -o ",batch_file[i,3],".trim", # fwd output
                    " ",rd,batch_file[i,3], # fwd read
                    " > ",batch_file[i,1],".trim.log")) # log files
      
    }
  
    
    ### Run model finding for loop on fwd reads only AND Size Filter HERE (ONLY 1 THREAD POSSIBLE w/ --info-file)
    for (i in 1:nrow(batch_file)){
      
      # Indicate what program is doing
      cat(paste0("\nFinding and Trimming Models and Length Filter for: ",batch_file[i,1],"\n"))
      
      # Run cutadapt command 
      system(paste0("cutadapt -g file:",md, # specify model file
                    " -O ",mm," -e ",et, # specify trimming params
                    " --discard-untrimmed", # discard read pairs without model
                    " --minimum-length ",rl, # discard read if length < rl bp
                    " --info-file ",batch_file[i,1],".info", # print 
                    " -o ",batch_file[i,3],".clean", # fwd output
                    " ",batch_file[i,3],".trim", # fwd read
                    " > ",batch_file[i,1],".clean.log")) # log files
      
    }
  
  
    ### Identifying barcodes and creating barcode db file
    
    # Identify flanking region to use as barcode query
    bc_flank <- system(paste0("grep -o \"NN.*.NN\" ",md," | sed 's/N//g'"), intern = T)
    bc_flank <- unique(bc_flank)
    
    for (i in 1:nrow(batch_file)){
      
      # Indicate what program is doing
      cat(paste0("\nFinding Barcodes for: ",batch_file[i,1],"\n"))
      
      # Filter cutadapt info file for reads with model
      system(paste0("awk '{if ($3!=-1) print}' ",batch_file[i,1],".info > ",batch_file[i,1],".info.filt"))
      
      # pull filtered model hit file
      ca_info <- fread(paste0(batch_file[i,1],".info.filt"), header = F)
      
      # create barcode output file
      bc_out <- data.frame(Read = ca_info$V1, model = ca_info$V8, mod_len = ca_info$V4)
      
      # create tmp barcode storage data frame
      bc_tmp_store <- data.frame()
      
      # loop over possible barcode flanking sequences 
      for (j in 1:length(bc_flank)){
        
        # parse ca_info for those with exact primer match
        bc_flank_i_matches <- ca_info[grep(bc_flank[j], ca_info$V6),c(1,6)]
        
        # skip iteration if no flank sequence match
        if(nrow(bc_flank_i_matches) == 0){
          next
        }
        
        # Find upstream 20bp of target sequence and sub out target for nothing
        matches <- regexpr(paste0(bc_flank[j],".{0,20}"), bc_flank_i_matches$V6, perl = T)
        bc_flank_i_matches$barcodes <- sub(bc_flank[j], "", regmatches(bc_flank_i_matches$V6, matches))
        
        # Create final data frame w/ column for flank and match and cbind to outupt
        bc_flank_i_matches <- data.frame(bc_flank_i_matches[,-2], flank_seq = bc_flank[j])
        bc_tmp_store <- rbind(bc_tmp_store, bc_flank_i_matches)
        
      }
      
      # merge barcodes to output file and remove sequencing barcode line
      bc_out <- merge(bc_out, bc_tmp_store, by.x = "Read", by.y = "V1", all.x = T)
      bc_out$Read <- sub(pattern = " .*$", replacement = "", bc_out$Read, perl = T)
      
      # write barcodes out 
      write.table(bc_out, paste0(batch_file[i,1],".bc"), row.names = F)
      
    }  
  
    # Run bowtie mapping for loop for SINGLE END mapping
    for (i in 1:nrow(batch_file)){
      
      # Indicate what program is doing
      cat(paste0("\nMapping: ",batch_file[i,1],"\n"))
      
      # run bowtie mapping
      system(paste0("bowtie2 -x ",gd,"/bt2/All_Genomes", # specify genome database
                    " -p ",cpu, # specify bowtie options
                    " -U ",batch_file[i,3],".clean", # specify fwd reads
                    " -S ",batch_file[i,1],".sam", # specify sam file output
                    " 2> ",batch_file[i,1],".bowtie.log")) # specify log file output
      
      # Convert SAM to BAM format with sorting and indexing
      system(paste0("samtools view -S -b ",batch_file[i,1],".sam > ",batch_file[i,1],".bam; samtools sort ",batch_file[i,1],".bam -o ",batch_file[i,1],".bam.sorted; samtools index ",batch_file[i,1],".bam.sorted"))
      
      # Run read hit stats script
      system(paste0("python3 ",scripts,"/bam_se_stats.py ",gd,"/scaff2bin.txt ", batch_file[i,1],".bam.sorted > ",batch_file[i,1],".hits"))  #****REMOVE HARDCODE
      
      # Integrate hit reads with barcodes into combined output and write
      hit_dat <- fread(paste0(batch_file[i,1],".hits"), header = T, stringsAsFactors = F)
      bc_dat <- fread(paste0(batch_file[i,1],".bc"), header = T, stringsAsFactors = F)
      merge_dat <- merge(hit_dat, bc_dat, by = "Read", all.x = T)
      
      write.table(merge_dat, paste0(batch_file[i,1],".hits2"), row.names = F, quote = F, sep = "\t")
      
      cat(paste0("Finished reading a BAM, wrote output to ", paste0(batch_file[i,1],".hits2")))
      cat("\nSneak peak of unfiltered results:\n")
      cat(print(table(merge_dat$GENOME1)))
    } # End bowtie mapping for loop for SINGLE END mapping

  } ### End Trimming / Mapping Steps of Junction Mapping Workflow (SINGLE END) 

  cat("Junction mapping finished successfully.")
} ### End Junction Mapping Workflow




###############
############### LITE META-G WORKFLOW
###############


if (wf == "lm") {
  print("MetaG workflow not ready yet!")
  q(save="no")
}






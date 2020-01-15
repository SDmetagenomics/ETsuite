#!/usr/bin/env Rscript

### Load Test Data
#batch_file <- read.table("test_data/ETmapper/test_batch.txt")



### Check and Load Libraries
if ("data.table" %in% installed.packages() == F){
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



### Collect and Parse arguments
args <- commandArgs(trailingOnly = T)

# NOT IN Operator for Arg Parseing
'%notin%' <- Negate('%in%')

# Display help if no args or -h
if("-h" %in% args | !("-w" %in% args) | !("-d" %in% args) | !("-b" %in% args) | !("-g" %in% args) | length(args) == 0) {
  cat("
  
    ####### #######                                           
    #          #    #    #   ##   #####  #####  ###### #####  
    #          #    ##  ##  #  #  #    # #    # #      #    # 
    #####      #    # ## # #    # #    # #    # #####  #    # 
    #          #    #    # ###### #####  #####  #      #####  
    #          #    #    # #    # #      #      #      #   #  
    #######    #    #    # #    # #      #      ###### #    #  v0.03
      
    Usage: ETmapper.R -w [workflow] -d [read_dir] -b [batch_file] -g [genome_db] [options]

    Mandatory Arguments:
    
      -w: Workflow type (No Default)
          jm - Junction mapping
          lm - Lite metagenomics coverage
      -d: Directory containing read files (No Default)
      -b: Batch file with sample information (No Default)
      -g: Directory containing genome database (No Default)

    Adapter Filtering Options:
    
      -ad: Adapter sequence file (Default: db/ETseq_2_4_adap.fa)
      -am: Min length of adapter match (Default: 5)
      -qs: Min base quality (Default: 20)

    Model Identification Options:
    
      -md: Junction model sequence file (Default: db/models.fa)
      -mm: Min length of model match (Default: 25)
      -et: Model match error (Default: 0.02)
      -rl: Min final read length (Default: 40)
      
    Read Mapping Options:
    
      -X: Maximum insert length (Default: 500)
      
    Program Control:
    
      -o: Output directory (Will be created if not specified)
      -cpu: Number of cores (Default: 1)
      -h: Bring up this help menu\n\n")
  
  
  q(save="no")
}

# NOTE: Banner is generated with http://www.bagill.com/ascii-sig.php ; Banner Font


## Mandatory Arguments

# Work Flow Type
wf <- args[which(args == "-w") + 1]
if(wf %notin% c("jm","lm")){
  cat(paste0("\n",wf, " is not a kown workflow...exiting"))
  q(save="no")
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



## Program Control Options

# number of cores to use where applicable
cpu <- 1
if("-cpu" %in% args){
  cpu <- as.numeric(args[which(args == "-cpu") + 1])
}

# Output directory (Will be created if not specified)
out_dir <- "ET_mapper"
if("-o" %in% args){
  out_dir <- args[which(args == "-o") + 1]
}
dir.create(out_dir, recursive = T)


###### HOW TO FIX NON-ARGUMENT ENTRY BUG i.e. ETmapper.R -w [workflow] -d -b [batch_file]
# for important arguments check that: 
# 
# which(argument) + 1 != c(all_other_arguments)
# 
# If somebody forgets to enter something then the thing +1 from the argument in the vector will
# be another argument 



#### BEGIN FUNCTION DEFINE ####


## Function 1: Clean up files and make output structure
clean.up <- function(){

  # Clean jm Workflow
  if (wf == "jm"){
  
    # say waht is happening
    cat("\n\nCleaning Up...\n")
  
    # remove un-needed files
    system(paste0("rm ",out_dir,"/*.tmpbam ",out_dir,"/*.info ",out_dir,"/*.sam ",out_dir,"/*.tmphits"))

    # create subdirectories
    dir.create(paste0(out_dir,"/logs"), recursive = T)
    dir.create(paste0(out_dir,"/hits"), recursive = T)
    dir.create(paste0(out_dir,"/map"), recursive = T)
    dir.create(paste0(out_dir,"/raw"), recursive = T)

    # move files to right places
    system(paste0("mv ",out_dir,"/*.log ",out_dir,"/logs/"))
    system(paste0("mv ",out_dir,"/*.hits ",out_dir,"/hits/"))
    system(paste0("mv ",out_dir,"/*.sorted.bam ",out_dir,"/*.sorted.bam.bai ",out_dir,"/map/"))
    system(paste0("mv ",out_dir,"/*.trim ",out_dir,"/*.clean ",out_dir,"/*.clean2 ",out_dir,"/*.bc ",out_dir,"/*info.filt ",out_dir,"/raw/"), ignore.stderr = T) 

  }
  
  # Clean lm Workflow
  if (wf == "lm"){
    
    # say waht is happening
    cat("\n\nCleaning Up...\n")
    
    # remove un-needed files
    system(paste0("rm ",out_dir,"/*.tmpbam ",out_dir,"/*.sam "))
    
    # create subdirectories
    dir.create(paste0(out_dir,"/logs"), recursive = T)
    dir.create(paste0(out_dir,"/hits"), recursive = T)
    dir.create(paste0(out_dir,"/map"), recursive = T)
    dir.create(paste0(out_dir,"/raw"), recursive = T)
    
    # move files to right places
    system(paste0("mv ",out_dir,"/*.log ",out_dir,"/logs/"))
    system(paste0("mv ",out_dir,"/*.mghits ",out_dir,"/hits/"))
    system(paste0("mv ",out_dir,"/*.sorted.bam ",out_dir,"/*.sorted.bam.bai ",out_dir,"/map/"))
    system(paste0("mv ",out_dir,"/*.trim ",out_dir,"/raw/"), ignore.stderr = T) 
    
  }
  
}


## Function 2: Pull run statistics from log files
pull.run.stats <- function(){
  
  # Pull stats jm workflow
  if (wf == "jm"){
    
    
    # make function to count Tn-primers
    count.primer <- function(){
      
      # make vector to hold total Tn-primer matches
      primer_counts <- c(0)
      
      # loop through possible Tn-primers in raw Fwd reads, return counts for each, and add to primer_coutns vector
      for (j in 1:length(bc_flank)){
        tmp_counts <- as.numeric(system(paste0("grep -c '",bc_flank[j],"' ",rd,batch_file[i,3]), intern = T))
        primer_counts <- primer_counts + tmp_counts
      }
      
      # return primer counts
      primer_counts
    }
    
    
    # Create dataframe to hold output
    jm_workflow_stats <- data.frame(batch_file,
                                    Total_Reads = 0, #5,4
                                    R1_adap = 0, #6,5
                                    R1_adap_frac = 0,
                                    R2_adap = 0, #8,7
                                    R2_adap_frac = 0,
                                    Tn_Primer = 0, #10,9
                                    Tn_Primer_Frac = 0,
                                    Model_Keep = 0, #12,11
                                    Model_Keep_Frac = 0,
                                    R2_Model = 0, #14,13
                                    R2_Model_Frac = 0,
                                    Good_Keep = 0, #16,15
                                    Good_Keep_Frac = 0,
                                    Raw_Map = 0, #18,17
                                    Raw_Map_Frac = 0)
    
    # Pull stats from jm trimming logs
    for (i in 1:nrow(batch_file)){
      
      # Say what is happening
      cat(paste0("\nAggregating ETmapper stats for: ",batch_file[i,1],"\n"))
      
      # Pull stats from logs
      if(paired_end_data == TRUE){
        jm_workflow_stats[i,5] <- as.numeric(system(paste0("grep 'Total read pairs processed:' ",out_dir,"/logs/",batch_file[i,1],".trim.log ","| awk '{print $5}' | sed 's/,//g'"), intern = T))
        jm_workflow_stats[i,6] <- as.numeric(system(paste0("grep 'Read 1 with adapter:' ",out_dir,"/logs/",batch_file[i,1],".trim.log ","| awk '{print $5}' | sed 's/,//g'"), intern = T))
        jm_workflow_stats[i,8] <- as.numeric(system(paste0("grep 'Read 2 with adapter:' ",out_dir,"/logs/",batch_file[i,1],".trim.log ","| awk '{print $5}' | sed 's/,//g'"), intern = T))
        jm_workflow_stats[i,10] <- count.primer()
        jm_workflow_stats[i,12] <- as.numeric(system(paste0("grep 'Pairs written' ",out_dir,"/logs/",batch_file[i,1],".clean.log ","| awk '{print $5}' | sed 's/,//g'"), intern = T))
        jm_workflow_stats[i,14] <- as.numeric(system(paste0("grep 'Read 2 with adapter:' ",out_dir,"/logs/",batch_file[i,1],".clean2.log ","| awk '{print $5}' | sed 's/,//g'"), intern = T))
        jm_workflow_stats[i,16] <- as.numeric(system(paste0("grep 'Pairs written' ",out_dir,"/logs/",batch_file[i,1],".clean2.log ","| awk '{print $5}' | sed 's/,//g'"), intern = T))
        jm_workflow_stats[i,18] <- as.numeric(system(paste0("sed 1d ",out_dir,"/hits/",batch_file[i,1],".hits | wc -l"), intern = T))
      }
      
      ### NOT FUNCTIONAL YET
      # if(paired_end_data == FALSE){
      #   lm_workflow_stats[i,4] <- system(paste0("grep 'Total reads processed:' ",out_dir,"/logs/",batch_file[i,1],".trim.log ","| awk '{print $4}' | sed 's/,//g'"), intern = T)
      #   lm_workflow_stats[i,5] <- system(paste0("grep 'Reads with adapters:' ",out_dir,"/logs/",batch_file[i,1],".trim.log ","| awk '{print $4}' | sed 's/,//g'"), intern = T)
      #   lm_workflow_stats[i,7] <- NA
      #   lm_workflow_stats[i,9] <- system(paste0("grep 'Reads written' ",out_dir,"/logs/",batch_file[i,1],".trim.log ","| awk '{print $5}' | sed 's/,//g'"), intern = T)
      #   lm_workflow_stats[i,11] <- system(paste0("sed 1d ",out_dir,"/hits/",batch_file[i,1],".mghits | wc -l"), intern = T)
      # }
      
    }
    
    # Calculate fraction columns
    jm_workflow_stats$R1_adap_frac <- jm_workflow_stats$R1_adap / jm_workflow_stats$Total_Reads
    jm_workflow_stats$R2_adap_frac <- jm_workflow_stats$R2_adap / jm_workflow_stats$Total_Reads
    jm_workflow_stats$Tn_Primer_Frac <- jm_workflow_stats$Tn_Primer / jm_workflow_stats$Total_Reads
    jm_workflow_stats$Model_Keep_Frac <- jm_workflow_stats$Model_Keep / jm_workflow_stats$Total_Reads
    jm_workflow_stats$R2_Model_Frac <- jm_workflow_stats$R2_Model / jm_workflow_stats$Model_Keep
    jm_workflow_stats$Good_Keep_Frac <- jm_workflow_stats$Good_Keep / jm_workflow_stats$Total_Reads
    jm_workflow_stats$Raw_Map_Frac <- jm_workflow_stats$Raw_Map / jm_workflow_stats$Total_Reads
    
    # Output data.frame
    jm_workflow_stats
    
  } # END JM BRANCH
  
  
  # Pull stats lm workflow
  if (wf == "lm"){
    
    # Create dataframe to hold output
    lm_workflow_stats <- data.frame(batch_file,
                                    Total_Reads = 0,
                                    R1_adap = 0,
                                    R1_adap_frac = 0,
                                    R2_adap = 0,
                                    R2_adap_frac = 0,
                                    Good_Keep = 0,
                                    Good_Keep_Frac = 0,
                                    Raw_Map = 0,
                                    Raw_Map_Frac = 0)
    
    # Pull stats from trimming logs
    for (i in 1:nrow(batch_file)){
      
      # Say what is happening
      cat(paste0("\nAggregating ETmapper stats for: ",batch_file[i,1],"\n"))
      
      # Pull stats from logs
      if(paired_end_data == TRUE){
        lm_workflow_stats[i,5] <- as.numeric(system(paste0("grep 'Total read pairs processed:' ",out_dir,"/logs/",batch_file[i,1],".trim.log ","| awk '{print $5}' | sed 's/,//g'"), intern = T))
        lm_workflow_stats[i,6] <- as.numeric(system(paste0("grep 'Read 1 with adapter:' ",out_dir,"/logs/",batch_file[i,1],".trim.log ","| awk '{print $5}' | sed 's/,//g'"), intern = T))
        lm_workflow_stats[i,8] <- as.numeric(system(paste0("grep 'Read 2 with adapter:' ",out_dir,"/logs/",batch_file[i,1],".trim.log ","| awk '{print $5}' | sed 's/,//g'"), intern = T))
        lm_workflow_stats[i,10] <- as.numeric(system(paste0("grep 'Pairs written' ",out_dir,"/logs/",batch_file[i,1],".trim.log ","| awk '{print $5}' | sed 's/,//g'"), intern = T))
        lm_workflow_stats[i,12] <- as.numeric(system(paste0("sed 1d ",out_dir,"/hits/",batch_file[i,1],".mghits | wc -l"), intern = T))
      }
      
      ### STILL UNTESTED
      if(paired_end_data == FALSE){
        lm_workflow_stats[i,4] <- system(paste0("grep 'Total reads processed:' ",out_dir,"/logs/",batch_file[i,1],".trim.log ","| awk '{print $4}' | sed 's/,//g'"), intern = T)
        lm_workflow_stats[i,5] <- system(paste0("grep 'Reads with adapters:' ",out_dir,"/logs/",batch_file[i,1],".trim.log ","| awk '{print $4}' | sed 's/,//g'"), intern = T)
        lm_workflow_stats[i,7] <- NA
        lm_workflow_stats[i,9] <- system(paste0("grep 'Reads written' ",out_dir,"/logs/",batch_file[i,1],".trim.log ","| awk '{print $5}' | sed 's/,//g'"), intern = T)
        lm_workflow_stats[i,11] <- system(paste0("sed 1d ",out_dir,"/hits/",batch_file[i,1],".mghits | wc -l"), intern = T)
      }
    
    }
    
    # Calculate fraction columns
    lm_workflow_stats$R1_adap_frac <- lm_workflow_stats$R1_adap / lm_workflow_stats$Total_Reads
    lm_workflow_stats$R2_adap_frac <- lm_workflow_stats$R2_adap / lm_workflow_stats$Total_Reads
    lm_workflow_stats$Good_Keep_Frac <- lm_workflow_stats$Good_Keep / lm_workflow_stats$Total_Reads
    lm_workflow_stats$Raw_Map_Frac <- lm_workflow_stats$Raw_Map / lm_workflow_stats$Total_Reads
    
    # Output data.frame
    lm_workflow_stats
    
  } # END LM BRANCH
    
}


## Function X

#### END FUNCTION DEFINE ####





#### BEGIN WORKFLOWS ####





###############
############### JUNCTION MAPPING WORKFLOW
###############


if (wf == "jm"){
  
  ### Create Correct Dir Structure
  
  out_dir <- paste0(out_dir,"/jm")
  dir.create(out_dir, recursive = T)
  
  
  ### Create Log File - WORKING!
  cat(
  paste0("ETmapper v0.03 Summary    Created: ", date()),"\n\n",
  "Program Parameters:\n",
  paste0("Workflow type is: ", wf),"\n",
  paste0("Total Samples: ",nrow(batch_file)),"\n",
  paste0("Adapter Trim DB: ", ad,"\n"),
  paste0("Model DB: ", md,"\n"),
  file = paste0(out_dir,"/summary.txt"))
  
  
  ### Determine if reads provided are se or pr end and initalize program
  if(ncol(batch_file) == 3){
    paired_end_data <- FALSE
    cat("\nInput Reads Identifed as Single End...Begining Analysis\n\n")
    cat(" Paired End = FALSE", file = paste0(out_dir,"/summary.txt"), append = T)
  }
  
  if(ncol(batch_file) == 4){
    paired_end_data <- TRUE
    cat("\nInput Data Identifed as Paired End...Begining Analysis\n\n")
    cat(" Paired End = TRUE", file = paste0(out_dir,"/summary.txt"), append = T)
  } 
  
  # If tests fail exit with error
  if (exists("paired_end_data") == FALSE){
    cat("ERROR: Batch file not formatted correctly")
    q(save="no")
  }
  
    
  

  #### PAIRED END BRANCH ####
  
  
  ### Begin Trimming / Mapping Steps of Junction Mapping Workflow (PAIRED END)
  if(paired_end_data == TRUE){

    ### Run adapter/flanking sequence trimming AND Quality Score Filtering
    for (i in 1:nrow(batch_file)){
    
      # Indicate what program is doing
      cat(paste0("\nTrimming Adapters/Flanks and Qscore Filter for: ",batch_file[i,1],"\n"))
    
      # run cutadapt command
      system(paste0("cutadapt -a file:",ad," -A file:",ad, # specify adapter types and file
                    " -j ",cpu," -O ",am," -q ",qs, # specify trimming params
                    " -o ",out_dir,"/",batch_file[i,3],".trim", # fwd output
                    " -p ",out_dir,"/",batch_file[i,4],".trim", # rev output
                    " ",rd,batch_file[i,3], # fwd read
                    " ",rd,batch_file[i,4], # rev read
                    " > ",out_dir,"/",batch_file[i,1],".trim.log")) # log files

    }
  
  
    ### Identifying and Trimming Models (ONLY 1 THREAD POSSIBLE w/ --info-file)
  
    # Run model finding for loop on fwd reads only but include reverse 
    for (i in 1:nrow(batch_file)){
    
      # Indicate what program is doing
      cat(paste0("\nFinding and Trimming Models for: ",batch_file[i,1],"\n"))
    
      # Run cutadapt command 
      system(paste0("cutadapt -g file:",md, # specify model file
                    " -j ",cpu," -O ",mm," -e ",et, # specify trimming params
                    " --discard-untrimmed", # discard read pairs without model
                    " --info-file ",out_dir,"/",batch_file[i,1],".info", # print 
                    " -o ",out_dir,"/",batch_file[i,3],".clean", # fwd output
                    " -p ",out_dir,"/",batch_file[i,4],".clean", # rev output
                    " ",out_dir,"/",batch_file[i,3],".trim", # fwd read
                    " ",out_dir,"/",batch_file[i,4],".trim", # rev read
                    " > ",out_dir,"/",batch_file[i,1],".clean.log")) # log files
    
    }
  
  
    ### Identifying barcodes and creating barcode db file
  
    # Identify flanking region to use as barcode query
    bc_flank <- system(paste0("grep -o \"NN.*.NN\" ",md," | sed 's/N//g'"), intern = T)
    bc_flank <- unique(bc_flank)
  
    for (i in 1:nrow(batch_file)){
  
      # Indicate what program is doing
      cat(paste0("\nFinding Barcodes for: ",batch_file[i,1],"\n"))
    
      # Filter cutadapt info file for reads with model
      system(paste0("awk '{if ($3!=-1) print}' ",out_dir,"/",batch_file[i,1],".info > ",out_dir,"/",batch_file[i,1],".info.filt"))
  
      # pull filtered model hit file
      ca_info <- fread(paste0(out_dir,"/",batch_file[i,1],".info.filt"), header = F)

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
      bc_out$Read <- sub(pattern = " .*$", replacement = "", bc_out$Read, perl = T) # kills anything after space in read name (i.e. 1:0:AGAAC, ect)
    
      # write barcodes out 
      write.table(bc_out, paste0(out_dir,"/",batch_file[i,1],".bc"), row.names = F, quote = F, sep = "\t")
  
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
                    " -o ",out_dir,"/",batch_file[i,3],".clean2", # fwd output
                    " -p ",out_dir,"/",batch_file[i,4],".clean2", # rev output
                    " ",out_dir,"/",batch_file[i,3],".clean", # fwd read
                    " ",out_dir,"/",batch_file[i,4],".clean", # rev read
                    " > ",out_dir,"/",batch_file[i,1],".clean2.log")) # log files
    
    }
  
  
    ### Run bowtie mapping 

    # Run bowtie mapping for loop for PAIRED END mapping
    for (i in 1:nrow(batch_file)){

      # Indicate what program is doing
      cat(paste0("\nMapping: ",batch_file[i,1],"\n"))

      # run bowtie mapping
      system(paste0("bowtie2 -x ",gd,"/bt2/All_Genomes", # specify genome database
                    " -p ",cpu," -X ",isl, # specify bowtie options
                    " -1 ",out_dir,"/",batch_file[i,3],".clean2", # specify fwd reads
                    " -2 ",out_dir,"/",batch_file[i,4],".clean2", # specify rev reads
                    " -S ",out_dir,"/",batch_file[i,1],".sam", # specify sam file output
                    " 2> ",out_dir,"/",batch_file[i,1],".bowtie.log")) # specify log file output
    
      # Convert SAM to BAM format with sorting and indexing
      system(paste0("samtools view -S -b ",out_dir,"/",batch_file[i,1],".sam > ",out_dir,"/",batch_file[i,1],".tmpbam; 
                    samtools sort ",out_dir,"/",batch_file[i,1],".tmpbam -o ",out_dir,"/",batch_file[i,1],".sorted.bam; 
                    samtools index ",out_dir,"/",batch_file[i,1],".sorted.bam"))
                  
      # Run read hit stats script
      system(paste0("python3 ",scripts,"/bam_pe_stats.py ",gd,"/scaff2bin.txt ",out_dir,"/",batch_file[i,1],".sorted.bam > ",out_dir,"/",batch_file[i,1],".tmphits"))
    
      # Integrate hit reads with barcodes into combined output and write
      hit_dat <- fread(paste0(out_dir,"/",batch_file[i,1],".tmphits"), header = T, stringsAsFactors = F)
      bc_dat <- fread(paste0(out_dir,"/",batch_file[i,1],".bc"), header = T, stringsAsFactors = F)
      merge_dat <- merge(hit_dat, bc_dat, by = "Read", all.x = T)
    
      write.table(merge_dat, paste0(out_dir,"/",batch_file[i,1],".hits"), row.names = F, quote = F, sep = "\t")
    
      cat(paste0("Finished reading a BAM, wrote output to ", paste0(batch_file[i,1],".hits")))
      cat("\nSneak peak of unfiltered results:\n")
      print(table(merge_dat[which(merge_dat$GENOME1 == merge_dat$GENOME2),]$GENOME1))
      cat("\n")
      
    }

    ### Clean up files and create out_dir structure
    clean.up()
    
    
    ### Build Metadata File
    jm_workflow_stats <- pull.run.stats()
    write.table(jm_workflow_stats, paste0(out_dir,"/jm_workflow_stats.txt"), row.names = F, quote = F, sep = "\t")
    
  } ### END Trimming / Mapping Steps of Junction Mapping Workflow (PAIRED END)
  
  
  
  
  #### SINGLE END BRANCH ####
  
  
  ### Begin Trimming / Mapping Steps of Junction Mapping Workflow (SINGLE END)
  if(paired_end_data == FALSE){
    
    
    ### Run adapter/flanking sequence trimming AND Quality Score Filtering
    for (i in 1:nrow(batch_file)){
      
      # Indicate what program is doing
      cat(paste0("\nTrimming Adapters/Flanks and Qscore Filter for: ",batch_file[i,1],"\n"))
      
      # run cutadapt command
      system(paste0("cutadapt -a file:",ad, # specify adapter types and file
                    " -j ",cpu," -O ",am," -q ",qs, # specify trimming params
                    " -o ",out_dir,"/",batch_file[i,3],".trim", # fwd output
                    " ",rd,batch_file[i,3], # fwd read
                    " > ",out_dir,"/",batch_file[i,1],".trim.log")) # log files
      
    }
  
    
    ### Run model finding for loop on fwd reads only AND Size Filter HERE (ONLY 1 THREAD POSSIBLE w/ --info-file)
    for (i in 1:nrow(batch_file)){
      
      # Indicate what program is doing
      cat(paste0("\nFinding and Trimming Models and Length Filter for: ",batch_file[i,1],"\n"))
      
      # Run cutadapt command 
      system(paste0("cutadapt -g file:",md, # specify model file
                    " -j ",cpu," -O ",mm," -e ",et, # specify trimming params
                    " --discard-untrimmed", # discard read pairs without model
                    " --minimum-length ",rl, # discard read if length < rl bp
                    " --info-file ",out_dir,"/",batch_file[i,1],".info", # print 
                    " -o ",out_dir,"/",batch_file[i,3],".clean", # fwd output
                    " ",out_dir,"/",batch_file[i,3],".trim", # fwd read
                    " > ",out_dir,"/",batch_file[i,1],".clean.log")) # log files
      
    }
  
  
    ### Identifying barcodes and creating barcode db file
    
    # Identify flanking region to use as barcode query
    bc_flank <- system(paste0("grep -o \"NN.*.NN\" ",md," | sed 's/N//g'"), intern = T)
    bc_flank <- unique(bc_flank)
    
    for (i in 1:nrow(batch_file)){
      
      # Indicate what program is doing
      cat(paste0("\nFinding Barcodes for: ",batch_file[i,1],"\n"))
      
      # Filter cutadapt info file for reads with model
      system(paste0("awk '{if ($3!=-1) print}' ",out_dir,"/",batch_file[i,1],".info > ",out_dir,"/",batch_file[i,1],".info.filt"))
      
      # pull filtered model hit file
      ca_info <- fread(paste0(out_dir,"/",batch_file[i,1],".info.filt"), header = F)
      
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
      write.table(bc_out, paste0(out_dir,"/",batch_file[i,1],".bc"), row.names = F, quote = F, sep = "\t")
      
    }  
  
    
    # Run bowtie mapping for loop for SINGLE END mapping
    for (i in 1:nrow(batch_file)){
      
      # Indicate what program is doing
      cat(paste0("\nMapping: ",batch_file[i,1],"\n"))
      
      # run bowtie mapping
      system(paste0("bowtie2 -x ",gd,"/bt2/All_Genomes", # specify genome database
                    " -p ",cpu, # specify bowtie options
                    " -U ",out_dir,"/",batch_file[i,3],".clean", # specify fwd reads
                    " -S ",out_dir,"/",batch_file[i,1],".sam", # specify sam file output
                    " 2> ",out_dir,"/",batch_file[i,1],".bowtie.log")) # specify log file output
      
      # Convert SAM to BAM format with sorting and indexing
      system(paste0("samtools view -S -b ",out_dir,"/",batch_file[i,1],".sam > ",out_dir,"/",batch_file[i,1],".tmpbam; 
                    samtools sort ",out_dir,"/",batch_file[i,1],".tmpbam -o ",out_dir,"/",batch_file[i,1],".sorted.bam; 
                    samtools index ",out_dir,"/",batch_file[i,1],".sorted.bam"))
      
      # Run read hit stats script
      system(paste0("python3 ",scripts,"/bam_se_stats.py ",gd,"/scaff2bin.txt ",out_dir,"/",batch_file[i,1],".sorted.bam > ",out_dir,"/",batch_file[i,1],".tmphits"))
      
      # Integrate hit reads with barcodes into combined output and write
      hit_dat <- fread(paste0(out_dir,"/",batch_file[i,1],".tmphits"), header = T, stringsAsFactors = F)
      bc_dat <- fread(paste0(out_dir,"/",batch_file[i,1],".bc"), header = T, stringsAsFactors = F)
      merge_dat <- merge(hit_dat, bc_dat, by = "Read", all.x = T)
      
      write.table(merge_dat, paste0(out_dir,"/",batch_file[i,1],".hits"), row.names = F, quote = F, sep = "\t")
    
      cat(paste0("Finished reading a BAM, wrote output to ", paste0(batch_file[i,1],".hits")))
      cat("\nSneak peak of unfiltered results:\n")
      print(table(merge_dat$GENOME))
      cat("\n")

    } # End bowtie mapping for loop for SINGLE END mapping
    
    ### Clean up files and create out_dir structure
    clean.up()

  } ### End Trimming / Mapping Steps of Junction Mapping Workflow (SINGLE END) 
  
  cat("Junction mapping workflow finished successfully :-)\n\n")

} ### End Junction Mapping Workflow





###############
############### LITE META-G WORKFLOW
###############


if (wf == "lm") {

  ### Create Correct Dir Structure
  
  out_dir <- paste0(out_dir,"/lm")
  dir.create(out_dir, recursive = T)
  
  
  ### Create Log File - WORKING!
  cat(
    paste0("ETmapper v0.03 Summary    Created: ", date()),"\n\n",
    "Program Parameters:\n",
    paste0("Workflow type is: ", wf),"\n",
    paste0("Total Samples: ",nrow(batch_file)),"\n",
    file = paste0(out_dir,"/summary.txt"))
  
  
  ### Determine if reads provided are se or pr end and initalize program
  if(ncol(batch_file) == 3){
    paired_end_data <- FALSE
    cat("\nInput Reads Identifed as Single End...Begining Analysis\n\n")
    cat(" Paired End = FALSE", file = paste0(out_dir,"/summary.txt"), append = T)
  }
  
  if(ncol(batch_file) == 4){
    paired_end_data <- TRUE
    cat("\nInput Data Identifed as Paired End...Begining Analysis\n\n")
    cat(" Paired End = TRUE", file = paste0(out_dir,"/summary.txt"), append = T)
  } 
  
  # If tests fail exit with error
  if (exists("paired_end_data") == FALSE){
    cat("ERROR: Batch file not formatted correctly")
    q(save="no")
  }
  
  
  
  
  #### PAIRED END BRANCH ####
  
  
  ### Begin Mapping for Lite Metagenomics Workflow (PAIRED END)
  if(paired_end_data == TRUE){
    
    ### Run adapter/flanking sequence trimming AND Quality Score Filtering
    for (i in 1:nrow(batch_file)){
      
      # Indicate what program is doing
      cat(paste0("\nTrimming Adapters, Qscore Filter, Size Filter for: ",batch_file[i,1],"\n"))
      
      # run cutadapt command
      system(paste0("cutadapt -a file:",ad," -A file:",ad, # specify adapter types and file
                    " -j ",cpu," -O ",am," -q ",qs, # specify trimming params
                    " --minimum-length ",rl, # discard read pairs if both reads are < rl bp
                    " --pair-filter=both", # At least 1 read must be >= rl
                    " -o ",out_dir,"/",batch_file[i,3],".trim", # fwd output
                    " -p ",out_dir,"/",batch_file[i,4],".trim", # rev output
                    " ",rd,batch_file[i,3], # fwd read
                    " ",rd,batch_file[i,4], # rev read
                    " > ",out_dir,"/",batch_file[i,1],".trim.log")) # log files
    
    }
    
    
    # Run bowtie mapping for loop for PAIRED END mapping
    for (i in 1:nrow(batch_file)){
      
      # Indicate what program is doing
      cat(paste0("\nMapping: ",batch_file[i,1],"\n"))
      
      # run bowtie mapping
      system(paste0("bowtie2 -x ",gd,"/bt2/All_Genomes", # specify genome database
                    " -p ",cpu," -X ",isl, # specify bowtie options
                    " -1 ",out_dir,"/",batch_file[i,3],".trim", # specify fwd reads
                    " -2 ",out_dir,"/",batch_file[i,4],".trim", # specify rev reads
                    " -S ",out_dir,"/",batch_file[i,1],".sam", # specify sam file output
                    " 2> ",out_dir,"/",batch_file[i,1],".bowtie.log")) # specify log file output
      
      # Convert SAM to BAM format with sorting and indexing
      system(paste0("samtools view -S -b ",out_dir,"/",batch_file[i,1],".sam > ",out_dir,"/",batch_file[i,1],".tmpbam; 
                    samtools sort ",out_dir,"/",batch_file[i,1],".tmpbam -o ",out_dir,"/",batch_file[i,1],".sorted.bam; 
                    samtools index ",out_dir,"/",batch_file[i,1],".sorted.bam"))
      
      # Run read hit stats script
      system(paste0("python3 ",scripts,"/bam_pe_stats.py ",gd,"/scaff2bin.txt ",out_dir,"/",batch_file[i,1],".sorted.bam > ",out_dir,"/",batch_file[i,1],".mghits"))
      
      # Write completed task
      cat(paste0("Finished reading a BAM, wrote output to ", paste0(batch_file[i,1],".mghits")))
      
    }
    
    ### Clean up files and create out_dir structure
    clean.up()
    
    
    ### Build Metadata File
    lm_workflow_stats <- pull.run.stats()
    write.table(lm_workflow_stats, paste0(out_dir,"/lm_workflow_stats.txt"), row.names = F, quote = F, sep = "\t")
    
  }
  
  cat("Lite metagenomics workflow finished successfully :-)\n\n")
  
}






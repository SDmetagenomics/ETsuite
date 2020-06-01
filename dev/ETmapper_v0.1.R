#!/usr/bin/env Rscript

### Set Version
etmapper_version <- c("v0.10")                                             ######**** Added VERSION NUMBER VARIABLE

### Load Test Data
# rd <- normalizePath("../Studies/20_1_21_BR_35_ET/et_reads/")
# batch_file <- fread("~/Desktop/test_batch_file.txt", sep = "\t", header = T)
#  ca_info <- fread("../Studies/20_2_27_BR_30_ET_KlebCurve_NT/ET_mapper_dual/jm/raw/BR_30_ET_1.info.filt")
#  md <- normalizePath("db/ETseq6_VcCas_mariner_v2.fa")
# # out_dir <-"~/Desktop"
# # bl <- 20
# # fe <- 1


### Check Libraries
if ("data.table" %in% installed.packages() == F){
  print("Please install R package data.table. Program quitting...")
  q(save="no")
}

if ("stringr" %in% installed.packages() == F){
  print("Please install R package stringr. Program quitting...")
  q(save="no")
}


### Load Libraries
library(data.table)
library(stringr)


### Set up path variables for associated scripts and databases

# Get relative path of ETmapper install directory
initial.options <- commandArgs(trailingOnly = F)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

# Set up default database paths (note: normalizePath creates a full path from relative path)    #### NEED TO SUPPRESS WARNIGNS IF THESE PATHS ARE NOT FOUND / PUT REASONABLE DEFAULTS SO THERE ARE NOT WARNINGS
ad <- normalizePath(paste0(script.basename,"/../db/ETseq6_adap.fa")) #adapter sequences
md <- normalizePath(paste0(script.basename,"/../db/ETseq_newprimers_mariner.fa")) #model sequences
scripts <- normalizePath(paste0(script.basename,"/../scripts")) #accesory scripts directory          ######**** CHANGED DEFAULT MODEL + ADAPTERS



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
    #######    #    #    # #    # #      #      ###### #    #  v0.04
      
    Usage: ETmapper.R -w [workflow] -d [read_dir] -b [batch_file] -g [genome_db] [options]

    Mandatory Arguments:
    
      -w: Workflow type (No Default)
          jm - Junction mapping
          lm - Metagenomics coverage
      -d: Directory containing read files (No Default)
      -b: Batch file with sample information (No Default)
      -g: Directory containing genome database (No Default)

    Adapter Filtering Options:
    
      -ad: Adapter sequence file (Default: db/ETseq6_adap.fa)
      -am: Min length of adapter match (Default: 5)
      -qs: Min base quality (Default: 20)

    Model/Barcode Identification Options:
    
      -md: Junction model sequence file (Default: db/ETseq_newprimers_mariner.fa)
      -mm: Min length of model match (Default: 25)
      -et: Model match error (Default: 0.02)
      -fe: Barcode flanking sequence match error (Default: 1)
      -bl: Barcode length (Default: 20)
      -rl: Min final read length (Default: 40)
      
    Read Mapping / Filtering Options:
      
      -N: Read length (Default: 150)
      -M: Min post trim fwd read size for PE mapping (Default: 20)
      -X: Maximum insert length (Default: 500)
      -Q: Min MapQ Score (Default: 20) 
      -E: Max number of mismatches (Default: 5; SE only)
      
    Program Control:
    
      -o: Output directory (Default: ET_Seq)
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
bf <- normalizePath(bf)
batch_file <- fread(bf, sep = "\t", header = T)                        ######**** BATCH FILE NOW HAS HEADER

# Genome Database Directory                                            ######**** ALL DATABASE / BATCH FILE PATHS ARE NOW NORMALIZED 
gd <- args[which(args == "-g") + 1]
gd <- normalizePath(gd)


## Adapter Filtering Arguments

# Adapter Database File
ad <- ad
if("-ad" %in% args){
  ad <- args[which(args == "-ad") + 1]
  ad <- normalizePath(ad)
}

# MMin Length of Adapter Match (Default: 5)
am <- 5
if("-am" %in% args){
  am <- as.numeric(args[which(args == "-am") + 1])
}

# Min base quality score (Default: 20)
qs <- 20
if("-qs" %in% args){
  qs <- as.numeric(args[which(args == "-qs") + 1])
}


## Model/Barcode Identification Options

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

# Barcode flanking sequence match error (Default: 1)
fe <- 1
if("-fe" %in% args){
  fe <- as.numeric(args[which(args == "-fe") + 1])
}

# Barcode length (Default: 20)
bl <- 20
if("-bl" %in% args){
  bl <- as.numeric(args[which(args == "-bl") + 1])
}

# Min final read length (Default: 40)
rl <- 40
if("-rl" %in% args){
  rl <- as.numeric(args[which(args == "-rl") + 1])
}


## Read Mapping Options

# Sequencing read length (Default: 150 bp)
sq_rd_len <- 150
if("-N" %in% args){
  sq_rd_len <- as.numeric(args[which(args == "-N") + 1])
}

# Min post trim fwd read size for PE mapping (Default: 20)
min_fwd_size <- 20
if("-M" %in% args){
  min_fwd_size <- as.numeric(args[which(args == "-M") + 1])
}

# Max insert length (Default: 500 bp)
isl <- 500
if("-X" %in% args){
  isl <- as.numeric(args[which(args == "-X") + 1])
}

# Min MapQ Score (Default: 20)
mq_cut <- 20
if("-Q" %in% args){
  mq_cut <- as.numeric(args[which(args == "-Q") + 1])                                            ######**** ADDED HIT FILTERING TO ETmapper SCRIPT 
}

# Max number of mismatches (Default: 5; SE only)
mm_cut <- 5
if("-E" %in% args){
  mm_cut <- as.numeric(args[which(args == "-E") + 1])
}


## Program Control Options

# number of cores to use where applicable
cpu <- 1
if("-cpu" %in% args){
  cpu <- as.numeric(args[which(args == "-cpu") + 1])
}

# Output directory (Will be created if not specified)
out_dir <- "ET_Seq"                                                                             ######**** CHANGED DEFAULT OUTPUT TO ET-SEQ 
if("-o" %in% args){
  out_dir <- args[which(args == "-o") + 1]
}
dir.create(out_dir, recursive = T)


## Set Data.Table Threads
setDTthreads(cpu)                                                                               ######**** SET Data.Table THREADS, CHANGING DATA FRAMES TO DATA.TABLES EVERYWHERE


###### HOW TO FIX NON-ARGUMENT ENTRY BUG i.e. ETmapper.R -w [workflow] -d -b [batch_file]
# for important arguments check that: 
# 
# which(argument) + 1 != c(all_other_arguments)
# 
# If somebody forgets to enter something then the thing +1 from the argument in the vector will
# be another argument 




#### BEGIN FUNCTION DEFINE ####


## Function 1: Identify Barcodes
# Write: bc_out
find.bc <- function(){
  
  ### TESTING:
  # ca_info <- fread("../Studies/20_2_27_BR_30_ET_KlebCurve_NT/ET_mapper_dual/jm/raw/BR_30_ET_1.info.filt")
  # md <- normalizePath("db/ETseq6_VcCas_mariner_v2.fa")
  # # out_dir <-"~/Desktop"
  #  bl <- 20
  #  fe <- 1
  
  
  # Filter cutadapt info file for reads with model
  system(paste0("awk '{if ($3!=-1) print}' ",out_dir,"/",batch_file$SAMPLE[i],".info > ",out_dir,"/",batch_file$SAMPLE[i],".info.filt"))
    
  # pull filtered model hit file
  ca_info <- fread(paste0(out_dir,"/",batch_file$SAMPLE[i],".info.filt"), header = F)
    
  # create barcode output file
  bc_out <- data.table(Read = ca_info$V1, model = ca_info$V8, mod_len = ca_info$V4)
    
  # create tmp barcode storage data frame
  bc_tmp_store <- data.table()
    
  # loop over each possible barcode flanking sequences
  for (j in 1:length(bc_flank)){
      
    # Store length of flank sequence j
    flank_length <- nchar(bc_flank[j])
      
    # parse ca_info for those hits with aproximate flank sequence match within error range
    bc_flank_i_matches <- ca_info[agrep(bc_flank[j], ca_info$V6, max = list(ins = 0, del = 0, sub = fe), ignore.case = T),c(1,6)]
      
    # skip iteration if no flank sequence match
    if(nrow(bc_flank_i_matches) == 0){
      next
    }
      
    # Use aproximate regex to identify the flank sequence and 20bp upstream (barcode) or whatever bl is
    matches <- aregexec(paste0(bc_flank[j],".{0,",bl,"}"), bc_flank_i_matches$V6, max = list(ins = 0, del = 0, sub = fe), ignore.case = T)
      
    # Create list of sequences aprox matching flank + bl
    match_seqs <- regmatches(bc_flank_i_matches$V6, matches)
      
    # Parse matches for barcodes and add to data frame
    bc_flank_i_matches$barcodes <- as.character(lapply(match_seqs, function(x) str_sub(x, start = -bl)))
      
    # Parse matches for flank sequence and add to data frame
    bc_flank_i_matches$flank_seq <- as.character(lapply(match_seqs, function(x) str_sub(x, end = flank_length)))
      
    # Create final data frame w/ column for flank and match and cbind to outupt
    bc_flank_i_matches <- data.table(bc_flank_i_matches[,-2])
    bc_tmp_store <- rbind(bc_tmp_store, bc_flank_i_matches)
      
  }
    
  # merge barcodes to output file and remove sequencing barcode line
  bc_out <- merge(bc_out, bc_tmp_store, by.x = "Read", by.y = "V1", all.x = T)
  bc_out$Read <- sub(pattern = " .*$", replacement = "", bc_out$Read, perl = T) # kills anything after space in read name (i.e. 1:0:AGAAC, ect)
    
  # write barcodes out 
  fwrite(bc_out, paste0(out_dir,"/",batch_file$SAMPLE[i],".bc"), row.names = F, quote = F, sep = "\t")

  # return bc_flank
  return(bc_flank)
  
}


## Function 2: Clean up files and make output structure
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


## Function 3: Pull run statistics from log files
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
    jm_workflow_stats <- data.table(batch_file,
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
      cat(paste0("\nAggregating ETmapper stats for: ",batch_file$SAMPLE[i],"\n"))
      
      # Pull stats from logs
      if(paired_end_data == TRUE){
        jm_workflow_stats[i,5] <- as.numeric(system(paste0("grep 'Total read pairs processed:' ",out_dir,"/logs/",batch_file$SAMPLE[i],".trim.log ","| awk '{print $5}' | sed 's/,//g'"), intern = T))
        jm_workflow_stats[i,6] <- as.numeric(system(paste0("grep 'Read 1 with adapter:' ",out_dir,"/logs/",batch_file$SAMPLE[i],".trim.log ","| awk '{print $5}' | sed 's/,//g'"), intern = T))
        jm_workflow_stats[i,8] <- as.numeric(system(paste0("grep 'Read 2 with adapter:' ",out_dir,"/logs/",batch_file$SAMPLE[i],".trim.log ","| awk '{print $5}' | sed 's/,//g'"), intern = T))
        jm_workflow_stats[i,10] <- count.primer()
        jm_workflow_stats[i,12] <- as.numeric(system(paste0("grep 'Pairs written' ",out_dir,"/logs/",batch_file$SAMPLE[i],".clean.log ","| awk '{print $5}' | sed 's/,//g'"), intern = T))
        jm_workflow_stats[i,14] <- as.numeric(system(paste0("grep 'Read 2 with adapter:' ",out_dir,"/logs/",batch_file$SAMPLE[i],".clean2.log ","| awk '{print $5}' | sed 's/,//g'"), intern = T))
        jm_workflow_stats[i,16] <- as.numeric(system(paste0("grep 'Pairs written' ",out_dir,"/logs/",batch_file$SAMPLE[i],".clean2.log ","| awk '{print $5}' | sed 's/,//g'"), intern = T))
        jm_workflow_stats[i,18] <- as.numeric(system(paste0("sed 1d ",out_dir,"/hits/",batch_file$SAMPLE[i],".hits | wc -l"), intern = T))
      }

      if(paired_end_data == FALSE){
        jm_workflow_stats[i,4] <- as.numeric(system(paste0("grep 'Total reads processed:' ",out_dir,"/logs/",batch_file$SAMPLE[i],".trim.log ","| awk '{print $4}' | sed 's/,//g'"), intern = T)) # Total Reads
        jm_workflow_stats[i,5] <- as.numeric(system(paste0("grep 'Reads with adapters:' ",out_dir,"/logs/",batch_file$SAMPLE[i],".trim.log ","| awk '{print $4}' | sed 's/,//g'"), intern = T)) # R1_adap
        jm_workflow_stats[i,7] <- NA # R2_adap
        jm_workflow_stats[i,9] <- count.primer() #Tn_Primer
        jm_workflow_stats[i,11] <- as.numeric(system(paste0("grep 'Reads with adapters:' ",out_dir,"/logs/",batch_file$SAMPLE[i],".clean.log ","| awk '{print $4}' | sed 's/,//g'"), intern = T)) #Model_Keep
        jm_workflow_stats[i,13] <- NA # R2_model
        jm_workflow_stats[i,15] <- as.numeric(system(paste0("grep 'Reads written' ",out_dir,"/logs/",batch_file$SAMPLE[i],".clean.log ","| awk '{print $5}' | sed 's/,//g'"), intern = T)) #Good_Keep
        jm_workflow_stats[i,17] <- as.numeric(system(paste0("sed 1d ",out_dir,"/hits/",batch_file$SAMPLE[i],".hits | wc -l"), intern = T)) #Raw_map 
      }
      
    }
    
    # Calculate fraction columns
    jm_workflow_stats$R1_adap_frac <- jm_workflow_stats$R1_adap / jm_workflow_stats$Total_Reads
    jm_workflow_stats$R2_adap_frac <- jm_workflow_stats$R2_adap / jm_workflow_stats$Total_Reads
    jm_workflow_stats$Tn_Primer_Frac <- jm_workflow_stats$Tn_Primer / jm_workflow_stats$Total_Reads
    jm_workflow_stats$Model_Keep_Frac <- jm_workflow_stats$Model_Keep / jm_workflow_stats$Total_Reads
    jm_workflow_stats$R2_Model_Frac <- jm_workflow_stats$R2_Model / jm_workflow_stats$Model_Keep
    jm_workflow_stats$Good_Keep_Frac <- jm_workflow_stats$Good_Keep / jm_workflow_stats$Total_Reads
    jm_workflow_stats$Raw_Map_Frac <- jm_workflow_stats$Raw_Map / jm_workflow_stats$Total_Reads
    
    # Output data.table
    return(jm_workflow_stats)
    
  } # END JM BRANCH
  
  
  # Pull stats lm workflow
  if (wf == "lm"){
    
    # Create dataframe to hold output
    lm_workflow_stats <- data.table(batch_file,
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
      cat(paste0("\nAggregating ETmapper stats for: ",batch_file$SAMPLE[i],"\n"))
      
      # Pull stats from logs
      if(paired_end_data == TRUE){
        lm_workflow_stats[i,5] <- as.numeric(system(paste0("grep 'Total read pairs processed:' ",out_dir,"/logs/",batch_file$SAMPLE[i],".trim.log ","| awk '{print $5}' | sed 's/,//g'"), intern = T))
        lm_workflow_stats[i,6] <- as.numeric(system(paste0("grep 'Read 1 with adapter:' ",out_dir,"/logs/",batch_file$SAMPLE[i],".trim.log ","| awk '{print $5}' | sed 's/,//g'"), intern = T))
        lm_workflow_stats[i,8] <- as.numeric(system(paste0("grep 'Read 2 with adapter:' ",out_dir,"/logs/",batch_file$SAMPLE[i],".trim.log ","| awk '{print $5}' | sed 's/,//g'"), intern = T))
        lm_workflow_stats[i,10] <- as.numeric(system(paste0("grep 'Pairs written' ",out_dir,"/logs/",batch_file$SAMPLE[i],".trim.log ","| awk '{print $5}' | sed 's/,//g'"), intern = T))
        lm_workflow_stats[i,12] <- as.numeric(system(paste0("sed 1d ",out_dir,"/hits/",batch_file$SAMPLE[i],".mghits | wc -l"), intern = T))
      }
      
      ### STILL UNTESTED
      if(paired_end_data == FALSE){
        lm_workflow_stats[i,4] <- system(paste0("grep 'Total reads processed:' ",out_dir,"/logs/",batch_file$SAMPLE[i],".trim.log ","| awk '{print $4}' | sed 's/,//g'"), intern = T)
        lm_workflow_stats[i,5] <- system(paste0("grep 'Reads with adapters:' ",out_dir,"/logs/",batch_file$SAMPLE[i],".trim.log ","| awk '{print $4}' | sed 's/,//g'"), intern = T)
        lm_workflow_stats[i,7] <- NA
        lm_workflow_stats[i,9] <- system(paste0("grep 'Reads written' ",out_dir,"/logs/",batch_file$SAMPLE[i],".trim.log ","| awk '{print $5}' | sed 's/,//g'"), intern = T)
        lm_workflow_stats[i,11] <- system(paste0("sed 1d ",out_dir,"/hits/",batch_file$SAMPLE[i],".mghits | wc -l"), intern = T)
      }
    
    }
    
    # Calculate fraction columns
    lm_workflow_stats$R1_adap_frac <- lm_workflow_stats$R1_adap / lm_workflow_stats$Total_Reads
    lm_workflow_stats$R2_adap_frac <- lm_workflow_stats$R2_adap / lm_workflow_stats$Total_Reads
    lm_workflow_stats$Good_Keep_Frac <- lm_workflow_stats$Good_Keep / lm_workflow_stats$Total_Reads
    lm_workflow_stats$Raw_Map_Frac <- lm_workflow_stats$Raw_Map / lm_workflow_stats$Total_Reads
    
    # Output data.table
    return(lm_workflow_stats)
    
  } # END LM BRANCH
    
}


## Function 4: Check post model trim read length for given read length
check.mod.short <- function(){
  
  ## Get max model length
  mod_lengths <- as.numeric(system(paste0("grep -v '^>' ",
                                          md,
                                          " | awk '{ print length }'"), intern = T))
  max_mod_length <- max(mod_lengths)
  
  ## Compare to sequencing read length
  max_fwd_size <- sq_rd_len - max_mod_length
  fwd_read_short <- max_fwd_size < min_fwd_size
  
  # Return output
  return(fwd_read_short)
  
}



#### END FUNCTION DEFINE ####





#### BEGIN WORKFLOWS ####





###############
############### JUNCTION MAPPING WORKFLOW
###############


if (wf == "jm"){
  
  
  ### Create Correct Dir Structure
  out_dir <- paste0(out_dir,"/jm")
  dir.create(out_dir, recursive = T)
  
  
  
  
  ### Create Log File
  cat(
  paste0("ETmapper v0.10 Summary    Created: ", date()),"\n\n",
  "Program Parameters:\n",
  paste0("Arguments: ", args),"\n",
  paste0("Workflow type is: ", wf),"\n",
  paste0("Total Samples: ",nrow(batch_file)),"\n",
  paste0("Batch File: ", bf,"\n"),                     ######**** LOG FILE NOW GIVES PATH OF BATCH FILE FOR FUTURE METADATA USE
  paste0("Adapter Trim File: ", ad,"\n"),
  paste0("Model File: ", md,"\n"),
  paste0("Genome Database: ", gd,"\n"),
  file = paste0(out_dir,"/run_log.txt"))         ######**** REMOVED MENTION OF PAIRED END INPUT DATA IN LOGS 
  
  
  
  
  ### Check batch file and start program
  
  ## If same number of read files for R1 and R2 start
  if((length(batch_file$R1) == length(batch_file$R2)) == T){
    cat("\nInput Data Correctly Identifed...Begining Analysis\n\n")
  } 
  
  ## If not same number indicate problem and quit
  if((length(batch_file$R1) == length(batch_file$R2)) == F){
    cat("ERROR: Batch file not formatted correctly")
    q(save="no")
  } 
  
  
  
  
  ### Identify Barcode Flanking Region from Model         ######**** MADE BC FLANK IDENTIFICATIION SOMETHING THAT IS DONE RIGHT AWAY 
  
  ## Identify flanking region to use as barcode query
  bc_flank <- system(paste0("grep -o \"^NN.*.NN\" ",md," | sed 's/N//g'"), intern = T)
  bc_flank <- unique(bc_flank)
  
  ## Identify flanks that are < fe+1 distance from each other and only keep 1 of those at that similarity level
  if (length(bc_flank) >= 2){
    bc_dist <- hclust(as.dist(adist(bc_flank)))
    bc_cut <- cutree(bc_dist, h = fe+1)
    bc_flank <- bc_flank[!duplicated(bc_cut)]
  }
  
  # test bc flank output
  cat(bc_flank)
  
  
  
  
  ### Check if fwd reads will be too short for PE mapping after trimming model and Log
  fwd_read_short <- check.mod.short()
  
  if(fwd_read_short == F){
    cat(paste0("Mapping Mode: PE\n"),
        file = paste0(out_dir,"/run_log.txt"),
        append = T)
  }
  
  if(fwd_read_short == T){
    cat(paste0("Mapping Mode: SE\n"),
        file = paste0(out_dir,"/run_log.txt"),
        append = T)
  }
  
  
  
  
  #### BEGIN MAIN PROGRAM LOOP
  
  for (i in 1:nrow(batch_file)){                  ######**** CHANGED PROGRAM SO NOW EACH FILE IS PROCESSED SEQUENTALLY 
    
    
  ### 1) Indicate Sample program is working on
    
    cat(paste0("\nNow Processing Sample: ",batch_file$SAMPLE[i],"\n"))
    
    
    
    
  ### 2) Run adapter/flanking sequence trimming AND Quality Score Filtering
    
    ## Indicate what program is doing
    cat("Trimming Adapters/Flanks and Qscore Filter...\n")
    
    ## run cutadapt command
    system(paste0("cutadapt -a file:",ad," -A file:",ad, # specify adapter types and file
                  " -j ",cpu," -O ",am," -q ",qs, # specify trimming params
                  " -o ",out_dir,"/",batch_file$SAMPLE[i],".R1.trim", # fwd output        ######**** CHANGED ALL INPUTS TO BE READ FROM batch_file$R1 and $R2, ALSO ALL OUTPUTS FROM TRIMMING NOW LEAD WITH SAMPLE NAME
                  " -p ",out_dir,"/",batch_file$SAMPLE[i],".R2.trim", # rev output
                  " ",rd,batch_file$R1[i], # fwd read
                  " ",rd,batch_file$R2[i], # rev read
                  " > ",out_dir,"/",batch_file$SAMPLE[i],".trim.log")) # log files

  
  
    
  ### 3) Identifying and Trimming Transposon Models in Fwd Reads only
    
    ## Indicate what program is doing
    cat("Identifying and Trimming Tn Models...\n")

    ## Run cutadapt command 
    system(paste0("cutadapt -g file:",md, # specify model file
                  " -j ",cpu," -O ",mm," -e ",et, # specify trimming params
                  " --discard-untrimmed", # discard read pairs without model
                  " --info-file ",out_dir,"/",batch_file$SAMPLE[i],".info", # print 
                  " -o ",out_dir,"/",batch_file$SAMPLE[i],".R1.clean", # fwd output
                  " -p ",out_dir,"/",batch_file$SAMPLE[i],".R2.clean", # rev output
                  " ",out_dir,"/",batch_file$SAMPLE[i],".R1.trim", # fwd read
                  " ",out_dir,"/",batch_file$SAMPLE[i],".R2.trim", # rev read
                  " > ",out_dir,"/",batch_file$SAMPLE[i],".clean.log")) # log files
    
    
    
    
  ### 4) Identifying barcodes and creating barcode db file
    
    ## Indicate what program is doing
    cat("Identifying Barcodes in Each Read...\n")
    
    ## Find barcodes function
    find.bc()
    
    
    
    
  ### 5) Trimming Rev reads and filtering for min length in at least 1 read of a pair

    ## Indicate what program is doing
    cat("Trimming Tn Models from Rev Reads and Filtering for Length...\n")
    
    # run cutadapt command
    system(paste0("cutadapt -A file:",md, # specify model file
                  " -O 10"," -j ",cpu," -e ",et, # specify trimming params
                  " --minimum-length ",rl, # discard read pairs if both reads are < rl bp
                  " --pair-filter=both", # At least 1 read must be >= rl
                  " -o ",out_dir,"/",batch_file$SAMPLE[i],".R1.clean2", # fwd output
                  " -p ",out_dir,"/",batch_file$SAMPLE[i],".R2.clean2", # rev output
                  " ",out_dir,"/",batch_file$SAMPLE[i],".R1.clean", # fwd read
                  " ",out_dir,"/",batch_file$SAMPLE[i],".R2.clean", # rev read
                  " > ",out_dir,"/",batch_file$SAMPLE[i],".clean2.log")) # log files
    
  
  
    
  ### 6) Run bowtie mapping 

    ## Indicate what program is doing
    cat("Performing mapping of junctions to genomes...\n")

      
    ## Perform PE or SE mapping based on if trimming model makes fwd reads too short

    ## IF fwd read length after model trim long enough
    if (fwd_read_short == FALSE){
        
      # run bowtie mapping using fwd and rev reads
      system(paste0("bowtie2 -x ",gd,"/bt2/All_Genomes", # specify genome database
                    " -p ",cpu," -X ",isl, # specify bowtie options
                    " -1 ",out_dir,"/",batch_file$SAMPLE[i],".R1.clean2", # specify fwd reads
                    " -2 ",out_dir,"/",batch_file$SAMPLE[i],".R2.clean2", # specify rev reads
                    " -S ",out_dir,"/",batch_file$SAMPLE[i],".sam", # specify sam file output
                    " 2> ",out_dir,"/",batch_file$SAMPLE[i],".bowtie.log")) # specify log file output
    
    }  
      
    ## IF fwd read length after model trim too short
    if (fwd_read_short == TRUE){
        
      # run bowtime mapping using rev reads only
      system(paste0("bowtie2 -x ",gd,"/bt2/All_Genomes", # specify genome database
                    " -p ",cpu, # specify bowtie options
                    " -U ",out_dir,"/",batch_file$SAMPLE[i],".R2.clean2", # specify rev reads
                    " -S ",out_dir,"/",batch_file$SAMPLE[i],".sam", # specify sam file output
                    " 2> ",out_dir,"/",batch_file$SAMPLE[i],".bowtie.log")) # specify log file output
      
    }
      
      
    
    
  ### 7) Convert SAM 2 BAM format
      
    ## Indicate what program is doing
    cat("Converting SAM 2 BAM...\n")
        
    # Convert SAM to BAM format with sorting and indexing
    system(paste0("samtools view -S -b --threads ",cpu," ",out_dir,"/",batch_file$SAMPLE[i],".sam > ",out_dir,"/",batch_file$SAMPLE[i],".tmpbam; 
                  samtools sort --threads ",cpu," ",out_dir,"/",batch_file$SAMPLE[i],".tmpbam -o ",out_dir,"/",batch_file$SAMPLE[i],".sorted.bam; 
                  samtools index ",out_dir,"/",batch_file$SAMPLE[i],".sorted.bam"))
    
    # Remove SAM file  
    system(paste0("rm ",out_dir,"/",batch_file$SAMPLE[i],".sam"))
      
    
    
    
  ### 8) Run read hit stats script for pe or se mapping output
    
    ## Indicate what program is doing
    cat("Generating Hit Table...\n")
    
    # If fwd read long enough and PE mapping done use bam_pe_stats.py
    if (fwd_read_short == FALSE){
      system(paste0("python3 ",scripts,"/bam_pe_stats.py ",gd,"/scaff2bin.txt ",out_dir,"/",batch_file$SAMPLE[i],".sorted.bam > ",out_dir,"/",batch_file$SAMPLE[i],".tmphits"))
    }
      
    # If fwd read too short and reverse read (SE) mapping done use bam_se_stats.py
    if (fwd_read_short == TRUE){
      system(paste0("python3 ",scripts,"/bam_se_stats.py ",gd,"/scaff2bin.txt ",out_dir,"/",batch_file$SAMPLE[i],".sorted.bam > ",out_dir,"/",batch_file$SAMPLE[i],".tmphits"))
    }
      
    # Integrate hit reads with barcodes into combined output
    hit_dat <- fread(paste0(out_dir,"/",batch_file$SAMPLE[i],".tmphits"), header = T, stringsAsFactors = F)
    bc_dat <- fread(paste0(out_dir,"/",batch_file$SAMPLE[i],".bc"), header = T, stringsAsFactors = F)
    merge_dat <- merge.data.table(hit_dat, bc_dat, by = "Read", all.x = T)
    merge_dat$barcodes <- ifelse(merge_dat$barcodes == "", NA, merge_dat$barcodes)
    
    
    
  ### 9) Filter Output Hits
    
    ## Indicate what program is doing
    cat("Filtering Hit Table...\n")
    
    ## Filtering if PE mapping was done
    if (fwd_read_short == FALSE){
      
      # Filter mappings based on genome and mapQ score
      merge_dat_filt <- subset(merge_dat, GENOME1 == GENOME2 & (MAPQ1 >= mq_cut | MAPQ2 >= mq_cut))
      
      # Filter out NA barcodes
      merge_dat_filt <- subset(merge_dat_filt, is.na(barcodes) == FALSE)
      
      # Add Tn Junction Postion
      merge_dat_filt$TNjunc <- ifelse(merge_dat_filt$STRAND1 == "+", merge_dat_filt$START1, merge_dat_filt$END1)
      
      # Remove GENOME2 Column and Change GNEOME1 to GENOME
      merge_dat_filt[,GENOME2:=NULL]
      setnames(merge_dat_filt, "GENOME1", "GENOME")
      
      # Calculate Losses from Filtering
      f_loss <- nrow(merge_dat_filt) / nrow(merge_dat)
      cat(paste0(round(f_loss * 100, digits = 2), "% Hits remaining after filtering...\n"))
      
    }
    
    
    ## Filtering if SE mapping was done
    if (fwd_read_short == TRUE){
      
      # Filter mappings based on genome and mapQ score
      merge_dat_filt <- subset(merge_dat, MAPQ >= mq_cut & NM <= mm_cut)
      
      # Filter out NA barcodes
      merge_dat_filt <- subset(merge_dat_filt, is.na(barcodes) == FALSE)
      
      # Add Tn Junction Postion
      merge_dat_filt$TNjunc <- ifelse(merge_dat_filt$STRAND == "+", merge_dat_filt$START, merge_dat_filt$END)
      
      # Calculate Losses from Filtering
      f_loss <- nrow(merge_dat_filt) / nrow(merge_dat)
      cat(paste0(round(f_loss * 100, digits = 2), "% Hits remaining after filtering...\n"))
      
    }
    
    
    
  ### 10) Write Output Data / Clean Shit
    
    # Write hit table out
    fwrite(merge_dat_filt, paste0(out_dir,"/",batch_file$SAMPLE[i],".hits"), row.names = F, quote = F, sep = "\t")
    
  } ### END OF MAIN PROGRAM FOR LOOP
    
  
  
    
  ### Clean up files and create out_dir structure
  #clean.up()
  
    
  ### Build Metadata File
  #jm_workflow_stats <- pull.run.stats()
  #write.table(jm_workflow_stats, paste0(out_dir,"/jm_workflow_stats.txt"), row.names = F, quote = F, sep = "\t")
    
} ### END Trimming / Mapping Steps of Junction Mapping Workflow (PAIRED END)
  
cat("\nJunction mapping workflow finished successfully :-)\n\n")  
  
                                                                                        ######**** REMOVED SINGLE END MAPPING WORKFLOW

  


  
  
  
### End Junction Mapping Workflow





###############
############### LITE META-G WORKFLOW
###############


if (wf == "lm") {

  ### Create Correct Dir Structure
  
  out_dir <- paste0(out_dir,"/lm")
  dir.create(out_dir, recursive = T)
  
  
  ### Create Log File - WORKING!
  cat(
    paste0("ETmapper v0.04 Summary    Created: ", date()),"\n\n",
    "Program Parameters:\n",
    paste0("Workflow type is: ", wf),"\n",
    paste0("Total Samples: ",nrow(batch_file)),"\n",
    paste0("Adapter Trim File: ", ad,"\n"),
    paste0("Genome Database: ", gd,"\n"),
    file = paste0(out_dir,"/run_log.txt"))
  
  
  ### Determine if reads provided are se or pr end and initalize program
  if(ncol(batch_file) == 3){
    paired_end_data <- FALSE
    cat("\nInput Reads Identifed as Single End...Begining Analysis\n\n")
    cat(" Paired End = FALSE", file = paste0(out_dir,"/run_log.txt"), append = T)
  }
  
  if(ncol(batch_file) == 4){
    paired_end_data <- TRUE
    cat("\nInput Data Identifed as Paired End...Begining Analysis\n\n")
    cat(" Paired End = TRUE", file = paste0(out_dir,"/run_log.txt"), append = T)
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
      cat(paste0("\nTrimming Adapters, Qscore Filter, Size Filter for: ",batch_file$SAMPLE[i],"\n"))
      
      # run cutadapt command
      system(paste0("cutadapt -a file:",ad," -A file:",ad, # specify adapter types and file
                    " -j ",cpu," -O ",am," -q ",qs, # specify trimming params
                    " --minimum-length ",rl, # discard read pairs if both reads are < rl bp
                    " --pair-filter=both", # At least 1 read must be >= rl
                    " -o ",out_dir,"/",batch_file[i,3],".trim", # fwd output
                    " -p ",out_dir,"/",batch_file[i,4],".trim", # rev output
                    " ",rd,batch_file[i,3], # fwd read
                    " ",rd,batch_file[i,4], # rev read
                    " > ",out_dir,"/",batch_file$SAMPLE[i],".trim.log")) # log files
    
    }
    
    
    # Run bowtie mapping for loop for PAIRED END mapping
    for (i in 1:nrow(batch_file)){
      
      # Indicate what program is doing
      cat(paste0("\nMapping: ",batch_file$SAMPLE[i],"\n"))
      
      # run bowtie mapping
      system(paste0("bowtie2 -x ",gd,"/bt2/All_Genomes", # specify genome database
                    " -p ",cpu," -X ",isl, # specify bowtie options
                    " -1 ",out_dir,"/",batch_file[i,3],".trim", # specify fwd reads
                    " -2 ",out_dir,"/",batch_file[i,4],".trim", # specify rev reads
                    " -S ",out_dir,"/",batch_file$SAMPLE[i],".sam", # specify sam file output
                    " 2> ",out_dir,"/",batch_file$SAMPLE[i],".bowtie.log")) # specify log file output
      
      # Convert SAM to BAM format with sorting and indexing
      system(paste0("samtools view -S -b ",out_dir,"/",batch_file$SAMPLE[i],".sam > ",out_dir,"/",batch_file$SAMPLE[i],".tmpbam; 
                    samtools sort ",out_dir,"/",batch_file$SAMPLE[i],".tmpbam -o ",out_dir,"/",batch_file$SAMPLE[i],".sorted.bam; 
                    samtools index ",out_dir,"/",batch_file$SAMPLE[i],".sorted.bam"))
      
      # Run read hit stats script
      system(paste0("python3 ",scripts,"/bam_pe_stats.py ",gd,"/scaff2bin.txt ",out_dir,"/",batch_file$SAMPLE[i],".sorted.bam > ",out_dir,"/",batch_file$SAMPLE[i],".mghits"))
      
      # Write completed task
      cat(paste0("Finished reading a BAM, wrote output to ", paste0(batch_file$SAMPLE[i],".mghits")))
      
    }
    
    ### Clean up files and create out_dir structure
    clean.up()
    
    ### Build Metadata File
    lm_workflow_stats <- pull.run.stats()
    write.table(lm_workflow_stats, paste0(out_dir,"/lm_workflow_stats.txt"), row.names = F, quote = F, sep = "\t")
    
  }
  
  
  ### Begin Mapping for Lite Metagenomics Workflow (PAIRED END)
  if(paired_end_data == FALSE){
    
    ### Run adapter/flanking sequence trimming AND Quality Score Filtering
    for (i in 1:nrow(batch_file)){
      
      # Indicate what program is doing
      cat(paste0("\nTrimming Adapters, Qscore Filter, Size Filter for: ",batch_file$SAMPLE[i],"\n"))
      
      # run cutadapt command
      system(paste0("cutadapt -a file:",ad," -A file:",ad, # specify adapter types and file
                    " -j ",cpu," -O ",am," -q ",qs, # specify trimming params
                    " --minimum-length ",rl, # discard read pairs if both reads are < rl bp
                    " --pair-filter=both", # At least 1 read must be >= rl
                    " -o ",out_dir,"/",batch_file[i,3],".trim", # fwd output
                    " -p ",out_dir,"/",batch_file[i,4],".trim", # rev output
                    " ",rd,batch_file[i,3], # fwd read
                    " ",rd,batch_file[i,4], # rev read
                    " > ",out_dir,"/",batch_file$SAMPLE[i],".trim.log")) # log files
      
    }
    
    ### Clean up files and create out_dir structure
    clean.up()
    
    ### Build Metadata File
    lm_workflow_stats <- pull.run.stats()
    write.table(lm_workflow_stats, paste0(out_dir,"/lm_workflow_stats.txt"), row.names = F, quote = F, sep = "\t")
    
  }
  
  cat("Lite metagenomics workflow finished successfully :-)\n\n")
  
}






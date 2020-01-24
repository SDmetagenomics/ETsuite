#!/usr/bin/env Rscript

### Set Version
etstats_version <- c("v0.03")

### Check and Load Libraries

# Check data.table
if ("data.table" %in% installed.packages() == FALSE){
  print("Please install R package data.table. Program quitting...")
  q(save="no")
}

# Check plyr
if ("dplyr" %in% installed.packages() == FALSE){
  print("Please install R package dplyr. Program quitting...")
  q(save="no")
}

# Check ggplot2
if ("ggplot2" %in% installed.packages() == FALSE){
  print("Please install R package ggplot2. Program quitting...")
  q(save="no")
}

# # Check ggforce
# if ("ggforce" %in% installed.packages() == FALSE){
#   print("Please install R package ggforce. Program quitting...")
#   q(save="no")
# }


### Load packages
suppressMessages(library(data.table, quietly = T))
suppressMessages(library(dplyr, quietly = T))
suppressMessages(library(ggplot2, quietly = T))
#suppressMessages(library(ggforce, quietly = T))


### Set up path variables for associated scripts and databases

# Get relative path of ETmapper install directory
initial.options <- commandArgs(trailingOnly = F)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

# Set up default database paths (note: normalizePath creates a full path from relative path)
scripts <- normalizePath(paste0(script.basename,"/../scripts")) #accesory scripts directory



### Collect and parse arguments
args <- commandArgs(trailingOnly = T)


### TESTING STUFF REMOVE LATER ###
print(args)
cat(paste0("\nTotal Arguments: ",length(args),"\n"))


## NOT IN Operator for Arg Parseing
'%notin%' <- Negate('%in%')

## Display help if no args or -h
if("-h" %in% args | !("-w" %in% args) | !("-d" %in% args) | length(args) == 0) {
  cat("
      ETstats v0.03
      
      Usage: ETstats.R -w [workflow] -d [Data_dir] [options]

      Mandatory Arguments:
      
      -w: Workflow type (No Default)
          ss - basic stats summary
          cs - comparative stats summary
          rs - read statistics
          ip - insertion plot
      -d: Directory with data for analysis (No Default; can be ETmapper or ETstats)
      
      Optional Arguments:
      
      -s: Sample info (No Default; This file format will depend on the workflow. See documentation)
      -g: Directory containing genome database (No Default) ## Remove this and get from ETmapper out
      -o: Output dir (Default: ET_stats)
      -p: Make plots (Default: FALSE)
      -r: Remove specific genome from analysis (must be quoted)
      
      Basic Stats Summary Options:
      
      -q: MapQ cutoff score (Default: 20)
      -m: Maximum read mismatches allowed (Default: 5; SE only)
      -C: Custom filter function (Overrides PE filters; must be quoted)
      -E: Turn on basic transformation efficency calculation (Default: OFF)
      -S: Spike in control organism for efficency calculation (Default: Bacteroides_thetaiotaomicron_VPI_5482_GCF_000011065.1)
      -h: Bring up this help menu
  
  
      Example Usage:
      
      ## Stats summary
      ETstats.R -w ss -d ./ETmapper_out -b sample_info.txt
      ## Stats summary custom filter
      ETstats.R -w ss -d ./ETmapper_out -b sample_info.txt -c 'GENOME1 == GENOME2 & MAPQ1 > 5'
      
    ")
  
  q(save="no")
}


## Mandatory Arguments

# Work Flow Type
wf <- args[which(args == "-w") + 1]
if(wf %notin% c("ss","rs","ip","cs")){
  cat(paste0("\n",wf," is not a known workflow...exiting"))
  q(save="no")
}

# ETmapper or ETstats output directory
dat_in <- args[which(args == "-d") + 1]



## Optional Arguments

# Sample Info
if("-s" %in% args){
si <- args[which(args == "-s") + 1]
sample_info <- read.table(si, sep = "\t", header = F)
}

# Genome Database Directory
#gd <- args[which(args == "-g") + 1]


# Output directory (Will be created if not specified)
out_dir <- "ET_stats"
if("-o" %in% args){
  out_dir <- args[which(args == "-o") + 1]
}
dir.create(out_dir, recursive = T)


# Make plots (Default: FALSE)
plot_res <- FALSE
if("-p" %in% args){
  plot_res <- TRUE
}

# Remove specific genome from analysis
if("-r" %in% args){
  rm_gen <- args[which(args == "-r") + 1]
  ### STRING SPLIT
}



## Basic stats summary options

# Mapq Cutoff (Default: 20)
mq_cut <- 20
if("-q" %in% args){
  mq_cut <- as.numeric(args[which(args == "-q") + 1])
}

# Mapq Cutoff (Default: 5)
mm_cut <- 5
if("-m" %in% args){
  mm_cut <- as.numeric(args[which(args == "-m") + 1])
}

# Custom filter function (Overrides PE filters)
if("-C" %in% args){
  custom_filt <- args[which(args == "-C") + 1]
}

# Calculate transformation efficency (Default: FALSE)
calc_eff <- FALSE
if("-E" %in% args){
  calc_eff <- TRUE
}

# Spike in control organism for efficency calculation (Default: Bacteroides_thetaiotaomicron_VPI_5482_GCF_000011065.1)
spike_in_org <- "Bacteroides_thetaiotaomicron_VPI_5482_GCF_000011065.1"
if("-S" %in% args){
  spike_in_org <- args[which(args == "-S") + 1]
}


# Load Testing Data
# etstats_version <- c("v0.03")
# wf <- "ss"
# dat_in <- "~/Desktop/ET_test/ET_mapper/"
# out_dir <- "~/Desktop/ET_test/ET_stats"
# mq_cut <- 20
# mm_cut <- 5
# plot_res <- FALSE
# calc_eff <- FALSE


# sample_info <- read.table(si, sep = "\t", header = F)
# jm_pe_example <- fread("test_data/ETmapper/paired_end_test/hits/JD_ZZ_2ndcycle_1.hits")
# lm_pe_example <- fread("~/Desktop/Meta_2020_01_14_1.mghits")
# se_example <- fread("test_data/ETmapper/single_end_test/hits/JD_ZZ_2ndcycle_1.hits")
# 
# tmp_hit_table <- se_example
# i = 1






#### BEGIN FUNCTION DEFINE ####


### Setup Functions

## Function 1: Make a basic log file that will be appended based on analysis
make.log.file <- function(){
  
  cat(paste0("ETstats ",etstats_version,"   Beginning Analysis...\n\n"))

  cat(
    paste0("ETstats ",etstats_version," Summary    Created: ", date()),"\n\n",
    "Program Parameters:\n",
    paste0("Workflow type is: ", wf),"\n",
    file = paste0(out_dir,"/run_log.txt"))
}

## Function 2: Identify what analysis has been run by ETmapper or ETstats already, store, and log
## Return: env_summary
check.env <- function(){
  
  # Env check for ss workflow
  if (wf == "ss"){
  
    # Make data.frame to hold environment variables
    env_summary <- data.frame(jm = FALSE,
                              jm_dir = "A1",
                              jm_gd = "A2",
                              lm = FALSE,
                              lm_dir = "B1",
                              lm_gd = "B2",
                              stringsAsFactors = F)
  
    # 1 check ETmapper out for jm folder, store location, and output to log
    jm_exist <- dir.exists(paste0(dat_in,"jm"))
    
    if (jm_exist == TRUE){
      cat("Junction Mapping Directory Found...\n")
      env_summary[1,1] <- TRUE
      env_summary[1,2] <- normalizePath(paste0(dat_in,"jm"))
      env_summary[1,3] <- system(paste0("grep 'Genome Database' ",dat_in,"jm/run_log.txt | awk '{print $3}'"), intern = T)
      cat(
        paste0("Found ETmapper jm Output: TRUE\n"),
        paste0("ETmapper jm Output Dir: ",env_summary[1,2],"\n"),
        paste0("ETmapper jm Genome Dir: ",env_summary[1,3],"\n"),
        file = paste0(out_dir,"/run_log.txt"), append = T)
    }
    if (jm_exist == FALSE){
      cat("Junction Mapping Directory Not Found...\n")
      env_summary[1,1] <- FALSE
      env_summary[1,2] <- NA
      env_summary[1,3] <- NA
      cat(
        paste0("Found ETmapper jm Output: FALSE\n"),
        paste0("ETmapper jm Output Dir: NA\n"),
        paste0("ETmapper jm Genome Dir: NA\n"),
        file = paste0(out_dir,"/run_log.txt"), append = T)
    }
    
    # 2 check ETmaper out for lm folder, store location, and output to log 
    lm_exist <- dir.exists(paste0(dat_in,"lm"))
    
    if (lm_exist == TRUE){
      cat("Lite Metagenomics Directory Found...\n")
      env_summary[1,4] <- TRUE
      env_summary[1,5] <- normalizePath(paste0(dat_in,"lm"))
      env_summary[1,6] <- system(paste0("grep 'Genome Database' ",dat_in,"lm/run_log.txt | awk '{print $3}'"), intern = T)
      cat(
        paste0("Found ETmapper lm Output: TRUE\n"),
        paste0("ETmapper lm Output Dir: ",env_summary[1,5],"\n"),
        paste0("ETmapper lm Genome Dir: ",env_summary[1,6],"\n"),
        file = paste0(out_dir,"/run_log.txt"), append = T)
    }
    if (lm_exist == FALSE){
      cat("Lite Metagenomics Directory Not Found...\n")
      env_summary[1,4] <- FALSE
      env_summary[1,5] <- NA
      env_summary[1,6] <- NA
      cat(
        paste0("Found ETmapper lm Output: FALSE\n"),
        paste0("ETmapper lm Output Dir: NA\n"),
        paste0("ETmapper lm Genome Dir: NA\n"),
        file = paste0(out_dir,"/run_log.txt"), append = T)
    }
    
    # If neither jm or lm are found throw error and quit
    if ((env_summary[1,1] == FALSE & env_summary[1,4] == FALSE) == TRUE){
      cat("Could not find output directories for jm or lm workflows in ETmapper output.\n
          Program Quitting...")
      q(save = "no")
    }
    
    # Output the env_summary data.frame
    return(env_summary)
  }
  
  
  # Env check for cs workflow
  if (wf == "cs"){
    
  }
  
}

## Function 3: Get list of files for analysis
## Return: input_data_names
get.files <- function(){
  
  # Get hit files for ss workflow
  if (wf =="ss"){
    
    # Get jm workflow hit file names
    if(env_summary$jm == TRUE){
      jm_file_list <- list.files(paste0(env_summary$jm_dir,"/hits/"))
    }else{jm_file_list <- NA}
    
    # Get lm workflow hit file names
    if(env_summary$lm == TRUE){
      lm_file_list <- list.files(paste0(env_summary$lm_dir,"/hits/"))
    }else{lm_file_list <- NA} 
    
    # Create list of file names
    input_data_names <- list()
    input_data_names$jm_file_list <- jm_file_list
    input_data_names$lm_file_list <- lm_file_list
    return(input_data_names)
  }
  
  if (wf == "cs"){
    
  }
  
}

## Function 4: Check if data processed as paired end or single end
## Return: pe_test_out
pe.test <- function(){
  
  if (wf == "ss"){
    
    # Generate df to hold pe.test output; logical jm =T/F,lm=T/F
    pe_test_out <- data.frame(jm = NA, lm = NA)
  
    # Check if jm data is present, identify if PE, update pe_test_out and log file
    if (env_summary[1,1] == TRUE){
      test_hits_file <- fread(paste0(env_summary[1,2],"/hits/",input_data_names$jm_file_list[1]), nrows = 5)
      
      if (sum(c("GENOME1", "GENOME2") %in% colnames(test_hits_file)) == 2){
        pe_test_out[1,1] <- TRUE
        cat("Workflow jm Data Paired End: TRUE\n", file = paste0(out_dir,"/run_log.txt"), append = T)
      }
      
      if (sum(c("GENOME1", "GENOME2") %in% colnames(test_hits_file)) == 1){
        pe_test_out[1,1] <- FALSE
        cat("Workflow jm Data Paired End: FALSE\n", file = paste0(out_dir,"/run_log.txt"), append = T)
      }
    }
    
    # Check if lm data is present, identify if PE, update pe_test_out and log file
    if (env_summary[1,4] == TRUE){
      test_hits_file <- fread(paste0(env_summary[1,5],"/hits/",input_data_names$lm_file_list[1]), nrows = 5)
      
      if (sum(c("GENOME1", "GENOME2") %in% colnames(test_hits_file)) == 2){
        pe_test_out[1,2] <- TRUE
        cat("Workflow lm Data Paired End: TRUE\n", file = paste0(out_dir,"/run_log.txt"), append = T)
      }
      
      if (sum(c("GENOME1", "GENOME2") %in% colnames(test_hits_file)) == 1){
        pe_test_out[1,2] <- FALSE
        cat("Workflow lm Data Paired End: FALSE\n", file = paste0(out_dir,"/run_log.txt"), append = T)
      }
    }
    
    # Return pe_test_out data.frame
    return(pe_test_out)
  }
  
  
  ## PE TEST FOR FUTURE WORKFLOW
  if (wf == "ss"){
  }
  
}  


### Analysis Functions

## Function 5: Summarize hits per genomes across samples for paired end data
jm.summary.pe <- function(){
  
  # Make df to hold final hit table
  all_hits <- data.table()
  filt_hits <- data.table()
  sample_names <- sub(".hits","",input_data_names$jm_file_list)
  
  # Build master df of all hits
  for (i in 1:length(sample_names)){
    
    # Say what is happening
    cat(paste0("Calculating Summary Stats for jm Sample: ",sample_names[i],"\n"))
    
    # Read in hit table
    tmp_hit_table <- fread(paste0(env_summary$jm_dir,"/hits/",sample_names[i],".hits"), sep = "\t")
    
    # Filter out all NA in GENOME1 column so that no reads mapping to GENOME 1 and 2 == NA get through
    tmp_hit_table <- tmp_hit_table[!(is.na(GENOME1))] ### THIS IS A SHITTY FIX 
    
    
    ### APPLY FILTERING (each filtering step is implemented individually for now so we can add options)
    
    # Create table to be filtered on
    tmp_filt_table <- tmp_hit_table
    
    # Custom filters if wated
    if(exists("custom_filt") == TRUE){
      
      # Hits to same genome
      tmp_filt_table <- subset(tmp_filt_table, eval(parse(text=custom_filt)))
      
    # Basic filters if no custom specified
    } else {
      
      # Hits to same genome
      tmp_filt_table <- subset(tmp_filt_table, GENOME1 == GENOME2)
      
      # Filter by mapq
      tmp_filt_table <- subset(tmp_filt_table, MAPQ1 >= mq_cut | MAPQ2 >= mq_cut)
      
    }

    # Make filtered data output folder and write tables out
    dir.create(paste0(out_dir,"/hits_filt/jm"), recursive = T)
    write.table(tmp_filt_table, file = paste0(out_dir,"/hits_filt/jm/",sample_names[i],".hits.filt"), quote = F, row.names = F, sep = "\t")
    
    # Sumarize raw data
    tmp_raw_summary <- data.table(SAMPLE = sample_names[i],
                                  tmp_hit_table %>% 
                                    group_by(GENOME1) %>% ## This may need to be fixed for pe if genomes are not equal (decide which read better for GENOME call)
                                    summarise(READ_RAW = n(),
                                              UNIQ_RAW = n_distinct(barcodes),
                                              BPR_RAW = n_distinct(barcodes) / n(),
                                              MMQ1_RAW = mean(MAPQ1,na.rm = T),
                                              MMM1_RAW = mean(NM1,na.rm = T),
                                              MLEN1_RAW = mean(LEN1,na.rm = T),
                                              MQUAL1_RAW = mean(QUAL1,na.rm = T),
                                              MQUALNM1_RAW = mean(QUALNM1,na.rm = T),
                                              MMQ2_RAW = mean(MAPQ2,na.rm = T),
                                              MMM2_RAW = mean(NM2,na.rm = T),
                                              MLEN2_RAW = mean(LEN2,na.rm = T),
                                              MQUAL2_RAW = mean(QUAL2,na.rm = T),
                                              MQUALNM2_RAW = mean(QUALNM2,na.rm = T)))

    
    # Sumarize filtered data
    tmp_filt_summary <- data.table(SAMPLE = sample_names[i],
                                   tmp_filt_table %>% 
                                     group_by(GENOME1) %>% ## This may need to be fixed for pe if genomes are not equal (decide which read better for GENOME call)
                                     summarise(READ_FLT = n(),
                                               UNIQ_FLT = n_distinct(barcodes),
                                               BPR_FLT = n_distinct(barcodes) / n(),
                                               MMQ1_FLT = mean(MAPQ1),
                                               MMM1_FLT = mean(NM1),
                                               MLEN1_FLT = mean(LEN1),
                                               MQUAL1_FLT = mean(QUAL1),
                                               MQUALNM1_FLT = mean(QUALNM1),
                                               MMQ2_FLT = mean(MAPQ2),
                                               MMM2_FLT = mean(NM2),
                                               MLEN2_FLT = mean(LEN2),
                                               MQUAL2_FLT = mean(QUAL2),
                                               MQUALNM2_FLT = mean(QUALNM2)))
    
    # Create table of all hits for all samples
    all_hits <- rbind(all_hits, tmp_raw_summary)
    
    # Create table of filt hits for all samples
    filt_hits <- rbind(filt_hits, tmp_filt_summary)
    
  }
  
  # Merge all and filtered hits and create output
  pe_hits_out <- merge(all_hits, filt_hits, by = c("SAMPLE", "GENOME1"), all.x = T)
  pe_hits_out_full <- pe_hits_out[order(SAMPLE, -UNIQ_FLT)]
  pe_hits_out_summary <- pe_hits_out_full[,c(1:4,16:18)]
  pe_hits_out_list <- list(pe_hits_out_full, pe_hits_out_summary)
  return(pe_hits_out_list)
  
}

## Function 6: Summarize hits per genomes across samples for single end data
jm.summary.se <- function(){
  
  # Make df to hold final hit table
  all_hits <- data.table()
  filt_hits <- data.table()
  
  # Build master df of all hits
  for (i in 1:nrow(sample_names)){
    
    # Read in hit table
    tmp_hit_table <- fread(paste0(dat_in,"hits/",sample_names[i,1],".hits"), sep = "\t")
    
    
    ### APPLY FILTERING 
    tmp_filt_table <- subset(tmp_hit_table, MAPQ >= mq_cut & NM <= mm_cut)
    
  
    # Sumarize raw data
    tmp_raw_summary <- data.table(SAMPLE = sample_names[i,1],
                                  GROUP = sample_names[i,2],
                                  tmp_hit_table %>% 
                                  group_by(GENOME) %>% 
                                  summarise(READ_RAW = n(),
                                            UNIQ_RAW = n_distinct(barcodes),
                                            BPR_RAW = n_distinct(barcodes) / n(),
                                            MMQ_RAW = mean(MAPQ),
                                            MMM_RAW = mean(NM),
                                            MLEN_RAW = mean(LEN)))
 
    # Sumarize filtered data
    tmp_filt_summary <- data.table(SAMPLE = sample_names[i,1],
                                   GROUP = sample_names[i,2],
                                   tmp_filt_table %>% 
                                   group_by(GENOME) %>% 
                                   summarise(READ_FLT = n(),
                                             UNIQ_FLT = n_distinct(barcodes),
                                             BPR_FLT = n_distinct(barcodes) / n(),
                                             MMQ_FLT = mean(MAPQ),
                                             MMM_FLT = mean(NM),
                                             MLEN_FLT = mean(LEN)))
    
    # Create table of all hits for all samples
    all_hits <- rbind(all_hits, tmp_raw_summary)
    
    # Create table of filt hits for all samples
    filt_hits <- rbind(filt_hits, tmp_filt_summary)
    
  }
  
  # Merge all and filtered hits and create output
  se_hits_out <- merge(all_hits, filt_hits, by = c("SAMPLE", "GROUP", "GENOME"), all.x = T)
  se_hits_out <- se_hits_out[order(SAMPLE, -UNIQ_FLT)]
  se_hits_out
  
}
  
## Function 7: Summarize Coverage per genomes across samples for paired end data
lm.summary.pe <- function(){
  
  # Load additional required data
  genome_stats <- fread(paste0(env_summary$lm_gd,"genome_stats.txt"), sep ="\t")
  
  # Make df to hold final hit table
  all_hits <- data.table()
  filt_hits <- data.table()
  sample_names <- sub(".mghits","",input_data_names$lm_file_list)
  
  
  # Build master df of all hits
  for (i in 1:length(sample_names)){
    
    # Say what is happening
    cat(paste0("Calculating Summary Stats for lm Sample: ",sample_names[i],"\n"))
    
    # Read in hit table
    tmp_hit_table <- fread(paste0(env_summary$lm_dir,"/hits/",sample_names[i],".mghits"), sep = "\t")
    
    
    ### APPLY FILTERING (each filtering step is implemented individually for now so we can add options)
    
    # Filter out all NA in GENOME1 column so that no reads mapping to GENOME 1 and 2 == NA get through
    tmp_hit_table <- tmp_hit_table[!(is.na(GENOME1))] ### THIS IS A SHITTY FIX 
    
    # Create table to be filtered on
    tmp_filt_table <- tmp_hit_table
    
    # Custom filters if wated
    if(exists("custom_filt") == TRUE){
      
      # Hits to same genome ##### THIS NEEDS TO BE ADDRESSED 
      tmp_filt_table <- subset(tmp_filt_table, eval(parse(text=custom_filt)))
      
      # Basic filters if no custom specified
    } else {
      
      # Apply All Quality Filters
      tmp_filt_table <- subset(tmp_filt_table,
                               (GENOME1 == GENOME2) &
                               (MAPQ1 >= mq_cut | MAPQ2 >= mq_cut) &
                               (STRAND1 != STRAND2))  
    }
    
    # Sumarize raw data
    tmp_raw_summary <- data.table(SAMPLE = sample_names[i],
                                  tmp_hit_table %>% 
                                    group_by(GENOME1) %>% ## This may need to be fixed for pe if genomes are not equal (decide which read better for GENOME call)
                                    summarise(READ_RAW = n(),
                                              TOTBP_RAW = sum(TOTAL_BP,na.rm = T),
                                              MMQ1_RAW = mean(MAPQ1,na.rm = T),
                                              MMM1_RAW = mean(NM1,na.rm = T),
                                              MLEN1_RAW = mean(LEN1,na.rm = T),
                                              MQUAL1_RAW = mean(QUAL1,na.rm = T),
                                              MQUALNM1_RAW = mean(QUALNM1,na.rm = T),
                                              MMQ2_RAW = mean(MAPQ2,na.rm = T),
                                              MMM2_RAW = mean(NM2,na.rm = T),
                                              MLEN2_RAW = mean(LEN2,na.rm = T),
                                              MQUAL2_RAW = mean(QUAL2,na.rm = T),
                                              MQUALNM2_RAW = mean(QUALNM2,na.rm = T)))
    
    
    # Sumarize filtered data
    tmp_filt_summary <- data.table(SAMPLE = sample_names[i],
                                   tmp_filt_table %>% 
                                     group_by(GENOME1) %>% ## This may need to be fixed for pe if genomes are not equal (decide which read better for GENOME call)
                                     summarise(READ_FLT = n(),
                                               TOTBP_FLT = sum(TOTAL_BP,na.rm = T),
                                               MMQ1_FLT = mean(MAPQ1),
                                               MMM1_FLT = mean(NM1),
                                               MLEN1_FLT = mean(LEN1),
                                               MQUAL1_FLT = mean(QUAL1),
                                               MQUALNM1_FLT = mean(QUALNM1),
                                               MMQ2_FLT = mean(MAPQ2),
                                               MMM2_FLT = mean(NM2),
                                               MLEN2_FLT = mean(LEN2),
                                               MQUAL2_FLT = mean(QUAL2),
                                               MQUALNM2_FLT = mean(QUALNM2)))
    
    
    # Merge in genome data and calculate coverage stats
    tmp_filt_summary <- merge(tmp_filt_summary, genome_stats[,c(1,3)], by.x = "GENOME1", by.y = "Genome", all.x = T)
    
    tmp_filt_summary <- data.table(tmp_filt_summary,
                                   RPK_FLT = tmp_filt_summary$READ_FLT / (tmp_filt_summary$Size / 1000),
                                   RPK_FLT_FRAC = (tmp_filt_summary$READ_FLT / (tmp_filt_summary$Size / 1000)) / sum(tmp_filt_summary$READ_FLT / (tmp_filt_summary$Size / 1000)),
                                   COV_FLT = tmp_filt_summary$TOTBP_FLT / tmp_filt_summary$Size,
                                   COV_FLT_FRAC = (tmp_filt_summary$TOTBP_FLT / tmp_filt_summary$Size) / sum(tmp_filt_summary$TOTBP_FLT / tmp_filt_summary$Size))
    
    # Create table of all hits for all samples
    all_hits <- rbind(all_hits, tmp_raw_summary)
    
    # Create table of filt hits for all samples
    filt_hits <- rbind(filt_hits, tmp_filt_summary)
    
  
  } ## End for loop for calculating each genome summary stats
  
  # Merge all and filtered hits and create output
  pe_hits_out <- merge(all_hits, filt_hits, by = c("SAMPLE", "GENOME1"), all.x = T)
  pe_hits_out_full <- pe_hits_out[order(SAMPLE, -TOTBP_FLT)]
  pe_hits_out_summary <- pe_hits_out_full[,c(1:4,15,16,28:31)]
  pe_hits_out_list <- list(pe_hits_out_full, pe_hits_out_summary)
  return(pe_hits_out_list)


  
}

## Function 8: Summarize Coverage per genomes across samples for single end data
lm.summary.se <- function(){
  
}

## Function 9: Calculate transformation efficency from combined jm/lm output 
## RETURN: eff_out
calc.trans.eff <- function(){
  
  # Say what is happening
  cat("\nCalculating Transformation Efficency for Samples...\n")
  
  # Create data frame to hold output
  eff_out <- data.table()
  
  # Loop over samples and pull individual data
  for (i in 1:nrow(sample_info)){
    
    # Get names for first set of corresponding samples
    tmp_jm_name <- as.character(sample_info[i,1])
    tmp_lm_name <- as.character(sample_info[i,2])
    
    # Pull data on only corresponding samples 
    tmp_jm_dat <- subset(ETstats_jm_out[[1]], SAMPLE == tmp_jm_name)
    tmp_lm_dat <- subset(ETstats_lm_out[[1]], SAMPLE == tmp_lm_name)
    
    # Merge lm data onto jm data 
    tmp_comb_dat <- merge(tmp_jm_dat, tmp_lm_dat, by.x = "GENOME1", by.y = "GENOME1", all.x = T)
    
    # Get Spike in Org Data and remove from tmp_comb_dat
    sio_tmp_dat <- subset(tmp_comb_dat, GENOME1 == spike_in_org)
    tmp_comb_dat <- tmp_comb_dat[-which(tmp_comb_dat$GENOME1 == spike_in_org),]
    
    # Calculate Efficency statistics for each genome
    sio_uniq_flt <- as.numeric(sio_tmp_dat$UNIQ_FLT)
    sio_rpk_frac <- as.numeric(sio_tmp_dat$RPK_FLT_FRAC)
    sio_cov_frac <- as.numeric(sio_tmp_dat$COV_FLT_FRAC)
    tmp_eff_dat <- data.table(SAMPLE = tmp_comb_dat$SAMPLE.x,
                              GENOME1 = tmp_comb_dat$GENOME1,
                              SIZE = tmp_comb_dat$Size,
                              RPK_FRAC = tmp_comb_dat$RPK_FLT_FRAC * 100,
                              COV_FRAC = tmp_comb_dat$COV_FLT_FRAC * 100,
                              SIO_UNIQ = sio_uniq_flt,
                              SIO_RPK_FRAC = sio_rpk_frac,
                              SIO_COV_FRAC = sio_cov_frac,
                              READS = tmp_comb_dat$READ_FLT.x,
                              UNIQ = tmp_comb_dat$UNIQ_FLT,
                              BPR = tmp_comb_dat$BPR_FLT,
                              RPK_EFF = ((tmp_comb_dat$UNIQ_FLT/tmp_comb_dat$RPK_FLT_FRAC) * 100)/(sio_uniq_flt/sio_rpk_frac), 
                              COV_EFF = ((tmp_comb_dat$UNIQ_FLT/tmp_comb_dat$COV_FLT_FRAC) * 100)/(sio_uniq_flt/sio_cov_frac))
    
    # Repalce NA in data with zeros
    tmp_eff_dat[is.na(tmp_eff_dat)] <- 0 
    
    # Add calculations to output data frame
    eff_out <- rbind(eff_out, tmp_eff_dat)
  
  }
  
  # Return data
  return(eff_out)
  
}

## Function 4: Basic Plotting of genomes (genomes v sample groups)
  

## Function 5: Read summary stats ### SPLIT THIS FUNCTION EVENTUALLY FOR READABLITY
read.stats <- function(){
  
  # Make df and list to hold summary and full hit tables
  #aggregate_hits <- data.table() --- Don't do this right now
  read_summary <- data.table()
  
  # loop over batch file depending on data type
  if (paired_end_data == TRUE){
    
    # Build master df of all hits
    for (i in 1:nrow(sample_names)){
      
      # Read in hit table
      tmp_hit_table <- fread(paste0(dat_in,"hits/",sample_names[i,1],".hits"), sep = "\t")
      
      # Remove genomes if in arguments
      if (exists("rm_gen") == TRUE){
        tmp_hit_table <- subset(tmp_hit_table, GENOME1 %notin% rm_gen[[1]] & GENOME2 %notin% rm_gen[[1]])
        
      }
        
      # Sumarize raw data on reads
      tmp_raw_summary <- data.table(SAMPLE = sample_names[i,1],
                                    GROUP = sample_names[i,2],
                                    tmp_hit_table %>% 
                                      summarise(READ_RAW = n(),
                                                UNIQ_RAW = n_distinct(barcodes),
                                                BPR_RAW = n_distinct(barcodes) / n(),
                                                MMQ1_RAW = mean(MAPQ1,na.rm = T),
                                                MMM1_RAW = mean(NM1,na.rm = T),
                                                MLEN1_RAW = mean(LEN1,na.rm = T),
                                                MQUAL1_RAW = mean(QUAL1,na.rm = T), # QUAL = number of bases above phred20
                                                MQUALNM1_RAW = mean(QUALNM1,na.rm = T), #QUALNM = number of mismatches ONLY at bases above phred20
                                                MANI1_RAW = mean(NM1 / LEN1, na.rm = T), # raw ANI = number of differences including errors
                                                MHQANI1_RAW = mean(QUALNM1 / QUAL1, na.rm = T), # high qual ANI = number of differences mostly not including errors
                                                MPQ1_RAW = mean(QUAL1 / LEN1, na.rm = T), # percent quality = the percentage of high quality bases on the read
                                                MMQ2_RAW = mean(MAPQ2,na.rm = T),
                                                MMM2_RAW = mean(NM2,na.rm = T),
                                                MLEN2_RAW = mean(LEN2,na.rm = T),
                                                MQUAL2_RAW = mean(QUAL2,na.rm = T), # QUAL = number of bases above phred20
                                                MQUALNM2_RAW = mean(QUALNM2,na.rm = T), # QUALNM = number of mismatches ONLY at bases above phred20
                                                MANI2_RAW = mean(NM2 / LEN2, na.rm = T), # raw ANI = number of differences including errors
                                                MHQANI2_RAW = mean(QUALNM2 / QUAL2, na.rm = T), # high qual ANI = number of differences mostly not including errors
                                                MPQ2_RAW = mean(QUAL2 / LEN2, na.rm = T))) # percent quality = the percentage of high quality bases on the read
      
      # Create taable of read summary statistics for each sample
      read_summary <- rbind(read_summary, tmp_raw_summary)
      
    } ## End looping over PE files
     
    
    # Plot Summary Data if plotting == T in args
    if (plot_res == TRUE){
      
      # Melt data for plotting
      read_summary_plot <- melt(read_summary)
      
      # Produce ggplot
      ggplot(read_summary_plot, aes(x = GROUP, y = value, color = GROUP)) +
        geom_point() +
        #geom_boxplot() +
        facet_wrap(. ~ variable, scales = "free_y")
      
      ggsave("read_stats.pdf", device = "pdf", width = 20, height = 20)
       
    } 
    
  } ## End paired end summmary analysis
  
  
  
  if (paired_end_data == FALSE){
   
    # Build master df of all hits
    for (i in 1:nrow(sample_names)){
      
      # Read in hit table
      tmp_hit_table <- fread(paste0(dat_in,"hits/",sample_names[i,1],".hits"), sep = "\t")
      
      # Remove genomes if in arguments
      if (exists("rm_gen") == TRUE){
        tmp_hit_table <- subset(tmp_hit_table, GENOME1 %notin% rm_gen[[1]])
        
      }
      
      # Sumarize raw data
      tmp_raw_summary <- data.table(SAMPLE = sample_names[i,1],
                                    GROUP = sample_names[i,2],
                                    tmp_hit_table %>% 
                                      summarise(READ_RAW = n(),
                                                UNIQ_RAW = n_distinct(barcodes),
                                                BPR_RAW = n_distinct(barcodes) / n(),
                                                MMQ_RAW = mean(MAPQ),
                                                MMM_RAW = mean(NM),
                                                MLEN_RAW = mean(LEN)))
      
      # Create taable of read summary statistics for each sample
      read_summary <- rbind(read_summary, tmp_raw_summary)
      
    } ## End looping over PE files
    
    
    # Plot Summary Data if plotting == T in args
    if (plot_res == TRUE){
      
      # Melt data for plotting
      read_summary_plot <- melt(read_summary)
      
      # Produce ggplot
      ggplot(read_summary_plot, aes(x = GROUP, y = value, color = GROUP)) +
        geom_point() +
        #geom_boxplot() +
        facet_wrap(. ~ variable, scales = "free_y")
      
      ggsave("read_stats.pdf", device = "pdf", width = 20, height = 20)
      
    }
    
  } ## End SE summmary analysis

  
  # Return Summary Data 
  read_summary
  
} ## End read.stats function


## Function 6: Comparative Statistics between groups (hits or read summaries)



#### END FUNCTION DEFINE ####




#### BEGIN WORKFLOWS ####


### Basic Stats Summary Workflow
if (wf == "ss"){
  
  ## Create Correct Dir Structure for Workflow
  out_dir <- paste0(out_dir,"/ss")
  dir.create(out_dir, recursive = T)
  
  ## Initalize program and create log file 
  make.log.file()
  
  ## Determine environment and important directories
  env_summary <- check.env()

  ## Create data.frame of input data file names 
  input_data_names <- get.files()
  
  ## Check if data is pe or se for all workflows run
  pe_test_out <- pe.test()
  
  ### Perform and Output Statistical Summaries for jm Workflow if present
  if(env_summary$jm == TRUE){
    
    if(pe_test_out$jm == TRUE){
      ETstats_jm_out <- jm.summary.pe()
    }
    
    if (pe_test_out$jm == FALSE){
      ETstats_jm_out <- jm.summary.se()
    }
    
    # Write Output
    write.table(ETstats_jm_out[[1]], file = paste0(out_dir,"/ETstats_jm_full.txt"), quote = F, row.names = F, sep = "\t")
    write.table(ETstats_jm_out[[2]], file = paste0(out_dir,"/ETstats_jm_summary.txt"), quote = F, row.names = F, sep = "\t")
    
  }
  
  ### Perform and Output Statistical Summaries for lm Workflow if present
  if(env_summary$lm == TRUE){
    
    if(pe_test_out$lm == TRUE){
      ETstats_lm_out <- lm.summary.pe()
    }
    
    if (pe_test_out$lm == FALSE){
      ETstats_lm_out <- lm.summary.se()
    }
    
    # Write Output
    write.table(ETstats_lm_out[[1]], file = paste0(out_dir,"/ETstats_lm_full.txt"), quote = F, row.names = F, sep = "\t")
    write.table(ETstats_lm_out[[2]], file = paste0(out_dir,"/ETstats_lm_summary.txt"), quote = F, row.names = F, sep = "\t")
    
  }
  
  ### Perform Efficency Calculation if requested and sample_info present
  if(calc_eff == TRUE & exists("sample_info") == TRUE){
  
      eff_out <- calc.trans.eff()
      write.table(eff_out, file = paste0(out_dir,"/ETstats_eff_summary.txt"), quote = F, row.names = F, sep = "\t")
  }
  
  cat("\nss Workflow completed successfully :-)\n")
  
}  
  
    
 



# if (wf == "rs"){
#  
#   # Check if data is pe or se
#   paired_end_data <- pe.test()
#   
#   # loop over batch file depending on data type
#   read_stats_out <- read.stats()
#   
#   # Write Output
#   write.table(read_stats_out, file = out_dir, quote = F, row.names = F, sep = "\t")
#  
#    
# }

# if (wf == "ip"){
#   
#   cat("ip does not yet exist")
# }







#!/usr/bin/env Rscript

### Set Version
etstats_version <- c("v0.04")

### Check and Load Libraries

# Check data.table
if ("data.table" %in% installed.packages() == FALSE){
  print("Please install R package data.table. Program quitting...")
  q(save="no")
}

# Check dplyr
if ("dplyr" %in% installed.packages() == FALSE){
  print("Please install R package dplyr. Program quitting...")
  q(save="no")
}

# Check plyr
if ("plyr" %in% installed.packages() == FALSE){
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
      ETstats v0.04
      
      Usage: ETstats.R -w [workflow] -d [Data_dir] [options]

      Mandatory Arguments:
      
      -w: Workflow type (No Default)
          ss - basic stats summary
          ts - targeted insertion summary
          cs - comparative stats summary
          rs - read statistics
          ip - insertion plot
      -d: Directory with data for analysis (No Default; can be ETmapper or ETstats)
      
      Optional Arguments:
      
      -s: Sample info (No Default; This file format will depend on the workflow. See documentation)
      -o: Output dir (Default: ET_stats)
      -p: Make plots (Default: FALSE)
      -r: Remove specific genome from analysis (must be quoted; Non Functional)
      
      Basic Stats Summary Options:
      
      -q: MapQ cutoff score (Default: 20)
      -m: Maximum read mismatches allowed (Default: 5; SE only)
      -C: Custom filter function (Overrides PE filters; must be quoted)
      -R: Multiread threshold for unique barcode detection (Default: 2)
      -E: Turn on basic transformation efficency calculation (Default: OFF)
      -S: Spike in control organism for efficency calculation (Default: Bacteroides_thetaiotaomicron_VPI_5482_GCF_000011065.1)
      -D: Conjugal Donor organism (Default: Escherichia_coli_BW25113_GCF_000750555.1)
      -V: Vector sequences included in mapping database (Default: 'pHLL249 pHLL250'; must be quoted)
      -h: Bring up this help menu
      
      Targeted Insertion Summary Options:
      
      -W: Positive target window (Defalut: 200bp)
      -T: Targeted transposon models (Default: 'VcCas_model'; must be quoted)
  
      
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
}

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

# Number of good reads to count a barcode in a genome (Default: 2)
bc_rep <- 2
if("-R" %in% args){
  bc_rep <- as.numeric(args[which(args == "-R") + 1])
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

# Conjugal Donor organism (Default: Escherichia_coli_BW25113_GCF_000750555.1)
donor_org <- "Escherichia_coli_BW25113_GCF_000750555.1"
if("-D" %in% args){
  donor_org <- args[which(args == "-D") + 1]
}

# Vector sequences included in mapping database (Default: "pHLL249 pHLL250" )
vector_names <- c("pHLL249", "pHLL250")
if("-V" %in% args){
  vector_names <- args[which(args == "-V") + 1]
  vector_names <- unlist(strsplit(vector_names, split = " "))
}



## Basic stats summary options

# Positive target window (Defalut: 200bp)
window <- 200
if("-W" %in% args){
  window <- as.numeric(args[which(args == "-W") + 1])
}

# Targeted transposon models (Default: "VcCas_model")
Tn_names <- c("VcCas_model")
if("-T" %in% args){
  vector_names <- args[which(args == "-T") + 1]
  vector_names <- unlist(strsplit(vector_names, split = " "))
}



##### TESTING DATA #####

### TEST ALL
# '%notin%' <- Negate('%in%')
# etstats_version <- c("v0.04")
# plot_res <- FALSE
# calc_eff <- FALSE


### TEST SS - PE WORKFLOW
# wf <- "ss"
# spike_in_org <- "Bacteroides_thetaiotaomicron_VPI_5482_GCF_000011065.1"
# donor_org <- "Escherichia_coli_BW25113_GCF_000750555.1"
# vector_names <- c("pHLL249", "pHLL250")
# mq_cut <- 20
# mm_cut <- 5
# bc_rep <- 2
# dat_in <- "~/Dropbox/Banfield_Lab_Files/Projects/mCAFE/mCAFE_Project/Studies/20_2_27_BR_30_ET_KlebCurve_NT/ET_mapper_dual/"
# out_dir <- "~/Desktop/ET_test/ET_stats"
# tmp_hit_table <- fread("~/Dropbox/Banfield_Lab_Files/Projects/mCAFE/mCAFE_Project/Studies/20_4_6_Kleb_Std_Curve_Final/ET_mapper/jm/hits/BR_1_1_Kcurve_1.hits")


### TEST SS - SE WORKFLOW
# wf <- "ss"
# spike_in_org <- "Bacteroides_thetaiotaomicron_VPI_5482_GCF_000011065.1"
# donor_org <- "Escherichia_coli_BW25113_GCF_000750555.1"
# vector_names <- c("pHLL249", "pHLL250")
# mq_cut <- 20
# mm_cut <- 5
# bc_rep <- 2
# dat_in <- "~/Dropbox/Banfield_Lab_Files/Projects/mCAFE/mCAFE_Project/Studies/20_3_12_BR_9_Barseq/ET_mapper_et0.1/"
# out_dir <- "~/Desktop/ET_test/ET_stats"



### TEST TS - SE WORKFLOW
# wf <- "ts"
# bc_rep <- 2
# window <- 200
# Tn_names <- c("VcCas_model")
# si <- "~/Desktop/ET_test/target_loci.txt"
# sample_info <- fread(si, sep = "\t")
# dat_in <- "~/Desktop/ET_test/ET_stats"
# out_dir <- "~/Desktop/ET_test/ET_stats"
# i <- 1
# j <- 1
# raw_hit_table <- fread("~/Desktop/ET_test/ET_stats/ss/hits_filt/jm/BarSeq_8_T30.hits.filt")
# sample_names <- "BarSeq_8_T30"

##### TESTING DATA #####




#### BEGIN FUNCTION DEFINE ####


### Setup Functions

## Function 1: Make a basic log file that will be appended based on analysis
make.log.file <- function(){
  
  cat(paste0("ETstats ",etstats_version,"   Beginning Analysis...\n\n"))

  cat(
    paste0("ETstats ",etstats_version," Summary    Created: ", date()),"\n\n",
    "Program Parameters:\n", ### Still need to get program params in here 
    paste0("Workflow type is: ", wf),"\n",
    file = paste0(out_dir,"/run_log.txt"))
}

## Function 2: Identify what analysis has been run by ETmapper or ETstats already, store, and log
## Return: env_summary
check.env <- function(){
  
  ## Env check for ss workflow
  if (wf == "ss"){
  
    # Make data.frame to hold environment variables
    env_summary <- data.frame(jm = FALSE,
                              jm_dir = "A1",
                              jm_gd = "A2",
                              jm_pe_map = "A3",
                              lm = FALSE,
                              lm_dir = "B1",
                              lm_gd = "B2",
                              lm_pe_map = "B3",
                              stringsAsFactors = F)
  
    # 1 check ETmapper out for jm folder, store location, and output to log
    jm_exist <- dir.exists(paste0(dat_in,"jm"))
    
    if (jm_exist == TRUE){
      cat("Junction Mapping Directory Found...\n")
      env_summary$jm <- TRUE
      env_summary$jm_dir <- normalizePath(paste0(dat_in,"jm"))
      env_summary$jm_gd <- system(paste0("grep 'Genome Database' ",dat_in,"jm/run_log.txt | awk '{print $3}'"), intern = T)
      env_summary$jm_pe_map <- as.logical(system(paste0("grep 'Paired End Map' ",dat_in,"jm/run_log.txt | awk '{print $4}'"), intern = T))
      cat(
        paste0("Found ETmapper jm Output: TRUE\n"),
        paste0("ETmapper jm Output Dir: ",env_summary$jm_dir,"\n"),
        paste0("ETmapper jm Genome Dir: ",env_summary$jm_gd,"\n"),
        paste0("ETmapper jm mapping PE: ",env_summary$jm_pe_map,"\n"),
        file = paste0(out_dir,"/run_log.txt"), append = T)
    }
    if (jm_exist == FALSE){
      cat("Junction Mapping Directory Not Found...\n")
      env_summary$jm <- FALSE
      env_summary$jm_dir <- NA
      env_summary$jm_gd <- NA
      env_summary$jm_pe_map <- NA
      cat(
        paste0("Found ETmapper jm Output: FALSE\n"),
        paste0("ETmapper jm Output Dir: NA\n"),
        paste0("ETmapper jm Genome Dir: NA\n"),
        paste0("ETmapper jm mapping PE: NA\n"),
        file = paste0(out_dir,"/run_log.txt"), append = T)
    }
    
    # 2 check ETmaper out for lm folder, store location, and output to log 
    lm_exist <- dir.exists(paste0(dat_in,"lm"))
    
    if (lm_exist == TRUE){
      cat("Lite Metagenomics Directory Found...\n")
      env_summary$lm <- TRUE
      env_summary$lm_dir <- normalizePath(paste0(dat_in,"lm"))
      env_summary$lm_gd <- system(paste0("grep 'Genome Database' ",dat_in,"lm/run_log.txt | awk '{print $3}'"), intern = T)
      env_summary$lm_pe_map <- as.logical(system(paste0("grep 'Paired End Map' ",dat_in,"lm/run_log.txt | awk '{print $4}'"), intern = T))
      cat(
        paste0("Found ETmapper lm Output: TRUE\n"),
        paste0("ETmapper lm Output Dir: ",env_summary$lm_dir,"\n"),
        paste0("ETmapper lm Genome Dir: ",env_summary$lm_gd,"\n"),
        paste0("ETmapper lm mapping PE: ",env_summary$lm_pe_map,"\n"),
        file = paste0(out_dir,"/run_log.txt"), append = T)
    }
    if (lm_exist == FALSE){
      cat("Lite Metagenomics Directory Not Found...\n")
      env_summary$lm <- FALSE
      env_summary$lm_dir <- NA
      env_summary$lm_gd <- NA
      env_summary$lm_pe_map <- NA
      cat(
        paste0("Found ETmapper lm Output: FALSE\n"),
        paste0("ETmapper lm Output Dir: NA\n"),
        paste0("ETmapper lm Genome Dir: NA\n"),
        paste0("ETmapper lm mapping PE: NA\n"),
        file = paste0(out_dir,"/run_log.txt"), append = T)
    }
    
    # If neither jm or lm are found throw error and quit
    if ((env_summary$jm == FALSE & env_summary$lm == FALSE) == TRUE){
      cat("Could not find output directories for jm or lm workflows in ETmapper output.\n
          Program Quitting...")
      q(save = "no")
    }
    
    # Output the env_summary data.frame
    return(env_summary)
  }
  
  
  
  ## Env check for ts workflow
  if (wf == "ts"){

    # Make data.frame to hold environment variables
    env_summary <- data.frame(ss = FALSE,
                              ss_dir = "A1",
                              ss_hits = "A2",
                              ss_mode = "A3",
                              stringsAsFactors = F)
    
    # Check ETmapper out for jm folder, store location, and output to log
    ss_exist <- dir.exists(paste0(dat_in,"/ss"))
    
    if (ss_exist == TRUE){
      cat("Summary Stats Dir Found...\n")
      env_summary$ss <- TRUE
      env_summary$ss_dir <- normalizePath(paste0(dat_in,"/ss"))
      env_summary$ss_hits <- normalizePath(paste0(dat_in,"/ss/hits_filt/jm"))
      env_summary$ss_mode <- system(paste0("grep 'Summary Mode (jm)' ",dat_in,"/ss/run_log.txt | awk '{print $4}'"), intern = T)
      cat(
        paste0("ETstats ss Output Dir: ",env_summary$ss_dir,"\n"),
        paste0("ETstats ss Filt Hits Dir: ",env_summary$ss_hits,"\n"),
        paste0("ETstats ss Mode: ",env_summary$ss_mode,"\n"),
        file = paste0(out_dir,"/run_log.txt"), append = T)
    }
    if (ss_exist == FALSE){
      cat("Could not find output directory for ETstats ss.\n
          Program Quitting...")
      q(save = "no")
    }
    
    # Output the env_summary data.frame
    return(env_summary)
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
  
  if (wf == "ts"){
    
    # Get ss workflow hit file names
    ss_file_list <- list.files(paste0(env_summary$ss_dir,"/hits_filt/jm/"), pattern = ".filt")
    return(ss_file_list)  
      
  }
  
}



### Analysis Functions

## Function 1: Summarize hits per genomes across samples for paired end data
jm.summary.pe <- function(){
  
  ### Set up data strucutres
  
  # Make df to hold final hit table
  filt_hits <- data.table()
  sample_names <- sub(".hits","",input_data_names$jm_file_list)
  
  # Create directories to store filtered and bad hits
  dir.create(paste0(out_dir,"/hits_filt/jm"), recursive = T)
  
  # Log this analysis 
  cat(paste0("Summary Mode (jm): PE\n"),
      file = paste0(out_dir,"/run_log.txt"), append = T)
  
  ### Run Filtering and Summarization on each sample and combine
  
  # Build master df of all hits
  for (i in 1:length(sample_names)){
    
    # Say what is happening
    cat(paste0("Calculating Summary Stats for jm Sample: ",sample_names[i],"\n"))
    
    # Read in hit table
    tmp_hit_table <- fread(paste0(env_summary$jm_dir,"/hits/",sample_names[i],".hits"), sep = "\t")
    
    
    ### APPLY FILTERING (each filtering step is implemented individually for now so we can add options)
    
    ## Custom filters if wated
    if(exists("custom_filt") == TRUE){
      
      # Custom filter evaluated from arguments
      tmp_filt_table <- subset(tmp_hit_table, eval(parse(text=custom_filt)))
      
    ## Basic quality filters if no custom specified
    } else {
      
      # Hits to same genome
      tmp_filt_table <- subset(tmp_hit_table, GENOME1 == GENOME2)
      
      # Filter by mapq
      tmp_filt_table <- subset(tmp_filt_table, MAPQ1 >= mq_cut | MAPQ2 >= mq_cut)
      
      # Filter out NA barcodes
      tmp_filt_table <- subset(tmp_filt_table, is.na(barcodes) == FALSE)
      
    }
    
    ## Multiple Barcode Hit Filtering: BC in > 1 genome filtered + written to disk 
    bc2genome <- data.table(tmp_filt_table %>%
                              group_by(barcodes) %>%
                              summarise(genomes = n_distinct(GENOME1),
                                        reads = n()))
    
    # vector holding barcodes in > 1 genomes for filtering
    bad_bc <- as.character(subset(bc2genome, genomes > 1)$barcodes)
    
    # subset and filter bad barcodes
    tmp_bad_bc <- subset(tmp_filt_table, barcodes %in% bad_bc)
    tmp_filt_table <- subset(tmp_filt_table, barcodes %notin% bad_bc)
    
    
    # Write filtered read and barcode data to disk 
    write.table(tmp_filt_table, file = paste0(out_dir,"/hits_filt/jm/",sample_names[i],".hits.filt"), quote = F, row.names = F, sep = "\t")
    write.table(tmp_bad_bc, file = paste0(out_dir,"/hits_filt/jm/",sample_names[i],".hits.badBC"), quote = F, row.names = F, sep = "\t")
    write.table(bc2genome, file = paste0(out_dir,"/hits_filt/jm/",sample_names[i],".hits.bc2gen"), quote = F, row.names = F, sep = "\t")
    
    
    
    ### Summarize Hit Data 
    
    ## Hit Summary Breakdown - Summarize hits
    tmp_filt_summary <- data.table(SAMPLE = sample_names[i],
                                   tmp_filt_table %>% 
                                     group_by(GENOME1,model) %>% ## This may need to be fixed for pe if genomes are not equal (decide which read better for GENOME call)
                                     summarise(RDS = n(),
                                               UBC = n_distinct(barcodes),
                                               RDSnR = sum(plyr::count(barcodes)[which(plyr::count(barcodes)[,2] >= bc_rep),2]),
                                               UBCnR = sum(plyr::count(barcodes)[,2] >= bc_rep),
                                               BCPR = n_distinct(barcodes) / n()))
                                               
    ## Hit Summary Breakdown - Calculate fractional Abundances
    tmp_filt_summary <- data.table(tmp_filt_summary,
                                   RDS_FRC = tmp_filt_summary$RDS/sum(tmp_filt_summary$RDS),
                                   UBC_FRC = tmp_filt_summary$UBC/sum(tmp_filt_summary$UBC),
                                   RDSnR_FRC = tmp_filt_summary$RDSnR/sum(tmp_filt_summary$RDSnR),
                                   UBCnR_FRC = tmp_filt_summary$UBCnR/sum(tmp_filt_summary$UBCnR))
    
    ## Hit Summary Breakdown - Normalize to spike in org #### POSSIBLE BUG IF SIO HITS MULTIPLE MODELS...USED MAX TO SELECT MOST ABUNDANT
    # store fractional sio data for current sample
    tmp_sio <- subset(tmp_filt_summary, GENOME1 == spike_in_org) 
    
    # normalize within sample data by dividing fractional abundances by sio fractional abundance
    tmp_filt_summary <- data.frame(tmp_filt_summary,
                                   RDS_NRM = tmp_filt_summary$RDS_FRC/max(tmp_sio$RDS_FRC),
                                   UBC_NRM = tmp_filt_summary$UBC_FRC/max(tmp_sio$UBC_FRC),
                                   RDS_NRM = tmp_filt_summary$RDSnR_FRC/max(tmp_sio$RDSnR_FRC),
                                   UBCnR_NRM = tmp_filt_summary$UBCnR_FRC/max(tmp_sio$UBCnR_FRC))
    
    
    
    ## Hit Summary Breakdown - Add counts of bad BC to genomes
    
    # count reads with bad bc going to each genome / model
    bad_bc_counts <- plyr::count(tmp_bad_bc, vars = c("GENOME1", "model"))
    
    # merge into tmp_filt_summary
    tmp_filt_summary <- merge(tmp_filt_summary, bad_bc_counts, by = c("GENOME1", "model"), all.x = T)

    # Create table of filt hits for all samples
    filt_hits <- rbind(filt_hits, tmp_filt_summary)
    
  }
  
  # Create ordered output with correct columns, complete observations, and zeros
  genomes <- unique(filt_hits$GENOME1)
  models <- unique(filt_hits$model)
  all_cat <- expand.grid(SAMPLE = sample_names, GENOME1 = genomes, model = models,
                         KEEP.OUT.ATTRS = F,
                         stringsAsFactors = F)
  
  filt_hits <-  merge(all_cat, filt_hits, by = c("SAMPLE", "GENOME1", "model"), all.x = T)
  filt_hits <- filt_hits[order(filt_hits$SAMPLE, -filt_hits$RDS, filt_hits$GENOME1, filt_hits$model),]
  filt_hits[is.na(filt_hits)] <- 0
  colnames(filt_hits)[c(2,3,17)] <- c("GENOME", "MODEL", "BADBCR") 
  
  # return output
  return(filt_hits)
  
} # *** UP TO DATE

## Function 2: Summarize hits per genomes across samples for single end data
jm.summary.se <- function(){
  
  ### Set up data strucutres
  
  # Make df to hold final hit table
  filt_hits <- data.table()
  sample_names <- sub(".hits","",input_data_names$jm_file_list)
  
  # Create directories to store filtered and bad hits
  dir.create(paste0(out_dir,"/hits_filt/jm"), recursive = T)
  
  # Log this analysis 
  cat(paste0("Summary Mode (jm): SE\n"),
      file = paste0(out_dir,"/run_log.txt"), append = T)
  
  ### Run Filtering and Summarization on each sample and combine
  
  # Build master df of all hits
  for (i in 1:length(sample_names)){
    
    # Say what is happening
    cat(paste0("Calculating Summary Stats for jm Sample: ",sample_names[i],"\n"))
    
    # Read in hit table
    tmp_hit_table <- fread(paste0(env_summary$jm_dir,"/hits/",sample_names[i],".hits"), sep = "\t")
    
    
    
    ### APPLY FILTERING (each filtering step is implemented individually for now so we can add options)
    
    # Filter by mapq and mismatch cutoff
    tmp_filt_table <- subset(tmp_hit_table, MAPQ >= mq_cut & NM <= mm_cut)
    
    # Filter out NA barcodes
    tmp_filt_table <- subset(tmp_filt_table, is.na(barcodes) == FALSE) 
  
    
    ## Multiple Barcode Hit Filtering: BC in > 1 genome filtered + written to disk 
    bc2genome <- data.table(tmp_filt_table %>%
                              group_by(barcodes) %>%
                              summarise(genomes = n_distinct(GENOME),
                                        reads = n()))
    
    
    # vector holding barcodes in > 1 genomes for filtering
    bad_bc <- as.character(subset(bc2genome, genomes > 1)$barcodes)
    
    # subset and filter bad barcodes
    tmp_bad_bc <- subset(tmp_filt_table, barcodes %in% bad_bc)
    tmp_filt_table <- subset(tmp_filt_table, barcodes %notin% bad_bc)
    
    
    # Write filtered read and barcode data to disk 
    write.table(tmp_filt_table, file = paste0(out_dir,"/hits_filt/jm/",sample_names[i],".hits.filt"), quote = F, row.names = F, sep = "\t")
    write.table(tmp_bad_bc, file = paste0(out_dir,"/hits_filt/jm/",sample_names[i],".hits.badBC"), quote = F, row.names = F, sep = "\t")
    write.table(bc2genome, file = paste0(out_dir,"/hits_filt/jm/",sample_names[i],".hits.bc2gen"), quote = F, row.names = F, sep = "\t")

    
    
    ### Summarize Hit Data 
    
    ## Hit Summary Breakdown - Summarize hits
    tmp_filt_summary <- data.table(SAMPLE = sample_names[i],
                                   tmp_filt_table %>% 
                                     group_by(GENOME,model) %>% ## This may need to be fixed for pe if genomes are not equal (decide which read better for GENOME call)
                                     summarise(RDS = n(),
                                               UBC = n_distinct(barcodes),
                                               UBCnR = sum(plyr::count(barcodes)[,2] >= bc_rep),
                                               BCPR = n_distinct(barcodes) / n()))
    
    ## Hit Summary Breakdown - Calculate fractional Abundances
    tmp_filt_summary <- data.table(tmp_filt_summary,
                                   RDS_FRC = tmp_filt_summary$RDS/sum(tmp_filt_summary$RDS),
                                   UBC_FRC = tmp_filt_summary$UBC/sum(tmp_filt_summary$UBC),
                                   UBCnR_FRC = tmp_filt_summary$UBCnR/sum(tmp_filt_summary$UBCnR))
    
    ## Hit Summary Breakdown - Normalize to spike in org #### POSSIBLE BUG IF SIO HITS MULTIPLE MODELS...USED MAX TO SELECT MOST ABUNDANT 
    # store fractional sio data for current sample
    tmp_sio <- subset(tmp_filt_summary, GENOME == spike_in_org) 
    
    # normalize within sample data by dividing fractional abundances by sio fractional abundance
    tmp_filt_summary <- data.frame(tmp_filt_summary,
                                   RDS_NRM = tmp_filt_summary$RDS_FRC/max(tmp_sio$RDS_FRC),
                                   UBC_NRM = tmp_filt_summary$UBC_FRC/max(tmp_sio$UBC_FRC),
                                   UBCnR_NRM = tmp_filt_summary$UBCnR_FRC/max(tmp_sio$UBCnR_FRC))
    
    ## Hit Summary Breakdown - Add counts of bad BC to genomes
    
    # count reads with bad bc going to each genome / model
    bad_bc_counts <- plyr::count(tmp_bad_bc, vars = c("GENOME", "model"))
    
    # merge into tmp_filt_summary
    tmp_filt_summary <- merge(tmp_filt_summary, bad_bc_counts, by = c("GENOME", "model"), all.x = T)
    
    # Create table of filt hits for all samples
    filt_hits <- rbind(filt_hits, tmp_filt_summary)

    
  }
  
  # Create ordered output with correct columns, complete observations, and zeros
  genomes <- unique(filt_hits$GENOME)
  models <- unique(filt_hits$model)
  all_cat <- expand.grid(SAMPLE = sample_names, GENOME = genomes, model = models,
                         KEEP.OUT.ATTRS = F,
                         stringsAsFactors = F)
  
  filt_hits <-  merge(all_cat, filt_hits, by = c("SAMPLE", "GENOME", "model"), all.x = T)
  filt_hits <- filt_hits[order(filt_hits$SAMPLE, -filt_hits$RDS, filt_hits$GENOME, filt_hits$model),]
  filt_hits[is.na(filt_hits)] <- 0
  colnames(filt_hits)[c(2,3,14)] <- c("GENOME", "MODEL", "BADBCR") 
  
  # return output
  return(filt_hits)
  
} # *** UP TO DATE
  
## Function 3: Summarize Coverage per genomes across samples for paired end data
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


  
} # *** NOT CHECKED

## Function 4: Summarize Coverage per genomes across samples for single end data
lm.summary.se <- function(){
  
} # *** NOT CHECKED

## Function 5: Calculates summary statistics for each total sample across all hits
## RETURN: samp_hit_smry
hit.summary <- function(){
  
  # Say what is happening
  cat("\nCalculating Overall Sample Summaries...\n")
  
  # Gather data and sample number
  tot_reads <- fread(paste0(env_summary$jm_dir,"/jm_workflow_stats.txt"))[,c(1,5)]
  raw_dat <- ETstats_jm_out
  samp_names <- unique(raw_dat$SAMPLE)
  
  # Create output data frame
  samp_hit_summary <- data.table()
  
  # Loop over each sample and generate output
  for (i in 1:length(samp_names)){
    
    # Make df to store tmp data
    tmp_hit_summary <- data.table()
    
    # Subset all data for specific sample
    tmp_sample_dat <- subset(raw_dat, SAMPLE == samp_names[i])
    
    # Get total reads for that sample by matching names between tot_reads and sample subset
    tmp_tot_reads <- as.numeric(tot_reads[which(tot_reads$V1 == samp_names[i]), 2])
    
    # Subset data for vector, donor, and sio if present
    tmp_vec_dat <- subset(tmp_sample_dat, GENOME %in% vector_names)
    tmp_donor_dat <- subset(tmp_sample_dat, GENOME == donor_org)
    tmp_sio_dat <- subset(tmp_sample_dat, GENOME == spike_in_org)

    
    # Temporary output data frame
    tmp_hit_summary <- data.table(Sample = samp_names[i],
                                  Total_Reads = tmp_tot_reads,
                                  Map = sum(tmp_sample_dat$RDS),
                                  Map_Frac = (sum(tmp_sample_dat$RDS) / tmp_tot_reads),
                                  Unq_BC = sum(tmp_sample_dat$UBC),
                                  Vec_Frac = (sum(tmp_vec_dat$RDS) / sum(tmp_sample_dat$RDS)),
                                  Vec_Unq = sum(tmp_vec_dat$UBC),
                                  Donor_Frac = (sum(tmp_donor_dat$RDS) / sum(tmp_sample_dat$RDS)),
                                  Donor_Unq = sum(tmp_donor_dat$UBC),
                                  SIO_Frac = (sum(tmp_sio_dat$RDS) / sum(tmp_sample_dat$RDS)),
                                  SIO_Unq = sum(tmp_sio_dat$UBC))
    
    tmp_hit_summary$Com_Reads <- (tmp_hit_summary$Map - sum(tmp_vec_dat$RDS,tmp_donor_dat$RDS,tmp_sio_dat$RDS))/tmp_hit_summary$Map
    tmp_hit_summary$Com_Unq <- tmp_hit_summary$Unq_BC - (tmp_hit_summary$Vec_Unq + tmp_hit_summary$Donor_Unq + tmp_hit_summary$SIO_Unq)
    tmp_hit_summary$Com_UnqvEff <- tmp_hit_summary$Com_Unq / tmp_hit_summary$Total_Reads
    tmp_hit_summary$Com_ReadsvEff <- tmp_hit_summary$Com_Reads / tmp_hit_summary$Total_Reads
    
    # Bind single sample data into final df
    samp_hit_summary <- rbind(samp_hit_summary, tmp_hit_summary)
    
  }
  
  # Return final df of all sample-level hit summaries
  return(samp_hit_summary)
} # *** UP TO DATE

## Function 6: Calculate frequency of on/off target hits to specific sites
## RETURN: 
ts.summary.se <- function(){
  
  ### Set up data strucutres
  
  # Make dfs to hold final hit table and summary counts
  targ_hits_summary <- data.table()
  basic_read_summary <- data.table()
  
  # get sample names
  sample_names <- sub(".hits.filt","",input_data_names)
  
  # Create directories to store filtered and bad hits
  dir.create(paste0(out_dir,"/hits_filt/"), recursive = T)
  
  # Name columns of protospacer target table
  colnames(sample_info) <- c("GENOME", "SCAF", "TARG", "STRAND", "START", "END", "NEUT")
  
  # Add start & end windows for target site based on strand  to target table ( + window bp past end of site)
  sample_info$STARTW <- ifelse(sample_info$STRAND == "+",sample_info$START, sample_info$START - window) 
  sample_info$ENDW <- ifelse(sample_info$STRAND == "+", sample_info$END + window, sample_info$END)
  

  
  ### Identify each target(s) in each filtered hits file
  
  # Build master df of all hits
  for (i in 1:length(sample_names)){
    
    # Make df to hold all sample hits
    sample_targ_hits <- data.table()
    
    # Say what is happening
    cat(paste0("Identifying targeted inserts for sample: ",sample_names[i],"\n"))
    
    # Read in filtered hit table from ETstats - ss output
    raw_hit_table <- fread(paste0(env_summary$ss_hits,"/",sample_names[i],".hits.filt"), sep = "\t")
  
    # Subset out reads connected to targeted trasnposons and store numbers identified in df
    sample_model_hits <- subset(raw_hit_table, model %in% Tn_names)
    Targ_Tn_counts <- plyr::count(sample_model_hits,vars = "model")
    Reg_Tn_count <- nrow(raw_hit_table) - sum(Targ_Tn_counts$freq)
    
    # Store basic read count summary
    tmp_read_summary <- data.table(SAMPLE = sample_names[i],
                                   MOD_RDS = sum(Targ_Tn_counts$freq),
                                   REG_RDS = Reg_Tn_count,
                                   TOT_RDS = nrow(raw_hit_table))
    
    # Merge to basic read summary for each sample
    basic_read_summary <- rbind(basic_read_summary, tmp_read_summary)
    
    
    ## Run for loop to look at each target gene in each organism
    for (j in 1:nrow(sample_info)){
      
      # store filters for target iteration j
      filters <- sample_info[j,]
      
      # Say what is happening
      cat(paste0("   Quantifying ",filters$GENOME," ",filters$TARG,"\n")) ### THIS MAY HAVE INCOMPATIBILITY WITH PE DATA
      
      # Subset raw hits by genome and scaffold that match target j
      tmp_targ_hits <- subset(sample_model_hits, GENOME == filters$GENOME & SCAF == filters$SCAF) ### THIS MAY HAVE INCOMPATIBILITY WITH PE DATA
      
      # Establish if putative inserts are same or reverse orientation to spacer (ORE = Oreintation relative to spacer) ### THIS IS SPECFIC FOR REVERSE READS
      tmp_targ_hits$ORE <- ifelse(tmp_targ_hits$STRAND == filters$STRAND, "REV", "SAME")
      
      # Establish single position for end of reverse read closest to model window (IC = Insertion Coordinate) ### THIS MAY HAVE INCOMPATIBILITY WITH PE DATA
      tmp_targ_hits$IC <- ifelse((filters$STRAND == "+" & tmp_targ_hits$STRAND == "-") | (filters$STRAND == "-" & tmp_targ_hits$STRAND == "-"),
                                 tmp_targ_hits$START, tmp_targ_hits$END) 
      
      # Identify and subset reads with IC bounded by target window
      tmp_targ_hits <- subset(tmp_targ_hits, between(IC, filters$STARTW, filters$ENDW) == TRUE)
      
      # Check if no targets were hit and if no hits skip to next target
      if (nrow(tmp_targ_hits) == 0){
        next
      }
      
      # Tidy up hits for target iteration j
      tmp_targ_hits <- data.table(TARG = filters$TARG, tmp_targ_hits)
      
      # bind to sample_targ_hits table
      sample_targ_hits <- rbind(sample_targ_hits, tmp_targ_hits)
      
    }
    
    
    
    ### Write targeted insertions to disk
    write.table(sample_targ_hits, file = paste0(out_dir,"/hits_filt/",sample_names[i],".targ.hits"), quote = F, row.names = F, sep = "\t")
   
 
    
    ### Summarize Hit Data for sample
    
    ## Check if no targets were hit, and if so skip add filler data and to next sample
    if(nrow(sample_targ_hits) == 0){
      next()
    }
    
    
    # Hit Summary Breakdown - Summarize hits
    tmp_hit_summary <- data.table(SAMPLE = sample_names[i],
                                   sample_targ_hits %>% 
                                     group_by(GENOME,TARG,ORE) %>%
                                     summarise(RDS = n(),
                                               UBC = n_distinct(barcodes),
                                               UBCnR = sum(plyr::count(barcodes)[,2] >= bc_rep)))
    
    # Merge tmp_hit_summary into master hit summary across all samples
    targ_hits_summary <- rbind(targ_hits_summary, tmp_hit_summary)
    
    
    # Log this analysis ### DO SOME LOGGING BRO
    #cat(paste0("Summary Mode (jm): PE\n"),
    #    file = paste0(out_dir,"/run_log.txt"), append = T)
    
  }
  
  ### Create Final Outputs
  
  # Collect all possible genomes and targets
  genomes <- unique(sample_info$GENOME)
  targets <- unique(sample_info$TARG)
  orientation <- c("REV", "SAME")
  all_cat <- expand.grid(SAMPLE = sample_names, GENOME = genomes, TARG = targets, ORE = orientation ,
                         KEEP.OUT.ATTRS = F,
                         stringsAsFactors = F)
  
  # Merge target hit summary onto expanded grid
  targ_hits_summary <-  merge(all_cat, targ_hits_summary, by = c("SAMPLE", "GENOME", "TARG", "ORE"), all.x = T)
  
  # Merge basic read summary onto expanded grid
  targ_hits_summary <- merge(targ_hits_summary, basic_read_summary, by = "SAMPLE", all.x = T)
  
  # Order output and add zeros
  targ_hits_summary <- targ_hits_summary[order(targ_hits_summary$SAMPLE, -targ_hits_summary$RDS, targ_hits_summary$GENOME, targ_hits_summary$TARG, targ_hits_summary$ORE),]
  targ_hits_summary[is.na(targ_hits_summary)] <- 0
  
  # return output
  return(targ_hits_summary)
  
} # *** UP TO DATE

## Function 7: Calculate transformation efficency from combined jm/lm output 
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
    tmp_comb_dat <- merge(tmp_jm_dat, tmp_lm_dat, by.x = "GENOME1", by.y = "GENOME1", all.x = T) # This will cause a bug with SE samples
    
    # Get Spike in Org Data and remove from tmp_comb_dat
    sio_tmp_dat <- subset(tmp_comb_dat, GENOME1 == spike_in_org)
    tmp_comb_dat <- tmp_comb_dat[-which(tmp_comb_dat$GENOME1 == spike_in_org),]
    
    # Calculate Efficency statistics for each genome
    sio_uniq_flt <- as.numeric(sio_tmp_dat$UNIQ_FLT)
    sio_rpk_frac <- as.numeric(sio_tmp_dat$RPK_FLT_FRAC)
    sio_cov_frac <- as.numeric(sio_tmp_dat$COV_FLT_FRAC)
    tmp_eff_dat <- data.table(SAMPLE = tmp_comb_dat$SAMPLE.x,
                              GENOME1 = tmp_comb_dat$GENOME1, # This will cause a bug with SE samples 
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
  
} # *** NOT CHECKED

## Function 8: Read summary stats ### SPLIT THIS FUNCTION EVENTUALLY FOR READABLITY
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
  
} # *** NOT CHECKED

## Function X: Comparative Statistics between groups (hits or read summaries)

## Function X: Basic Plotting of genomes (genomes v sample groups)


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

  ## Create list of input data file names 
  input_data_names <- get.files()
  
  ## Perform and Output Statistical Summaries for jm Workflow if present
  if(env_summary$jm == TRUE){
    
    if(env_summary$jm_pe_map == TRUE){
      ETstats_jm_out <- jm.summary.pe()
    }
    
    if (env_summary$jm_pe_map == FALSE){
      ETstats_jm_out <- jm.summary.se()
    }
    
    # Calculate Sample-level hit summaries
    samp_hit_summary <- hit.summary()
      
    # Write Output
    write.table(ETstats_jm_out, file = paste0(out_dir,"/ETstats_jm_hit_summary.txt"), quote = F, row.names = F, sep = "\t")
    write.table(samp_hit_summary, file = paste0(out_dir,"/ETstats_jm_sample_summary.txt"), quote = F, row.names = F, sep = "\t")
    
  }
  
  ## Perform and Output Statistical Summaries for lm Workflow if present
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
  



### Targeted Insertion Summary Workflow 
if (wf == "ts"){
  
  ## Create Correct Dir Structure for Workflow
  out_dir <- paste0(out_dir,"/ts")
  dir.create(out_dir, recursive = T)
  
  ## Initalize program and create log file 
  make.log.file()
  
  ## Determine environment and important directories
  env_summary <- check.env()
  
  ## Create list of input data file names 
  input_data_names <- get.files()
  
  ## Perform and Output Statistical Summaries for ts Workflow
  targ_hits_summary <- ts.summary.se()
  
  ## Write Output
  write.table(targ_hits_summary, file = paste0(out_dir,"/ETstats_ts_targ_hits_summary.txt"), quote = F, row.names = F, sep = "\t")
  
  ## Congratulate
  cat("\nts Workflow completed successfully :-)\n")

  
  
}

### General Read Stats Summary Workflow ****NON FUNCTIONAL
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







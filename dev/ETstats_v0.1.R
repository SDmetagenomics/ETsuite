#!/usr/bin/env Rscript

### Set Version
etstats_version <- c("v0.10")



### Check and Load Libraries

## Check data.table
if ("data.table" %in% installed.packages() == FALSE){
  print("Please install R package data.table. Program quitting...")
  q(save="no")
}

## Check dplyr
if ("dplyr" %in% installed.packages() == FALSE){
  print("Please install R package dplyr. Program quitting...")
  q(save="no")
}


## Load packages
suppressMessages(library(data.table, quietly = T))
suppressMessages(library(dplyr, quietly = T))
suppressMessages(library(ggplot2, quietly = T))



### Set up path variables for associated scripts and databases

## Get relative path of ETmapper install directory
initial.options <- commandArgs(trailingOnly = F)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

## Set up default database paths (note: normalizePath creates a full path from relative path)
scripts <- normalizePath(paste0(script.basename,"/../scripts")) #accesory scripts directory



### Collect and parse arguments

## Create Args Variable
args <- commandArgs(trailingOnly = T)

## NOT IN Operator for Arg Parseing
'%notin%' <- Negate('%in%')

## Display help if no args or -h
if("-h" %in% args | !("-w" %in% args) | !("-d" %in% args) | length(args) == 0) {
  cat("
  
      ####### #######                                  
      #          #     ####  #####   ##   #####  ####  
      #          #    #        #    #  #    #   #      
      #####      #     ####    #   #    #   #    ####  
      #          #         #   #   ######   #        # 
      #          #    #    #   #   #    #   #   #    # 
      #######    #     ####    #   #    #   #    ####   v0.10
      
      Usage: ETstats.R -w [workflow] -d [Data_dir] [options]
      By: Spencer Diamond (sdiamond@berkeley.edu)
      
      Mandatory Arguments:
      
      -w: Workflow type (No Default)
          hc - hit count summary
          ts - targeted insertion summary
          cs - comparative stats summary

      -d: Directory with data for analysis (No Default; ET_Seq)

      Hit Count Summary Options:
      
      -H: Hamming distance used by bartender for clustering (Default: 3)
      -R: Multiread threshold for unique barcode detection (Default: 2)
      
      Targeted Insertion Summary Options:
      
      -W: Positive target window (Defalut: 200bp)
      -T: Targeted transposon models (Default: 'VcCas_model'; must be quoted)
      
      Program Control Options:
      
      -o: Output dir (Default: Data_dir/[workflow_type])
      -cpu: Number of cores (Default: 1)
      -h: Bring up this help menu
  
      
    ")
  
  q(save="no")
}




### Define Arguments

## Mandatory Arguments

# Work Flow Type
wf <- args[which(args == "-w") + 1]
if(wf %notin% c("hc","ts","cs")){
  cat(paste0("\n",wf," is not a known workflow...exiting"))
  q(save="no")
}

# ET_Seq working directory
et_dir <- args[which(args == "-d") + 1]
et_dir <- normalizePath(et_dir)



## Hit Count Summary Options

# Number of good reads to count a barcode in a genome (Default: 2)
h_dist <- 3
if("-H" %in% args){
  h_dist <- as.numeric(args[which(args == "-H") + 1])
}

# Number of good reads to count a barcode in a genome (Default: 2)
bc_rep <- 2
if("-R" %in% args){
  bc_rep <- as.numeric(args[which(args == "-R") + 1])
}



## Targeted Insertion Options

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



## Program Control Options

# Output directory (Will be created if not specified)
if("-o" %in% args){
  out_dir <- args[which(args == "-o") + 1]
  dir.create(out_dir, recursive = T)
} else{
  out_dir <- paste0(et_dir,"/",wf)
  dir.create(out_dir, recursive = T)
}

# number of cores to use where applicable
cpu <- 1
if("-cpu" %in% args){
  cpu <- as.numeric(args[which(args == "-cpu") + 1])
}


## Set data.table Threads
setDTthreads(cpu)





#### BEGIN FUNCTION DEFINE ####


### Setup Functions

## Function 1: Make a basic log file that will be appended based on analysis
make.log.file <- function(){
  
  cat(paste0("ETstats ",etstats_version,"   Beginning Analysis...\n\n"))

  cat(paste0("ETstats ",etstats_version," Summary    Created: ", date(),"\n\n"),
      paste0("Program Parameters:\n"),
      paste0("Arguments: "), args,"/n", 
      paste0("Workflow type is: ", wf),"\n",
      file = paste0(out_dir,"/run_log.txt"))
}

## Function 2: Identify the things you need to run a workflow in the ET_Seq Folder
## Return: env_summary
check.env <- function(){
  
  ## Env check for hc workflow
  if (wf == "hc"){
  
    # Make data.frame to hold environment variables
    env_summary <- data.frame(jm_dir = "A1",
                              jm_gd = "A2",
                              jm_map_mode = "A3",
                              batch_file = "A4",
                              stringsAsFactors = F)
  
    # 1 check ETmapper out for jm folder, store location, and output to log
    jm_exist <- dir.exists(paste0(et_dir,"/jm"))
    
    # IF it can find the jm folder then do this stuff:
    if (jm_exist == TRUE){
      
      # Say waht is happening 
      cat("Junction Mapping Directory Found...\n")
      
      # find all the directories of the files we need in the log file 
      env_summary$jm_dir <- normalizePath(paste0(et_dir,"/jm"))
      env_summary$jm_gd <- system(paste0("grep 'Genome Database' ",et_dir,"/jm/run_log.txt | awk '{print $3}'"), intern = T)
      env_summary$jm_map_mode <- system(paste0("grep 'Mapping Mode:' ",et_dir,"/jm/run_log.txt | awk '{print $3}'"), intern = T)
      env_summary$batch_file <- system(paste0("grep 'Batch File:' ",et_dir,"/jm/run_log.txt | awk '{print $3}'"), intern = T)
      
      # Log these locations in the run_log file
      cat(paste0("Found ETmapper jm Output: TRUE\n"),
          paste0("ETmapper jm Output Dir: ",env_summary$jm_dir,"\n"),
          paste0("ETmapper jm Batch File: ",env_summary$batch_file,"\n"),
          paste0("ETmapper jm Genome DB: ",env_summary$jm_gd,"\n"),
          paste0("ETmapper jm mapping PE: ",env_summary$jm_map_mode,"\n"),
          file = paste0(out_dir,"/run_log.txt"), append = T)
      
      
    }
    
    # If jm not found throw error and quit
    if (jm_exist == FALSE){
      cat("Could not find output directory for jm workflow in ETmapper output.\n
          Program Quitting...")
      q(save = "no")
    }
    
    # Output the env_summary data.frame
    return(env_summary)
    
  }
  
  
  ## Env check for ts workflow
  if (wf == "ts"){

    # Make data.frame to hold environment variables
    env_summary <- data.frame(hc = FALSE,
                              hc_dir = "A1",
                              hc_hits = "A2",
                              hc_mode = "A3",
                              stringsAsFactors = F)
    
    # Check ETmapper out for jm folder, store location, and output to log
    hc_exist <- dir.exists(paste0(et_dir,"/hc"))
    
    if (hc_exist == TRUE){
      cat("Summary Stats Dir Found...\n")
      env_summary$hc <- TRUE
      env_summary$hc_dir <- normalizePath(paste0(et_dir,"/hc"))
      env_summary$hc_hits <- normalizePath(paste0(et_dir,"/hc/hits_filt/jm"))
      env_summary$hc_mode <- system(paste0("grep 'Summary Mode (jm)' ",et_dir,"/hc/run_log.txt | awk '{print $4}'"), intern = T)
      cat(
        paste0("ETstats hc Output Dir: ",env_summary$hc_dir,"\n"),
        paste0("ETstats hc Filt Hits Dir: ",env_summary$hc_hits,"\n"),
        paste0("ETstats hc Mode: ",env_summary$hc_mode,"\n"),
        file = paste0(out_dir,"/run_log.txt"), append = T)
    }
    if (hc_exist == FALSE){
      cat("Could not find output directory for ETstats hc.\n
          Program Quitting...")
      q(save = "no")
    }
    
    # Output the env_summary data.frame
    return(env_summary)
  }
  
}


## Function 3: Clean up



#### END FUNCTION DEFINE ####






###############
############### HIT COUNT SUMMARY WORKFLOW
###############


if (wf == "hc"){
  
  ### Make Log File
  make.log.file()
  
  
  
  
  ### Load Needed Files to Begin Analysis
  
  ## Create env_summary table that has important directories 
  env_summary <- check.env()
  
  ## Load Batch File In
  batch_file <- fread(env_summary$batch_file, header = T, stringsAsFactors = F)
  
  ## Load Genome Stats Summary In
  genome_stats <- fread(paste0(env_summary$jm_gd,"/genome_stats.txt"), header = T, stringsAsFactors = F)
  
  ## Identify Hit File Names
  hit_tables <- list.files(env_summary$jm_dir,"/hits/*.hits")
  
  
  
  
  ### Create All Sample BC Table and Save Barcodes to Disk in loop
  
  ## Set up dfs to save output
  all_sample_bc <- data.table()
  
  
  ## Log Number of Samples For Processing
  num_samp <- length(hit_tables)
  
  cat(paste0("Number of Samples: ",num_samp),
      file = paste0(out_dir,"/run_log.txt"))

  
  ## Loop Over Each Hit Table, Collect Barcodes, and Summarize Hits
  for (i in 1:length(hit_tables)){
    
    # Get sample name
    sample <- sub(".hits","", hit_tables[i])
    
    # Say what is happening
    cat(paste0("Aggregating Barcodes for Sample: ",sample,"\n"))
    
    # Read in hit table
    tmp_hit_table <- fread(paste0(env_summary$jm_dir,"/hits/",raw_files[i]), header = T, sep = "\t")
    
    # Summarize hits by unique barcode, genome, junction, model combo
    bc2genome <- data.table(SAMPLE = sample,
                            tmp_hit_table %>%
                              group_by(barcodes, GENOME, TNjunc, model) %>%
                              summarise(reads = n()))
    
    # bind unique barcodes per sample to all_sample_bc
    all_sample_bc <- rbind(all_sample_bc, bc2genome)
    
    # Save all barcodes from sample to disk for clustering
    barcode_save <- data.table(tmp_filt_table$barcodes)
    fwrite(barcode_save, paste0(out_dir,"/",sample,".bc"), sep = ",", col.names = F)
    
  } 
  
  
  ## Remove un-needed variables
  rm(tmp_hit_table,bc2genome,barcode_save)
  
  
  
  
  ### Create file with all barcodes in Bartender format and Execute barcode clustering with Bartender
  
  ## Say what is happening
  cat("\nCreating Bartender Input File...\n")
  
  ## Concatenate all barcodes to single file
  system(paste0("cat ",out_dir,"/*.bc > ",out_dir,"/all_bc.txt"))
  
  ## Generate sequential numbering for all BCs
  system(paste0("nl ",out_dir,"/all_bc.txt | awk '{print $1}' > ",out_dir,"/lnumbs"))
  
  ## Merge barcodes and numbers
  system(paste0("paste -d , ",out_dir,"/all_bc.txt ",out_dir,"/lnumbs > ",out_dir,"/all_bc_BTinput.txt"))
  
  ## Remove unneeded files
  system(paste0("rm ",out_dir,"/all_bc.txt ",out_dir,"/lnumbs"))
  
  ## Run Bartender Clustering Command
  
  # Say what is happening
  cat("\nClustering Barcodes...\n")

  # Run Bartender
  system(paste0("bartender_single_com",
               " -f ",out_dir,"/all_bc_BTinput.txt",
               " -o ",out_dir,"/all_bc",
               " -l 4",
               " -t ",cpu,
               " -s 1",
               " -d ",h_dist,
               " > ",out_dir,"/all_bc_bartender.log"))
  
  # Say what is happening
  cat("\nBarcode Clustering Complete...\n")
  
  
  
  
  ### Load in clustered barcodes and calculate summary statistics
  
  ## Say what is happening
  cat("\nCreating Barcode Cluster Summaries...\n")
  
  ## Import bartender clustering data         ############################# THIS IS WHERE I STOPPED
  bc_clust <- fread(paste0(out_dir,"/all_bc_barcode.csv"))
  colnames(bc_clust) <- c("barcodes", "reads","clstID")
  
  
  ## Perform some summary statistics + plots on bc_clustering
  
  bc_clust_summary <- data.table(bc_clust %>%
                                   group_by(clstID) %>%
                                   summarise(SIZE = n(),
                                             RDS = sum(reads),
                                             MAX_RDS = max(reads),
                                             RST_RDS = RDS - MAX_RDS,
                                             PRIME_BC = barcodes[which.max(reads)],
                                             PRIME_FRC = MAX_RDS/RDS))
  
  bc_general_stats <- data.table(Total_Clusters = nrow(bc_clust_summary),
                                 Total_Variants = sum(bc_clust_summary$SIZE),
                                 Total_Reads = sum(bc_clust_summary$RDS),
                                 Reads_in_Prime = sum(bc_clust_summary$MAX_RDS),
                                 Frac_in_Prime = sum(bc_clust_summary$MAX_RDS) / sum(bc_clust_summary$RDS),
                                 Clust_5RDS = nrow(subset(bc_clust_summary, RDS >= 5)),
                                 Variants_5RDS = sum(subset(bc_clust_summary, RDS >= 5)$SIZE),
                                 Frac_Variants_5RDS = sum(subset(bc_clust_summary, RDS >= 5)$SIZE) / sum(bc_clust_summary$SIZE),
                                 Reads_5RDS = sum(subset(bc_clust_summary, RDS >= 5)$RDS),
                                 Frac_Reads_5RDS = sum(subset(bc_clust_summary, RDS >= 5)$RDS) / sum(bc_clust_summary$RDS))
  
  ## Output Summary Data
  fwrite(bc_general_stats, paste0(out_dir,"/all_bc_cluster_stats.txt", col.names = T, sep = "\t"))
  
  cat("\nDONE FOR NOW...\n")
  
  
}








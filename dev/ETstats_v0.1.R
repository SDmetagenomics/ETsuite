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
      -C: Number of reads to confirm a BC cluster (Default: 5)
      -F: Barcode cluster filtering level (Default: high)
      -S: Emperically determined read swap rate (Default: Auto)
      -R: Min reads in sample to count barcode (Default: 2)
      
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

# Hamming distance used by bartender for clustering (Default: 3)
h_dist <- 3
if("-H" %in% args){
  h_dist <- as.numeric(args[which(args == "-H") + 1])
}

# Number of reads to count a BC cluster (Default: 5)
bc_r_count <- 5
if("-C" %in% args){
  bc_r_count <- as.numeric(args[which(args == "-C") + 1])
}

# Barcode cluster filtering level (Default: high)
bc_filt_level <- "high"
if("-F" %in% args){
  bc_filt_level <- as.character(args[which(args == "-F") + 1])
}

# Swap rate calculation (Default: Auto)
swap_rate_method <- "Auto"
if("-S" %in% args){
  swap_rate_method <- as.numeric(args[which(args == "-S") + 1])
}

# Min reads in sample to count barcode (Default: 2)
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
      paste0("Arguments: "), args, 
      paste0("\nWorkflow type is: ", wf,"\n"),
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
  hit_tables <- list.files(paste0(env_summary$jm_dir,"/hits/"))
  
  
  
  
  ### Create All Sample BC Table and Save Barcodes to Disk in loop
  
  ## Set up dfs to save output
  all_sample_bc <- data.table()
  
  
  ## Log Number of Samples For Processing
  num_samp <- length(hit_tables)
  
  cat(paste0("Number of Samples: ",num_samp,"\n"),
      file = paste0(out_dir,"/run_log.txt"),
      append = T)

  
  ## Loop Over Each Hit Table, Collect Barcodes, and Summarize Hits
  for (i in 1:length(hit_tables)){
    
    # Get sample name
    sample <- sub(".hits","", hit_tables[i])
    
    # Say what is happening
    cat(paste0("Aggregating Barcodes for Sample: ",sample,"\n"))
    
    # Read in hit table
    tmp_hit_table <- fread(paste0(env_summary$jm_dir,"/hits/",hit_tables[i]), header = T, sep = "\t")
    
    # Summarize hits by unique barcode, genome, junction, model combo
    bc2genome <- data.table(SAMPLE = sample,
                            tmp_hit_table %>%
                              group_by(barcodes, GENOME, TNjunc, model) %>%
                              summarise(reads = n()))
    
    # bind unique barcodes per sample to all_sample_bc
    all_sample_bc <- rbind(all_sample_bc, bc2genome)
    
    # Save all barcodes from sample to disk for clustering
    barcode_save <- data.table(tmp_hit_table$barcodes)
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
  
  ## Say what is happening
  cat("\nClustering Barcodes...\n")

  ## Run Bartender
  system(paste0("bartender_single_com",
               " -f ",out_dir,"/all_bc_BTinput.txt",
               " -o ",out_dir,"/all_bc",
               " -l 4",
               " -t ",cpu,
               " -s 1",
               " -d ",h_dist,
               " > ",out_dir,"/all_bc_bartender.log"))
  
  ## Say what is happening
  cat("\nBarcode Clustering Complete...\n")
  
  ## Log Stuff About Bartender Options
  cat(paste0("Bartender Dist: ",h_dist,"\n"),
      file = paste0(out_dir,"/run_log.txt"),
      append = T)
  
  
  
  
  ### Load in clustered barcodes and calculate summary statistics
  
  ## Say what is happening
  cat("\nCreating Barcode Cluster Summaries...\n")
  
  ## Import bartender clustering data
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
  fwrite(bc_general_stats, paste0(out_dir,"/all_bc_cluster_stats.txt"), col.names = T, sep = "\t")
  fwrite(bc_clust_summary, paste0(out_dir,"/all_bc_cluster_summary.txt"), col.names = T, sep = "\t")

  
  
  
  ### Assignment of clusters to barcodes in all sample BC master list + Calculation of Barcode Cluster Purity
  
  ## Say what is happening
  cat("\nMerging BC Clusters with Samples and Calculating Purity...\n")
  
  ## Assign Clusters to all_sample_bc by merging in bc_clust based on barcode cluster *** THIS DATA SHOULD BE OUTPUT
  all_sample_bc_clust <- merge(all_sample_bc, bc_clust[,c("barcodes","clstID")], by = "barcodes", all.x = T)
  
  ## Summarize prevelance of genomes and mapping position by barcode cluster *** THIS DATA SHOULD BE OUTPUT 
  master_clust_assignment <- data.table(all_sample_bc_clust %>%
                                          group_by(clstID) %>%
                                          summarise(RDS = sum(reads), # Sum of all reads assigned to BC cluster
                                                    NUM_BC = n_distinct(barcodes), # Number of Barcode variants seen in BC cluster
                                                    NUM_SAMP = n_distinct(SAMPLE), # Number of Samples seen in BC cluster
                                                    NUM_GEN = n_distinct(GENOME), # Number of Genomes seen in BC cluster 
                                                    NUM_POS = n_distinct(TNjunc), # Number of unique mapping positions seen in BC cluster
                                                    NUM_MOD = n_distinct(model), # Number of unique transposon models seen in BC cluster
                                                    PRIME_SAMP = SAMPLE[which.max(reads)], # Sample in sample/genome/barcode/junction/model combination that has the most reads
                                                    PRIME_SAMP_RDS = sum(reads[which(SAMPLE == PRIME_SAMP)]), # Sum of all reads for above genome in cluster 
                                                    PRIME_SAMP_FRC = PRIME_SAMP_RDS / RDS, # Fraction of reads going to Prime Sample for a cluster
                                                    PRIME_GEN = GENOME[which.max(reads)], # Genome in sample/genome/barcode/junction/model combination that has the most reads
                                                    PRIME_GEN_RDS = sum(reads[which(GENOME == PRIME_GEN)]), # Sum of all reads for Prime Genome in cluster
                                                    PRIME_GEN_FRC = PRIME_GEN_RDS / RDS, # Fraction of reads in a cluster coming from Prime Genome
                                                    PRIME_POS = TNjunc[which.max(reads)], # Position in sample/genome/barcode/junction/model combination that has the most reads
                                                    PRIME_POS_RDS = sum(reads[which(TNjunc == PRIME_POS)]), # Sum of all reads for the Prime Position in cluster
                                                    PRIME_POS_FRC = PRIME_POS_RDS / RDS, # Fraction of reads in a clsuter coming from Prime Position
                                                    PRIME_POS3_RDS = sum(reads[which(TNjunc %in% (PRIME_POS-3):(PRIME_POS+3))]), # Sum of all reads within 3bp of Prime Position
                                                    PRIME_POS3_FRC = PRIME_POS3_RDS / RDS, # Fraction of reads in a clsuter coming from within 3bp of Prime Position
                                                    PRIME_MOD = model[which.max(reads)],
                                                    PRIME_MOD_RDS = sum(reads[which(model == PRIME_MOD)]),
                                                    PRIME_MOD_FRC = PRIME_MOD_RDS / RDS,
                                                    PRIME_GS_RDS = sum(reads[which(GENOME == PRIME_GEN & SAMPLE == PRIME_SAMP)]), ## Number of reads in cluster where genome is prime genome and sample is prime sample
                                                    PRIME_GS_FRC = PRIME_GS_RDS / RDS,
                                                    PRIME_GNS_RDS = sum(reads[which(GENOME == PRIME_GEN & SAMPLE != PRIME_SAMP)]))) ## Number of reads in cluster where genome is prime genome and sample is NOT prime sample
  
  ## Generate and output summary statistics on overall purity of barcode clusters
  master_clust_summary <- data.table(master_clust_assignment %>%
                                       summarise(Total_Clusters = n(),
                                                 SAMP_n2PLUS = sum(NUM_SAMP >= 2),
                                                 SAMP_n2PLUS_FRC = SAMP_n2PLUS / Total_Clusters,
                                                 GEN_n2PLUS = sum(NUM_GEN >= 2),
                                                 GEN_n2PLUS_FRC = GEN_n2PLUS / Total_Clusters,
                                                 POS_n2PLUS = sum(NUM_POS >= 2),
                                                 POS_n2PLUS_FRC = POS_n2PLUS / Total_Clusters,
                                                 MOD_n2PLUS = sum(NUM_MOD >= 2),
                                                 MOD_n2PLUS_FRC = MOD_n2PLUS / Total_Clusters,
                                                 FRC_SAMP_PURE = sum(PRIME_SAMP_FRC >= 0.75) / Total_Clusters,
                                                 FRC_GEN_PURE = sum(PRIME_GEN_FRC >= 0.75) / Total_Clusters,
                                                 FRC_POS_PURE = sum(PRIME_POS_FRC >= 0.75) / Total_Clusters,
                                                 FRC_POS3_PURE = sum(PRIME_POS3_FRC >= 0.75) / Total_Clusters,
                                                 FRC_MOD_PURE = sum(PRIME_MOD_FRC >= 0.75) / Total_Clusters,
                                                 FRC_GS_PURE = sum(PRIME_GS_FRC >= 0.75) / Total_Clusters,
                                                 FRC_METRIC1 = sum(PRIME_GEN_RDS >= bc_r_count & PRIME_GEN_FRC >= 0.75) / Total_Clusters,
                                                 FRC_METRIC2 = sum(PRIME_GEN_RDS >= bc_r_count & PRIME_GEN_FRC >= 0.75 & PRIME_MOD_FRC >= 0.75) / Total_Clusters,
                                                 FRC_METRIC3 = sum(PRIME_GEN_RDS >= bc_r_count & PRIME_GEN_FRC >= 0.75 & PRIME_POS3_FRC >= 0.75) / Total_Clusters))
  
  ## Output Summary Data
  fwrite(master_clust_summary, paste0(out_dir,"/all_bc_cluster_purity_stats.txt"), col.names = T, sep = "\t")
  fwrite(master_clust_assignment, paste0(out_dir,"/all_bc_cluster_master_assignment.txt"), col.names = T, sep = "\t")
  
  
  
  
  
  ### Filter Master Cluster Assignment and Calculate Barcode Swap Rate
  
  ## Say what is happening
  cat("\nFiltering Barcode Clusters...\n")
  
  ## Apply Cluster Filter 
  if(bc_filt_level == "high"){
    master_clust_assignment_filt <- subset(master_clust_assignment, PRIME_GEN_FRC >= 0.75 & PRIME_GEN_RDS >= bc_r_count & PRIME_MOD_FRC >= 0.75)
  }
  
  if(bc_filt_level == "medium"){
    master_clust_assignment_filt <- subset(master_clust_assignment, PRIME_GEN_FRC >= 0.75 & PRIME_GEN_RDS >= bc_r_count)
  }
  
  ## Caculate Filter Losses
  bc_filt_loss <- round((nrow(master_clust_assignment_filt) / nrow(master_clust_assignment)) * 100, digits = 2)
  cat(paste0(bc_filt_loss, "% of BC clusters remaining after filtering\n"))
  
  ## Log FIlter Used
  cat(paste0("Barcode Cluster Filter: ",bc_filt_level,"\n"),
      paste0(bc_filt_loss, "% of BC clusters remaining after filtering\n"),
      file = paste0(out_dir,"/run_log.txt"),
      append = T)
  
  
  
  
  
  ### Calculate barcode cluster swap rates
  # *** NOTE: Because a barcode cluster can only be assigned to one genome we calcualte the swap rate based on the counts of prime vs not-prime sample reads only for reads going to the primary genome.
  
  ## Say what is happening
  cat("\nCalculating and Applying BC Swap Rates...\n")
  
  ## Summary Table of Counts and Swap Rates for Each Org
  count_swap_summary <- data.table(master_clust_assignment_filt %>%
                                     group_by(PRIME_GEN) %>%
                                     summarise(ALL_RDS = sum(RDS),
                                               PRIME_GEN_READS = sum(PRIME_GEN_RDS),
                                               PRIME_SAMP_READS = sum(PRIME_SAMP_RDS),
                                               OUT_SAMP_READS = ALL_RDS - PRIME_SAMP_READS,
                                               PRIME_GS_READS = sum(PRIME_GS_RDS),
                                               PRIME_GNS_READS = sum(PRIME_GNS_RDS),
                                               S_RATE_ALL = OUT_SAMP_READS/ALL_RDS, # Swap rate based on all reads for a BC in prime sample vs out of prime sample
                                               S_RATE_PrimeOnly = PRIME_GNS_READS / PRIME_GEN_READS)) # Swap rate based on only PRIME GENOME reads for a BC PRIME GENOME in PRIME SAMPLE vs PRIME GENOME out PRIME SAMPLE *** This is the one we use
  
  ## Create data.frame of Genome Types (Target/Vector/Ect) Merge to count_swap_summary
  genome_type <- genome_stats[,c("Genome", "Type")]
  count_swap_summary <- merge(count_swap_summary, genome_type, by.x = "PRIME_GEN", by.y = "Genome", all.x = T)
  
  
  ## Determine How to Caluclate Swap Rate
  
  # IF swap reate Auto calculate from data
  if(swap_rate_method == "Auto"){
    
    # Subset only target organisms to calculate swap rate table
    count_swap_filt <- subset(count_swap_summary, Type == "Target")
    
    # Calculate in and out sample counts for things expected in only one sample 
    in_counts <- sum(count_swap_filt$PRIME_GS_READS)
    out_counts <- sum(count_swap_filt$PRIME_GNS_READS)
    
    # Add sd to mean to get high end swap rate
    s_rate <- out_counts / (in_counts + out_counts)
    
  }
  
  # IF number passed for swap rate use the number
  if(is.numeric(swap_rate_method) == T){
    
    # Say what is happening
    cat(paste0("\nProvided swap rate (",swap_rate_method,") applied...\n"))
    
    # Apply provided swap rate
    s_rate <- swap_rate_method
  }
  
  
  ## Log swap rate
  cat(paste0("Swap Rate Used: ",s_rate,"\n"),
      file = paste0(out_dir,"/run_log.txt"),
      append = T)
  
  ## Append number of reads per BC cluster expected to swap based on s_rate
  # Calculates the lower bound of the number of reads that could swap into any one sample based on calculated swap rate + 2 SD
  master_clust_assignment_filt <- data.table(master_clust_assignment_filt,
                                             SWAP_RATE = ceiling(master_clust_assignment_filt$PRIME_GEN_RDS*(s_rate / num_samp) + 2*(sqrt(master_clust_assignment_filt$PRIME_GEN_RDS*(1-s_rate)*s_rate))))
  
  ## Merge cluster assignments and swap rate (number of expected swaps) into bc_assign *** The merge will take the intersect and preserve only barcodes from clusters that pass filter
  bc_assign <- merge(bc_clust, master_clust_assignment_filt[,c("clstID", "PRIME_GEN", "PRIME_MOD", "SWAP_RATE")], by = "clstID")
  
  ## Write Output Data
  fwrite(count_swap_summary, paste0(out_dir,"/all_bc_swap_summary.txt"), col.names = T, sep = "\t")
  fwrite(master_clust_assignment_filt, paste0(out_dir,"/all_bc_cluster_master_assignment_filt.txt"), col.names = T, sep = "\t")
  fwrite(bc_assign, paste0(out_dir,"/all_bc_BC_ASSIGN.txt"), col.names = T, sep = "\t")
  
  
  

  ### Run For Loop to Analyze Each Individual Hit Table, Assign Clusters, and Count
  
  ## Say what is happening
  cat("\nCalculating Summary Counts for Each Sample...\n")
  
  ## Make df to hold final hit table + Gather All Sample Names
  filt_hits <- data.table()
  
  ## Loop Over Each Hit Table, Collect Barcodes, and Summarize Hits
  for (i in 1:length(hit_tables)){
    
    # Get sample name
    sample <- sub(".hits","", hit_tables[i])
    
    # Say what is happening
    cat(paste0("Processing Sample: ",sample,"\n"))
    
    # Read in hit table
    tmp_hit_table <- fread(paste0(env_summary$jm_dir,"/hits/",hit_tables[i]), header = T, sep = "\t")
    
    # Merge by bc string *** This should result in NAs for BCs that were removed
    tmp_hit_table <- merge(tmp_hit_table, bc_assign[,c("clstID","barcodes","PRIME_GEN","PRIME_MOD","SWAP_RATE")], by = "barcodes", all.x = T)
    
    # Remove reads where BC cluster did not pass filters or was not assigned
    tmp_hit_table <- subset(tmp_hit_table, is.na(clstID) == FALSE)
    
    # Remove reads that do not match BC PRIME_GEN & | PRIME_MOD depending on BC Cluster filter 
    if(bc_filt_level == "high"){
      tmp_hit_table <- subset(tmp_hit_table, GENOME == PRIME_GEN & model == PRIME_MOD)
    }
    if(bc_filt_level == "medium"){
      tmp_hit_table <- subset(tmp_hit_table, GENOME == PRIME_GEN)
    }

    # First aggregate number of reads by genome/model/BCcluster + include swap rate
    tmp_hit_summary <- data.table(SAMPLE = sample,
                                  tmp_hit_table %>% 
                                  group_by(GENOME,model,clstID,SWAP_RATE) %>%
                                  summarise(counts = n()))

    # Now calculate all base summary statistics that include counts of reads and unique barcode clusters
    tmp_hit_summary <- data.table(SAMPLE = sample,
                                  tmp_hit_summary %>% 
                                  group_by(GENOME,model) %>%
                                  summarise(RDS = sum(counts),
                                            UBC = n_distinct(clstID),
                                            RDSnR = sum(counts[counts >= bc_rep]),
                                            UBCnR = n_distinct(clstID[counts >= bc_rep]),
                                            RDSdC = sum(counts[counts >= 2 & counts > SWAP_RATE]),
                                            UBCdC = n_distinct(clstID[counts >= 2 & counts > SWAP_RATE]),
                                            BCPR = UBC / RDS,
                                            BCPRdC = UBCdC / RDSdC))
    
    # Create table of filt hits for all samples
    filt_hits <- rbind(filt_hits, tmp_hit_summary)
    
  } 
  
  
  
  
  ### Create ordered output with correct columns, complete observations, and zeros
  
  ## Say what is happening
  cat("\nProducing Output...\n")
  
  ## Create Expanded Grid
  genomes <- unique(genome_stats$Genome)
  models <- unique(filt_hits$model)
  sample_names <- unique(filt_hits$SAMPLE)
  all_cat <- expand.grid(SAMPLE = sample_names, GENOME = genomes, model = models,
                         KEEP.OUT.ATTRS = F,
                         stringsAsFactors = F)
  
  ## Grab all metadata from batch file
  sample_metadat <- batch_file[,c("SAMPLE", "DESC", "TF_TYPE", "EXP")]

  ## Merge Filtered Hits into Expanded Grid
  filt_hits <- merge(all_cat, filt_hits, by = c("SAMPLE", "GENOME", "model"), all.x = T)
  
  ## Merge Metadata into Expanded Grid
  filt_hits <- merge(filt_hits, sample_metadat, by = "SAMPLE", all.x = T)
  
  ## Order Output and add Zeros for NAs
  filt_hits <- filt_hits[order(filt_hits$SAMPLE, -filt_hits$RDS, filt_hits$GENOME, filt_hits$model),]
  filt_hits[is.na(filt_hits)] <- 0
  
  ## Identify Genomes that have all zero hits and remove to make a small version of filt_hits
  genome_sums <- aggregate(UBCdC ~ GENOME, filt_hits, sum)
  all_zero_genomes <- genome_sums$GENOME[which(genome_sums$UBCdC == 0)]
  filt_hits_small <- subset(filt_hits, GENOME %notin% all_zero_genomes)
  
  ## write outputs
  fwrite(filt_hits, paste0(out_dir,"/ETstats_hc_output.txt"), col.names = T, sep = "\t")
  fwrite(filt_hits_small, paste0(out_dir,"/ETstats_hc_output_small.txt"), col.names = T, sep = "\t")
  
  
}
  




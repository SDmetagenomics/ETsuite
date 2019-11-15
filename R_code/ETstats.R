#!/usr/bin/env Rscript

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

# Check ggforce
if ("ggforce" %in% installed.packages() == FALSE){
  print("Please install R package ggforce. Program quitting...")
  q(save="no")
}


### Load packages
suppressMessages(library(data.table, quietly = T))
suppressMessages(library(dplyr, quietly = T))
suppressMessages(library(ggplot2, quietly = T))
suppressMessages(library(ggforce, quietly = T))


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

print(args)
cat(paste0("\nTotal Arguments: ",length(args),"\n"))

## NOT IN Operator for Arg Parseing
'%notin%' <- Negate('%in%')

## Display help if no args or -h
if("-h" %in% args | !("-w" %in% args) | !("-d" %in% args) | !("-b" %in% args) | length(args) == 0) {
  cat("
      ETstats v0.02
      
      Usage: ETstats.R -w [workflow] -d [ETmapper_out_dir] -b [batch_file] [options]

      Mandatory Arguments:
      
      -w: Workflow type (No Default)
          ss - stats summary
          rs - read statistics
          ip - insertion plot
      -d: ETmapper output directory (No Default)
      -b: Batch file (No Default)
      
      Optional Arguments:
      
      -g: Directory containing genome database (No Default)
      -o: Output file (Default: ETstats_out.txt)
      -p: Make plots (Default: FALSE)
      -r: Remove specific genome from analysis (must be quoted)
      
      Stats Summary Options:
      
      -q: MapQ cutoff score (Default: 20)
      -m: Maximum read mismatches allowed (Default: 5; SE only)
      -C: Custom filter function (Overrides PE filters; must be quoted)
      -h: Bring up this help menu
  
  
      Example Usage:
      
      ## Stats summary
      ETstats.R -w ss -d ./ETmapper_out -b batch_file.txt
      ## Stats summary custom filter
      ETstats.R -w ss -d ./ETmapper_out -b batch_file.txt -c 'GENOME1 == GENOME2 & MAPQ1 > 5'
      
    ")
  
  q(save="no")
}


## Mandatory Arguments

# Work Flow Type
wf <- args[which(args == "-w") + 1]
if(wf %notin% c("ss","rs","ip")){
  cat(paste0("\n",wf," is not a kown workflow...exiting"))
  q(save="no")
}

# ETmapper output directory
em_out <- args[which(args == "-d") + 1]

# Batch File
bf <- args[which(args == "-b") + 1]
batch_file <- read.table(bf, sep = "\t", header = F)



## Optional Arguments

# Genome Database Directory
gd <- args[which(args == "-g") + 1]

# Output file / direcctory (eventually)
out_dir <- "ETstats_out.txt"
if("-o" %in% args){
  out_dir <- args[which(args == "-o") + 1]
}
#dir.create(out_dir, recursive = T)

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



## Stats summary options

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





# Load Testing Data
#batch_file <- read.table("test_data/ETmapper/test_batch.txt", header = F, sep = "\t", stringsAsFactors = F)
#em_out <- "test_data/ETmapper/full_paired_end_test/"
# pe_example <- fread("test_data/ETmapper/paired_end_test/hits/JD_ZZ_2ndcycle_1.hits")
# se_example <- fread("test_data/ETmapper/single_end_test/hits/JD_ZZ_2ndcycle_1.hits")
# 
# tmp_hit_table <- se_example
# i = 1
# mq_cut <- 20
# mm_cut <- 5





#### BEGIN FUNCTION DEFINE ####


## Function 1: Check if data processed as paired end or single end 
pe.test <- function(){

  if(ncol(batch_file) == 3){
    paired_end_data <- FALSE
    cat("Input Data Identifed as Single End...Begining Analysis\n\n")
    #cat(" Paired End = FALSE", file = paste0(out_dir,"/summary.txt"), append = T)
  }
  
  if(ncol(batch_file) == 4){
    paired_end_data <- TRUE
    cat("Input Data Identifed as Paired End...Begining Analysis\n\n")
    #cat(" Paired End = TRUE", file = paste0(out_dir,"/summary.txt"), append = T)
  } 
  
  # If tests fail exit with error
  if (exists("paired_end_data") == FALSE){
    cat("ERROR: Batch file not formatted correctly")
    q(save="no")
  }
  
  paired_end_data
  
}


## Function 2: Summarize hits per genomes across samples for paired end data
hit.summary.pe <- function(){
  
  # Make df to hold final hit table
  all_hits <- data.table()
  filt_hits <- data.table()
  
  # Build master df of all hits
  for (i in 1:nrow(batch_file)){
    
    # Read in hit table
    tmp_hit_table <- fread(paste0(em_out,"hits/",batch_file[i,1],".hits"), sep = "\t")
    
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

    
    # Sumarize raw data
    tmp_raw_summary <- data.table(SAMPLE = batch_file[i,1],
                                  GROUP = batch_file[i,2],
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
    tmp_filt_summary <- data.table(SAMPLE = batch_file[i,1],
                                   GROUP = batch_file[i,2],
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
  pe_hits_out <- merge(all_hits, filt_hits, by = c("SAMPLE", "GROUP", "GENOME1"), all.x = T)
  pe_hits_out <- pe_hits_out[order(SAMPLE, -UNIQ_FLT)]
  pe_hits_out
  
}


## Function 3: Summarize hits per genomes across samples for single end data
hit.summary.se <- function(){
  
  # Make df to hold final hit table
  all_hits <- data.table()
  filt_hits <- data.table()
  
  # Build master df of all hits
  for (i in 1:nrow(batch_file)){
    
    # Read in hit table
    tmp_hit_table <- fread(paste0(em_out,"hits/",batch_file[i,1],".hits"), sep = "\t")
    
    
    ### APPLY FILTERING 
    tmp_filt_table <- subset(tmp_hit_table, MAPQ >= mq_cut & NM <= mm_cut)
    
  
    # Sumarize raw data
    tmp_raw_summary <- data.table(SAMPLE = batch_file[i,1],
                                  GROUP = batch_file[i,2],
                                  tmp_hit_table %>% 
                                  group_by(GENOME) %>% 
                                  summarise(READ_RAW = n(),
                                            UNIQ_RAW = n_distinct(barcodes),
                                            BPR_RAW = n_distinct(barcodes) / n(),
                                            MMQ_RAW = mean(MAPQ),
                                            MMM_RAW = mean(NM),
                                            MLEN_RAW = mean(LEN)))
 
    # Sumarize filtered data
    tmp_filt_summary <- data.table(SAMPLE = batch_file[i,1],
                                   GROUP = batch_file[i,2],
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
  

## Function 4: Basic Plotting of genomes (genomes v sample groups)
  

## Function 5: Read summary stats ### SPLIT THIS FUNCTION EVENTUALLY FOR READABLITY
read.stats <- function(){
  
  # Make df and list to hold summary and full hit tables
  #aggregate_hits <- data.table() --- Don't do this right now
  read_summary <- data.table()
  
  # loop over batch file depending on data type
  if (paired_end_data == TRUE){
    
    # Build master df of all hits
    for (i in 1:nrow(batch_file)){
      
      # Read in hit table
      tmp_hit_table <- fread(paste0(em_out,"hits/",batch_file[i,1],".hits"), sep = "\t")
      
      # Remove genomes if in arguments
      if (exists("rm_gen") == TRUE){
        tmp_hit_table <- subset(tmp_hit_table, GENOME1 %notin% rm_gen[[1]] & GENOME2 %notin% rm_gen[[1]])
        
      }
        
      # Sumarize raw data on reads
      tmp_raw_summary <- data.table(SAMPLE = batch_file[i,1],
                                    GROUP = batch_file[i,2],
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
    for (i in 1:nrow(batch_file)){
      
      # Read in hit table
      tmp_hit_table <- fread(paste0(em_out,"hits/",batch_file[i,1],".hits"), sep = "\t")
      
      # Remove genomes if in arguments
      if (exists("rm_gen") == TRUE){
        tmp_hit_table <- subset(tmp_hit_table, GENOME1 %notin% rm_gen[[1]])
        
      }
      
      # Sumarize raw data
      tmp_raw_summary <- data.table(SAMPLE = batch_file[i,1],
                                    GROUP = batch_file[i,2],
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


if (wf == "ss"){
  
  # Check if data is pe or se
  paired_end_data <- pe.test()
  
  # loop over batch file depending on data type
  if (paired_end_data == TRUE){
    ETstats_out <- hit.summary.pe()
  }
  
  if (paired_end_data == FALSE){
    ETstats_out <- hit.summary.se()
  }
  
  # Write Output
  write.table(ETstats_out, file = out_dir, quote = F, row.names = F, sep = "\t")

}


if (wf == "rs"){
 
  # Check if data is pe or se
  paired_end_data <- pe.test()
  
  # loop over batch file depending on data type
  read_stats_out <- read.stats()
  
  # Write Output
  write.table(read_stats_out, file = out_dir, quote = F, row.names = F, sep = "\t")
 
   
}

if (wf == "ip"){
  
  cat("ip does not yet exist")
}







#! /usr/local/bin/Rscript

## TO DO
## 1. Add ability to set output dir
## 2. Create test suite
## 3. Add statistical Output with optional argument
## 4. This program should also do the length filtering on the mapped reads
## 5. The Program should be able to automatically find any accessory files it needs (ie: in its installed dir)




## Check and Load Libraries
if("data.table" %in% installed.packages() == FALSE){
  print("Please install R package data.table. Program quitting...")
  q(save="no")
}

library(data.table)



## Collect arguments
args <- commandArgs(trailingOnly = T)


## Display help if no args or -h
if("-h" %in% args | !("-i" %in% args) | length(args) == 0) {
  cat("
      Usage: ETstats_SR -i [hit_tables] [options]

      Mandatory arguments:
      -i: Batch file input with hit table paths, sample names, and group

      Optionnal arguments:
      -q: MapQ cutoff score (Default: 20)
      -m: Maximum read mismatches allowed (Default: 5)
      -h: Bring up this help menu\n\n")
  

  q(save="no")
}


## Parse Arguments

# Input file 
batch_file <- args[which(args == "-i") + 1]

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

###print(batch_file)
###print(mq_cut)
###print(mm_cut)


### Testing Batch File
### batch_file <- read.table("test_data/multiple_sample/test_batch.txt", header = F, sep = "\t", stringsAsFactors = F)
### i <- 1



## Load Batch File
batch_file <- read.table(batch_file, header = F, sep = "\t", stringsAsFactors = F)
colnames(batch_file) <- c("file", "sample_name", "group")



## ANALYSIS STEP 1 - Conduct Count Analysis for each Hit Table

# Create List Object to Hold All Samples 
Final_Results <- list()


for (i in 1:nrow(batch_file)){

  # Read in hit_table
  print(paste0("Processing: ",batch_file$sample_name[i]), quote = F)
  hit_table <- fread(batch_file$file[i], sep = "\t")
  colnames(hit_table) <- c("genome", "read", "mapq", "mis", "length", "pos")
  
  # capture all genome names in hit_table
  gen <- unique(hit_table$genome)
  
  # generate hit_table filtered for scores 
  filt_hit_table <- subset(hit_table, mapq >= mq_cut & mis <= mm_cut)

  
  # generate raw results dataframe
  raw_results <- data.frame()
  
  ## LOOP FOR RAW RESULTS
  for (j in 1:length(gen)){
    
    tmp_genome <- hit_table[hit_table$genome == gen[j],]
    
    # calculates stats for each genome - unique hits by unique positions 
    tmp_dat <- data.frame(Genome = as.character(gen[j]),
                          Reads_Raw = nrow(tmp_genome),
                          Uniq_Raw = length(unique(tmp_genome$pos)),
                          MMQ_Raw = mean(tmp_genome$mapq),
                          MMM_Raw = mean(tmp_genome$mis),
                          RF_Raw = nrow(tmp_genome) / nrow(hit_table))
    
    raw_results <- rbind(raw_results, tmp_dat)
  } 
  
  
  # generate filt results dataframe
  filt_results <- data.frame()
  
  ## LOOP FOR FILT RESULTS
  for (j in 1:length(gen)){
    
    tmp_genome <- filt_hit_table[filt_hit_table$genome == gen[j],]
    
    tmp_dat <- data.frame(Genome = as.character(gen[j]),
                          Reads_Filt = nrow(tmp_genome),
                          Uniq_Filt = length(unique(tmp_genome$pos)),
                          MMQ_Filt = mean(tmp_genome$mapq),
                          MMM_Filt = mean(tmp_genome$mis),
                          RF_Filt = nrow(tmp_genome) / nrow(filt_hit_table))
    
    filt_results <- rbind(filt_results, tmp_dat)
    
  } 
  
  
  ## Create dataframe to hold combined results for one sample and pass to list
  raw_results$Genome <- as.character(raw_results$Genome)
  filt_results$Genome <- as.character(filt_results$Genome)
  results_merged <- merge(raw_results, filt_results, by = "Genome")
  
  Final_Results[[i]] <- results_merged

}




## ANALYSIS STEP 2 - Create Final DataFrame with All Values 

## Collect All Unique Genome Names
uniq_gen <- vector()

for (i in 1:length(Final_Results)){
  tmp_gen <- Final_Results[[i]][[1]]
  uniq_gen <- c(uniq_gen, tmp_gen)
} 

uniq_gen <- unique(uniq_gen)


## Combine Metadata and count summaries (includes every genome per sample)
Hit_Summary <- data.frame()

for (i in 1:length(Final_Results)){
  
  tmp_meta <- data.frame(Genome = uniq_gen,
                         Sample = batch_file$sample_name[i],
                         Group = batch_file$group[i])
  
  tmp_merge <- merge(tmp_meta, Final_Results[[i]], by = "Genome", all.x = T)
  
  # Order dataframe by Uniq_Filt 
  tmp_merge <- tmp_merge[order(tmp_merge$Uniq_Filt, decreasing = T),]
  
  # Remove NAs and NaN from data
  #tmp_merge[is.na(tmp_merge)] <- 0
  
  # Concatenate rows into final master dataframe
  Hit_Summary <- rbind(Hit_Summary, tmp_merge)
  
}


## Write Count Summary Data
#out_dir <- getwd()
write.table(Hit_Summary, file = "Hit_Summary.txt", quote = F, row.names = F, sep = "\t")



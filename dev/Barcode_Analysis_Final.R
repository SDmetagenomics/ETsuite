library(data.table)
library(dplyr)
library(ggplot2)


### 
setDTthreads(6)
getDTthreads()


'%notin%' <- Negate('%in%')
input_dir <- "~/Desktop/ETseq_BC/Kleb_Curve_Experiment/ETmapper_raw_hits/"
output_dir <- "~/Desktop/ETseq_BC/Kleb_Curve_Experiment/barcode_clust/sample_bc/"
#input_dir <- "~/Desktop/ETseq_BC/NovaSeq_MultipleTransform/ETmapper_raw_hits/"
#output_dir <- "~/Desktop/ETseq_BC/NovaSeq_MultipleTransform/barcode_clust/sample_bc/"
setwd(input_dir)
raw_files <- list.files(input_dir)
#i = 6
mq_cut <- 20
spike_in_org <- "Bacteroides_thetaiotaomicron_VPI_5482_GCF_000011065.1"
bc_rep <- 2

### Set up dfs to save output
all_sample_bc <- data.table()
num_samp <- length(raw_files)


### Loop over files and identify all barcodes that are unique to organisims within samples,
### and make a master list of barcodes, the sample they came from, and what organism they were in
### across all samples 

for (i in 1:length(raw_files)){
  
  # Get sample name
  sample <- sub(".hits","", raw_files[i])
  
  # Say what is happening
  cat(paste0("Aggregating Barcodes from Sample: ",sample,"\n"))
  
  # Read in hit table
  tmp_hit_table <- fread(raw_files[i], sep = "\t")
  
  
  ### APPLY FILTERING (each filtering step is implemented individually for now so we can add options)
    
    # Hits to same genome
    tmp_filt_table <- subset(tmp_hit_table, GENOME1 == GENOME2)
    
    # Filter by mapq
    tmp_filt_table <- subset(tmp_filt_table, MAPQ1 >= mq_cut | MAPQ2 >= mq_cut)
    
    # Filter out NA barcodes
    tmp_filt_table <- subset(tmp_filt_table, is.na(barcodes) == FALSE)
    
    # Add Tn Junction Postion
    tmp_filt_table$TNjunc <- ifelse(tmp_filt_table$STRAND1 == "+", tmp_filt_table$START1, tmp_filt_table$END1)
    
  
  ### COLLECT SAMPLE SIZE INFORMATION
  
    
  ### IDENTIFY BARCODES IN READS AND WRITE TO DISK - all unique barcode/genome/junction combinations are parsed and counted
    
    # summarize by unique barcode, genome, junction combo
    bc2genome <- data.table(SAMPLE = sample,
                            tmp_filt_table %>%
                              group_by(barcodes, GENOME1, TNjunc) %>%
                              summarise(reads = n()))

    # bind unique barcodes per sample to all_sample_bc
    all_sample_bc <- rbind(all_sample_bc, bc2genome)
    
    # Save all barcodes from sample to disk for clustering
    barcode_save <- data.table(tmp_filt_table$barcodes)
    fwrite(barcode_save, paste0(output_dir,sample,".bc"), sep = ",", col.names = F)
    
    
} ### WE NOW HAVE all_sample_bc which has every possible combination of barcode/genome/TN_junction 



# rm to save mem
rm(barcode_save)
rm(tmp_filt_table)
rm(tmp_hit_table)
rm(bc2genome)

### Execute barcode clustering and load bc_clust data

# FOR KLEB CURVE
bc_clust <- fread("~/Desktop/ETseq_BC/Kleb_Curve_Experiment/barcode_clust/bartender/clust_d3_l4/all_bc_final_d3_l4_barcode.csv")
colnames(bc_clust) <- c("barcodes", "reads","clstID")

# FOR ET-SEQ - BIG
#bc_clust <- fread("~/Desktop/ETseq_BC/NovaSeq_MultipleTransform/barcode_clust/all_bc_d3_barcode.csv")
#colnames(bc_clust) <- c("barcodes", "reads","clstID")


## Perform some summary statistics + plots on bc_clustering

bc_clust_summary <- data.table(bc_clust %>%
                              group_by(clstID) %>%
                              summarise(SIZE = n(),
                                        RDS = sum(reads),
                                        MAX_RDS = max(reads),
                                        RST_RDS = RDS - MAX_RDS,
                                        PRIME_BC = barcodes[which.max(reads)],
                                        PRIME_FRC = MAX_RDS/RDS))

bc_general_stats <-t(data.table(Total_Clusters = nrow(bc_clust_summary),
                                Total_Variants = sum(bc_clust_summary$SIZE),
                                Total_Reads = sum(bc_clust_summary$RDS),
                                Reads_in_Prime = sum(bc_clust_summary$MAX_RDS),
                                Frac_in_Prime = sum(bc_clust_summary$MAX_RDS) / sum(bc_clust_summary$RDS),
                                Clust_5RDS = nrow(subset(bc_clust_summary, RDS >= 5)),
                                Variants_5RDS = sum(subset(bc_clust_summary, RDS >= 5)$SIZE),
                                Frac_Variants_5RDS = sum(subset(bc_clust_summary, RDS >= 5)$SIZE) / sum(bc_clust_summary$SIZE),
                                Reads_5RDS = sum(subset(bc_clust_summary, RDS >= 5)$RDS),
                                Frac_Reads_5RDS = sum(subset(bc_clust_summary, RDS >= 5)$RDS) / sum(bc_clust_summary$RDS)))


# # plot relationship between number of variants VS number of reads to most abundant variant 
# ggplot(bc_clust_summary, aes (x = SIZE, y = MAX_RDS)) +
#   geom_point(alpha = 0.2) +
#   scale_y_log10() +
#   scale_x_log10() +
#   xlab("Number of Sequence Variants") +
#   ylab("Reads to Most Abundant Variant")
# 
# 
# # plot relationship between rest of reads in cluster VS number of reads to most abundant variant
# ggplot(bc_clust_summary, aes (x = RST_RDS, y = MAX_RDS, color = SIZE)) +
#   geom_point(alpha = 0.5) +
#   scale_y_log10() +
#   scale_x_log10() 



### Assign Clusters to all_sample_bc by merging in bc_clust based on barcode cluster
all_sample_bc_clust <- merge(all_sample_bc, bc_clust[,c("barcodes","clstID")], by = "barcodes", all.x = T)

rm(all_sample_bc)

### Summarize prevelance of genomes and mapping position by barcode cluster 
master_clust_assignment <- data.table(all_sample_bc_clust %>%
                                        group_by(clstID) %>%
                                        summarise(RDS = sum(reads), # Sum of all reads assigned to BC cluster
                                                  NUM_BC = n_distinct(barcodes),
                                                  NUM_SAMP = n_distinct(SAMPLE), # Number of Samples seen in BC cluster
                                                  NUM_GEN = n_distinct(GENOME1), # Number of Genomes seen in BC cluster 
                                                  NUM_POS = n_distinct(TNjunc), # Number of unique mapping positions seen in BC cluster
                                                  PRIME_SAMP = SAMPLE[which.max(reads)], # Genome in sample, barcode, junction combination that has the most reads (This will probably work, but we are identifying based on a single combination)
                                                  PRIME_SAMP_RDS = sum(reads[which(SAMPLE == PRIME_SAMP)]), # Sum of all reads for above genome in cluster 
                                                  PRIME_SAMP_FRC = PRIME_SAMP_RDS / RDS, # Fraction that the sum of all reads for the primary genome make up of all reads in the cluster
                                                  PRIME_GEN = GENOME1[which.max(reads)],
                                                  PRIME_GEN_RDS = sum(reads[which(GENOME1 == PRIME_GEN)]),
                                                  PRIME_GEN_FRC = PRIME_GEN_RDS / RDS,
                                                  PRIME_POS = TNjunc[which.max(reads)],
                                                  PRIME_POS_RDS = sum(reads[which(TNjunc == PRIME_POS)]),
                                                  PRIME_POS_FRC = PRIME_POS_RDS / RDS,
                                                  MEAN_PURITY = (PRIME_GEN_FRC + PRIME_POS_FRC) / 2,
                                                  SWAP_RATE = ceiling(RDS*(.0004 / num_samp) + sqrt(RDS*(1-.0004)*.0004)))) # Calculates the lower bound of the number of reads that could barcode swap into any one sample based on 0.04% swap rate + 1 SD


master_clust_summary <-t(data.table(master_clust_assignment %>%
                                      summarise(Total_Clusters = n(),
                                                SAMP_n2PLUS = sum(NUM_SAMP >= 2),
                                                SAMP_n2PLUS_FRC = SAMP_n2PLUS / Total_Clusters,
                                                GEN_n2PLUS = sum(NUM_GEN >= 2),
                                                GEN_n2PLUS_FRC = GEN_n2PLUS / Total_Clusters,
                                                POS_n2PLUS = sum(NUM_POS >= 2),
                                                POS_n2PLUS_FRC = POS_n2PLUS / Total_Clusters,
                                                FRC_SAMP_PURE = sum(PRIME_SAMP_FRC >= 0.75) / Total_Clusters,
                                                FRC_GEN_PURE = sum(PRIME_GEN_FRC >= 0.75) / Total_Clusters,
                                                FRC_MEAN_PURE = sum(MEAN_PURITY >= 0.75) / Total_Clusters, 
                                                FRC_METRIC1 = sum(PRIME_GEN_RDS >= 5 & PRIME_GEN_FRC >= 0.75) / Total_Clusters,
                                                FRC_METRIC2 = sum(RDS >= 5 & MEAN_PURITY >= 0.75) / Total_Clusters)))
                                    
                                    
# rm to save mem
#rm(all_sample_bc_clust)

#bc_clust <- fread("~/Desktop/ETseq_BC/NovaSeq_MultipleTransform/barcode_clust/all_bc_d3_barcode.csv")
#colnames(bc_clust) <- c("barcodes", "reads","clstID")

filter_type <- "low"



### SELECT BC CLUSTER FILTER

if(filter_type == "high"){
  ## Strong Filter 1 ***THIS IS A PLACE WE COULD ADJUST PRIME_GEN_RDS FILTER
  master_clust_assignment_filt <- subset(master_clust_assignment, PRIME_GEN_FRC >= 0.75 & PRIME_GEN_RDS >= 5)
  bc_assign <- merge(bc_clust, master_clust_assignment_filt[,c("clstID", "PRIME_GEN","SWAP_RATE")])
}

if(filter_type == "med"){
  ## Med Filter 1
  master_clust_assignment_filt <- subset(master_clust_assignment, PRIME_GEN_FRC >= 0.75)
  bc_assign <- merge(bc_clust, master_clust_assignment_filt[,c("clstID", "PRIME_GEN","SWAP_RATE")])
}
  
if(filter_type == "low"){  
  ## Low Filter 1 - NO CLUSTERS REMOVED
  master_clust_assignment_filt <- master_clust_assignment
  bc_assign <- merge(bc_clust, master_clust_assignment_filt[,c("clstID", "PRIME_GEN","SWAP_RATE")])
}

## Check correct values merged
length(unique(bc_clust$clstID))
length(unique(master_clust_assignment_filt$clstID)) == length(unique(bc_assign$clstID))

# rm to save mem
# rm(bc_clust)
# rm(master_clust_assignment_filt)

### RUN LOOP FOR SAMPLES, ASSIGN BC CLUSTERS, AND COUNT

# Make df to hold final hit table
filt_hits <- data.table()
sample_names <- sub(".hits","",raw_files)

# loop over samples
for (i in 1:length(raw_files)){
  
  # Say what is happening
  sample <- sub(".hits","", raw_files[i])
  cat(paste0("Calculating Summary Stats for Sample: ",sample,"\n"))
  
  # Read in hit table
  tmp_hit_table <- fread(raw_files[i], sep = "\t")
  tot_reads <- nrow(tmp_hit_table)
  
  # Add Tn Junction Postion
  #tmp_hit_table$TNjunc <- ifelse(tmp_hit_table$STRAND1 == "+", tmp_filt_table$START1, tmp_filt_table$END1)
  
  
  ### APPLY FILTERING STEP 1 (each filtering step is implemented individually for now so we can add options)
  
  # Hits to same genome
  tmp_filt_table <- subset(tmp_hit_table, GENOME1 == GENOME2)
  
  # Filter by mapq
  tmp_filt_table <- subset(tmp_filt_table, MAPQ1 >= mq_cut | MAPQ2 >= mq_cut)
  
  # Filter out NA barcodes
  tmp_filt_table <- subset(tmp_filt_table, is.na(barcodes) == FALSE)
  
  # Remove tmp_hit_table to save mem
  rm(tmp_hit_table)
  
  # Calculate basic filter loss 1
  f1_loss <- nrow(tmp_filt_table) / tot_reads
  
  
  ### MERGE IN BARCODE CLUSTER NUMBERS

  # merge by bc
  tmp_filt_table <- merge(tmp_filt_table, bc_assign[,c("clstID","barcodes","PRIME_GEN","SWAP_RATE")], by = "barcodes", all.x = T)
  
  
  ### APPLY FILTERING STEP 2 - if high or med
  
  # remove reads where BC cluster was not identified
  tmp_filt_table <- subset(tmp_filt_table, is.na(clstID) == FALSE)
  
  if(filter_type %in% c("high", "med")){
    # MED + HIGH Filter 2: remove reads where GENOME != PRIME_GEN
    tmp_filt_table <- subset(tmp_filt_table, GENOME1 == PRIME_GEN)
  }
  
  # Calculate basic filter loss 2
  f2_loss <- nrow(tmp_filt_table) / tot_reads
  
  
  
  ### SUMMARIZE HIT DATA 
  
  ## Hit Summary Breakdown - Summarize hits
  
  # First aggregate number of reads by genome/model/BCcluster + include swap rate
  tmp_filt_summary <- data.table(SAMPLE = sample,
                                 tmp_filt_table %>% 
                                   group_by(GENOME1,model,clstID,SWAP_RATE) %>%
                                   summarise(counts = n()))
  
  # Now calculate all base summary statistics that include counts of reads and unique barcode clusters
  tmp_filt_summary <- data.table(SAMPLE = sample,
                                 tmp_filt_summary %>% 
                                   group_by(GENOME1,model) %>%
                                   summarise(RDS = sum(counts),
                                             UBC = n_distinct(clstID),
                                             RDSnR = sum(counts[counts >= bc_rep]),
                                             UBCnR = n_distinct(clstID[counts >= bc_rep]),
                                             RDSdC = sum(counts[counts >= 2 & counts >= SWAP_RATE]),
                                             UBCdC = n_distinct(clstID[counts >= 2 & counts >= SWAP_RATE]),
                                             BCPR = UBC / RDS))
  
  ## Hit Summary Breakdown - Calculate fractional Abundances
  tmp_filt_summary <- data.table(tmp_filt_summary,
                                 RDS_FRC = tmp_filt_summary$RDS/sum(tmp_filt_summary$RDS),
                                 UBC_FRC = tmp_filt_summary$UBC/sum(tmp_filt_summary$UBC),
                                 RDSnR_FRC = tmp_filt_summary$RDSnR/sum(tmp_filt_summary$RDSnR),
                                 UBCnR_FRC = tmp_filt_summary$UBCnR/sum(tmp_filt_summary$UBCnR),
                                 RDSdC_FRC = tmp_filt_summary$RDSdC/sum(tmp_filt_summary$RDSdC),
                                 UBCdC_FRC = tmp_filt_summary$UBCdC/sum(tmp_filt_summary$UBCdC))

  
  ## Hit Summary Breakdown - Normalize to spike in org #### POSSIBLE BUG IF SIO HITS MULTIPLE MODELS...USED MAX TO SELECT MOST ABUNDANT
  # store fractional sio data for current sample
  tmp_sio <- subset(tmp_filt_summary, GENOME1 == spike_in_org) 
  
  # normalize within sample data by dividing fractional abundances by sio fractional abundance
  tmp_filt_summary <- data.frame(tmp_filt_summary,
                                 RDS_NRM = tmp_filt_summary$RDS_FRC/max(tmp_sio$RDS_FRC),
                                 UBC_NRM = tmp_filt_summary$UBC_FRC/max(tmp_sio$UBC_FRC),
                                 RDSnR_NRM = tmp_filt_summary$RDSnR_FRC/max(tmp_sio$RDSnR_FRC),
                                 UBCnR_NRM = tmp_filt_summary$UBCnR_FRC/max(tmp_sio$UBCnR_FRC),
                                 RDSdC_NRM = tmp_filt_summary$RDSdC_FRC/max(tmp_sio$RDSdC_FRC),
                                 UBCdC_NRM = tmp_filt_summary$UBCdC_FRC/max(tmp_sio$UBCdC_FRC))
  
  
  
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


# colnames(filt_hits)[c(2,3,17)] <- c("GENOME", "MODEL", "BADBCR") 

fwrite(filt_hits, paste0("~/Desktop/ETstats_jm_hit_summary_BCclust_",filter_type,".txt"), sep = "\t")




  
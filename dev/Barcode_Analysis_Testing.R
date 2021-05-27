library(data.table)
library(dplyr)
library(ggplot2)

'%notin%' <- Negate('%in%')
#input_dir <- "~/Desktop/ETseq_BC/Kleb_Curve_Experiment/ETmapper_raw_hits/"
input_dir <- "~/Desktop/ETseq_BC/NovaSeq_MultipleTransform/ETmapper_raw_hits/"
output_dir <- "~/Desktop/ETseq_BC/NovaSeq_MultipleTransform/"
setwd(input_dir)
raw_files <- list.files(input_dir)
#i = 6
mq_cut <- 20
spike_in_org <- "Bacteroides_thetaiotaomicron_VPI_5482_GCF_000011065.1"
bc_rep <- 2

### Set up dfs to save output
#filter_losses <- data.table()
all_sample_bc <- data.table()



### Loop over files and identify all barcodes that are unique to organisims within samples,
### and make a master list of barcodes, the sample they came from, and what organism they were in
### across all samples 

for (i in 1:length(raw_files)){
  
  # Say what is happening
  cat(paste0("Calculating Summary Stats for Sample: ",raw_files[i],"\n"))
  
  # Read in hit table
  tmp_hit_table <- fread(raw_files[i], sep = "\t")
  
  
  ### APPLY FILTERING (each filtering step is implemented individually for now so we can add options)
    
    # Hits to same genome
    tmp_filt_table <- subset(tmp_hit_table, GENOME1 == GENOME2)
    
    # Filter by mapq
    tmp_filt_table <- subset(tmp_filt_table, MAPQ1 >= mq_cut | MAPQ2 >= mq_cut)
    
    # Filter out NA barcodes
    tmp_filt_table <- subset(tmp_filt_table, is.na(barcodes) == FALSE)
    
    # Calculate basic filter loss
    sample <- sub(".hits","", raw_files[i])
    f1_loss <- nrow(tmp_filt_table) / nrow(tmp_hit_table)
    
    # Add Tn Junction Postion
    tmp_filt_table$TNjunc <- ifelse(tmp_filt_table$STRAND1 == "+", tmp_filt_table$START1, tmp_filt_table$END1)
    
  ### IN SAMPLE BC FILTER: BCs in > 1 genome IN SAMPLE filtered + written to disk 
    
    # Summarize barcodes by read counts and number of genomes attached to them
    # bc2genome <- data.table(SAMPLE = sample,
    #                         tmp_filt_table %>%
    #                           group_by(barcodes) %>%
    #                           summarise(genomes = n_distinct(GENOME1),
    #                                     reads = n()))
    
    bc2genome <- data.table(SAMPLE = sample,
                            tmp_filt_table %>%
                              group_by(barcodes, GENOME1, TNjunc) %>%
                              summarise(reads = n()))

    # vector holding barcodes in > 1 genomes for filtering
    # bad_bc <- as.character(subset(bc2genome, genomes > 1)$barcodes)

    # subset and filter bad barcodes
    #tmp_bad_bc <- subset(tmp_filt_table, barcodes %in% bad_bc)
    # tmp_filt_table <- subset(tmp_filt_table, barcodes %notin% bad_bc)
    # f2_loss <- nrow(tmp_filt_table) / nrow(tmp_hit_table)
    
  ### Generate list of unique barcodes in a given sample attached to their genome names   
    
    # Generate bc2genome_name table
    # bc2genome_name <- tmp_filt_table[!duplicated(tmp_filt_table$barcodes),c("GENOME1", "barcodes")]
    
    # Get All barcodes for a sample with their genomes after filtering
    # all_filt_bc <- bc2genome[bc2genome$barcodes %in% tmp_filt_table$barcodes]
      
    # Merge in genome names
    # all_filt_bc <- merge(all_filt_bc, bc2genome_name, by = "barcodes")
    
    # bind unique barcodes per sample to all_sample_bc
    all_sample_bc <- rbind(all_sample_bc, bc2genome)
    
    # Save all barcodes from sample to disk for clustering
    
    barcode_save <- data.table(tmp_filt_table$barcodes)
    fwrite(barcode_save, paste0("~/Desktop/ETseq_BC/",sample,".bc"), sep = ",", col.names = F)
    
    
}

# rm to save mem
rm(barcode_save)


### Execute barcode clustering and load bc_clust data

bc_clust <- fread("~/Desktop/ETseq_BC/barcode_clust/bartender/clust_d3_l4/all_bc_final_d3_l4_barcode.csv")
#bc_clust <- fread("~/Desktop/ETseq_BC/barcode_clust/starcode/all_bc_final_sc_d3.out")
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


# plot relationship between number of variants VS number of reads to most abundant variant 
ggplot(bc_clust_summary, aes (x = SIZE, y = MAX_RDS)) +
  geom_point(alpha = 0.2) +
  scale_y_log10() +
  scale_x_log10() +
  xlab("Number of Sequence Variants") +
  ylab("Reads to Most Abundant Variant")


# plot relationship between rest of reads in cluster VS number of reads to most abundant variant
ggplot(bc_clust_summary, aes (x = RST_RDS, y = MAX_RDS, color = SIZE)) +
  geom_point(alpha = 0.5) +
  scale_y_log10() +
  scale_x_log10() 



### Assign Clusters to all_sample_bc by merging in bc_clust based on barcode cluster
all_sample_bc_clust <- merge(all_sample_bc, bc_clust[,c("barcodes","clstID")], by = "barcodes", all.x = T)



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
                                                  MEAN_PURITY = (PRIME_GEN_FRC + PRIME_POS_FRC) / 2))


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
                                                FRC_METRIC1 = sum(PRIME_GEN_RDS >= 3 & PRIME_GEN_FRC >= 0.75) / Total_Clusters,
                                                FRC_METRIC2 = sum(RDS >= 3 & MEAN_PURITY >= 0.75) / Total_Clusters)))
                                    
                                    


filter_type <- "low"



### SELECT BC CLUSTER FILTER

if(filter_type == "high"){
  ## Strong Filter 1 ***THIS IS A PLACE WE COULD ADJUST PRIME_GEN_RDS FILTER
  master_clust_assignment_filt <- subset(master_clust_assignment, PRIME_GEN_FRC >= 0.75 & PRIME_GEN_RDS >= 5)
  bc_assign <- merge(bc_clust, master_clust_assignment_filt[,c("clstID", "PRIME_GEN")])
}

if(filter_type == "med"){
  ## Med Filter 1
  master_clust_assignment_filt <- subset(master_clust_assignment, PRIME_GEN_FRC >= 0.75)
  bc_assign <- merge(bc_clust, master_clust_assignment_filt[,c("clstID", "PRIME_GEN")])
}
  
if(filter_type == "low"){  
  ## Low Filter 1 - NO CLUSTERS REMOVED
  master_clust_assignment_filt <- master_clust_assignment
  bc_assign <- merge(bc_clust, master_clust_assignment_filt[,c("clstID", "PRIME_GEN")])
}

## Check correct values merged
length(unique(bc_clust$clstID))
length(unique(master_clust_assignment_filt$clstID)) == length(unique(bc_assign$clstID))




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
  
  # Add Tn Junction Postion
  tmp_filt_table$TNjunc <- ifelse(tmp_filt_table$STRAND1 == "+", tmp_filt_table$START1, tmp_filt_table$END1)
  
  
  ### APPLY FILTERING STEP 1 (each filtering step is implemented individually for now so we can add options)
  
  # Hits to same genome
  tmp_filt_table <- subset(tmp_hit_table, GENOME1 == GENOME2)
  
  # Filter by mapq
  tmp_filt_table <- subset(tmp_filt_table, MAPQ1 >= mq_cut | MAPQ2 >= mq_cut)
  
  # Filter out NA barcodes
  tmp_filt_table <- subset(tmp_filt_table, is.na(barcodes) == FALSE)
  
  # Calculate basic filter loss 1
  f1_loss <- nrow(tmp_filt_table) / nrow(tmp_hit_table)
  
  
  ### MERGE IN BARCODE CLUSTER NUMBERS

  # merge by bc
  tmp_filt_table <- merge(tmp_filt_table, bc_assign[,c("clstID","barcodes","PRIME_GEN")], by = "barcodes", all.x = T)
  
  
  ### APPLY FILTERING STEP 2 - if high or med
  
  # remove reads where BC cluster was not identified
  tmp_filt_table <- subset(tmp_filt_table, is.na(clstID) == FALSE)
  
  if(filter_type %in% c("high", "med")){
    # MED + HIGH Filter 2: remove reads where GENOME != PRIME_GEN
    tmp_filt_table <- subset(tmp_filt_table, GENOME1 == PRIME_GEN)
  }
  
  # Calculate basic filter loss 2
  f2_loss <- nrow(tmp_filt_table) / nrow(tmp_hit_table)
  
  
  
  ### SUMMARIZE HIT DATA 
  
  ## Hit Summary Breakdown - Summarize hits
  tmp_filt_summary <- data.table(SAMPLE = sample,
                                 tmp_filt_table %>% 
                                   group_by(GENOME1,model) %>% ## This may need to be fixed for pe if genomes are not equal (decide which read better for GENOME call)
                                   summarise(RDS = n(),
                                             UBC = n_distinct(clstID),
                                             RDSnR = sum(plyr::count(clstID)[which(plyr::count(clstID)[,2] >= bc_rep),2]),
                                             UBCnR = sum(plyr::count(clstID)[,2] >= bc_rep),
                                             BCPR = n_distinct(clstID) / n()))
  
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
                                 RDSnR_NRM = tmp_filt_summary$RDSnR_FRC/max(tmp_sio$RDSnR_FRC),
                                 UBCnR_NRM = tmp_filt_summary$UBCnR_FRC/max(tmp_sio$UBCnR_FRC))
  
  
  
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















### Summmarize statistics about barcodes across all samples 

# Summarize barcodes by read counts and number of genomes attached to them
all_sample_bc_summary <- data.table(all_sample_bc %>%
                                      group_by(barcodes) %>%
                                      summarise(tot_samples = n_distinct(SAMPLE),
                                                tot_genomes = n_distinct(GENOME1),
                                                tot_reads = sum(reads),
                                                min_reads = min(reads),
                                                max_reads = max(reads),
                                                med_reads = median(reads)))

all_sample_bc_summary$SorG_Multi <- ifelse(all_sample_bc_summary$tot_samples > 1 | all_sample_bc_summary$tot_genomes > 1, TRUE, FALSE)

summary_of_summary <- data.table(all_sample_bc_summary %>%
                                   group_by(SorG_Multi) %>%
                                   summarise(count = n(),
                                             tot_rds = sum(tot_reads)))




# Plot stuff from above
plyr::count(all_sample_bc_summary$SorG_Multi)

ggplot(all_sample_bc_summary, aes(x = as.factor(genomes), y = log10(reads + 1))) +
  geom_boxplot()

ggplot(all_sample_bc_summary, aes(x = as.factor(samples), y = log10(reads + 1))) +
  geom_boxplot()

ggplot(all_sample_bc_summary, aes(x = SorG_Multi, y = log10(tot_reads + 1))) +
  geom_boxplot()








### Hit Summary Breakdown - Summarize hits
bc_rep <- 3

tmp_filt_summary <- data.table(SAMPLE = sample,
                               tmp_filt_table %>%
                                 group_by(GENOME1,model) %>% ## This may need to be fixed for pe if genomes are not equal (decide which read better for GENOME call)
                                 summarise(RDS = n(),
                                           UBC = n_distinct(barcodes),
                                           RDSnR = sum(plyr::count(barcodes)[which(plyr::count(barcodes)[,2] >= bc_rep),2]),
                                           UBCnR = sum(plyr::count(barcodes)[,2] >= bc_rep),
                                           BCPR = n_distinct(barcodes) / n()))







### For each sample identify the non-B. theta barcodes, then look if those barcodes overlap with barcodes
### from any other sample (including B. theta), leaving out the sample that is currently being analyzed


sample_non_sio <- subset(all_sample_bc, SAMPLE == "BR_1_14_Kcurve_14" & GENOME1 != "Bacteroides_thetaiotaomicron_VPI_5482_GCF_000011065.1")

rest <- subset(all_sample_bc, SAMPLE != "BR_1_14_Kcurve_14")

in_sample_overlap <- sample_non_sio[sample_non_sio$barcodes %in% rest$barcodes,]

plyr::count(in_sample_overlap$GENOME1)

out_sample_overlap <- rest[rest$barcodes %in% sample_non_sio$barcodes,]

plyr::count(out_sample_overlap$SAMPLE)

plyr::count(out_sample_overlap$GENOME1)

foo <- merge(out_sample_overlap, in_sample_overlap, by = "barcodes", all = T)

nrow(subset(foo, GENOME1.x == GENOME1.y))





### MANUAL SHIT - Barcode Clustering

foo <- subset(tmp_filt_table, GENOME1 == "Klebsiella_michiganensis_M5a1_GCF_002090195.1")
foo <- subset(tmp_filt_table, barcodes == "GTCAACCTATGGGGATTTGT")

bar <- plyr::count(foo$barcodes)
bar <- bar[order(bar$x),]

bc_dist <- adist(foo$barcodes)
plot(bc_dist, cex = 0.5)

bc_cut <- cutree(bc_dist, h = fe+1)
bc_flank <- bc_flank[!duplicated(bc_cut)]

nex <- adist(bar$x)


### asdadkj
## Are reads actually chimeras???
## Are barcodes sample swapped...only for non-b.theta but compare to samples with b.theta
## Manually inspect barcodes to look for may bc that are super similar in sequence...are we seeing something like polymerase slippage that gives a lot of BC in one sample
## STEP1: 
  
  
  
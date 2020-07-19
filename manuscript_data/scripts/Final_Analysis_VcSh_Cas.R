library(data.table)
library(dplyr)
library(ggplot2)


### NOTE ON GENOMIC TARGETS
## DUPLICATION REGION IN E. COLI GENOME: 749,903 bp --> 750,380 bp
## VcCas Protospacer: 5'-CCCTTTCGCCAGCTGGCGTAATAGCGAAGAGG-3' (- Strand)
## VcCas Protospacer Coordinates: 335,714 bp <-- 335,745 bp
## ShCas Protospacer: 


### Create Final VcCas Container
final_VcCas_dat <- data.table()

### Read In Data

## Barcode References 
#all_bc_ref <- fread("~/Dropbox/Banfield_Lab_Files/Projects/mCAFE/mCAFE_Project/Studies/19_10_21_JD_BR_VCas/ET_Seq/hc/all_bc_cluster_master_assignment_filt.txt")
tmp_hit_table <- fread("~/Dropbox/Banfield_Lab_Files/Projects/mCAFE/mCAFE_Project/Studies/19_10_21_JD_BR_VCas/ET_Seq/hc/VcCas_rep3.hits")


### Perform Basic MapQ Filtering
#tmp_hit_table <- subset(tmp_hit_table, MAPQ1 >= 3 & MAPQ2 >= 3)


### Aggregate number of reads by genome/model/BCcluster + include swap rate
tmp_hit_summary <- data.table(tmp_hit_table %>% 
                                group_by(GENOME,model,clstID,SWAP_RATE) %>%
                                summarise(counts = n()))

# Identify barcode clusters that pass filtering criteria 
tmp_hit_summary <- data.frame(tmp_hit_summary,
                              PASS_FILTER = ifelse(tmp_hit_summary$counts > 2 & tmp_hit_summary$counts > tmp_hit_summary$SWAP_RATE, TRUE, FALSE))

#tmp_hit_summary <- data.frame(tmp_hit_summary,
#                              PASS_FILTER = ifelse(tmp_hit_summary$counts > 2, TRUE, FALSE))

good_clust <- subset(tmp_hit_summary, PASS_FILTER == TRUE)$clstID


# Filter Reads for only those associated with good clusters 
tmp_hit_filt <- subset(tmp_hit_table, clstID %in% good_clust)




### Identify Confidenct Unique Hit Locations (Clusters) and Strandedness 

### Summarize Reads comming from Identical Positions and Directions 
TNjunc_summary <- data.table(tmp_hit_filt %>%
                                group_by(clstID, TNjunc, STRAND1) %>%
                                summarise(reads = n()))

### Identify good mapping positions based on clusters that overwelmingly go to a certian position
TNjunc_stats <- data.table(TNjunc_summary %>% 
                            group_by(clstID) %>%
                            summarise(PRIME_POS = TNjunc[which.max(reads)],
                                      PRIME_POS_RDS = sum(reads[which(TNjunc == PRIME_POS)]),
                                      PRIME_POS_FRC = PRIME_POS_RDS / sum(reads),
                                      PRIME_DIR = STRAND1[which.max(reads)],
                                      PRIME_DIR_RDS = sum(reads[which(STRAND1 == PRIME_DIR)]),
                                      PRIME_DIR_FRC = PRIME_DIR_RDS / sum(reads),
                                      PRIME_POS3_RDS = sum(reads[which(TNjunc %in% (PRIME_POS-3):(PRIME_POS+3))]),
                                      PRIME_POS3_FRC = PRIME_POS3_RDS / sum(reads)))


### Filter Out Unique Hits that dont have <= 75% of positions (within 3 bp) or directionality matching
TNjunc_Final <- subset(TNjunc_stats, PRIME_DIR_FRC >= 0.75 & PRIME_POS3_FRC >= 0.75)
#TNjunc_Final <- TNjunc_stats

### Create Final Hit Table

sample <- "VcCas_3"
final_tmp <- data.table(SAMPLE = sample,
                        clstID = TNjunc_Final$clstID,
                        TNjunc = TNjunc_Final$PRIME_POS,
                        R1_DIR = TNjunc_Final$PRIME_DIR)


### Aggregate Samples 
final_VcCas_dat <- rbind(final_VcCas_dat, final_tmp)





### Add In metadata
final_VcCas_dat$ON_TARG <- ifelse((final_VcCas_dat$TNjunc >= 335514 & final_VcCas_dat$TNjunc <= 335713), TRUE, FALSE)
final_VcCas_dat$ORI <- ifelse((final_VcCas_dat$TNjunc >= 335514 & final_VcCas_dat$TNjunc <= 335713 & final_VcCas_dat$R1_DIR == "-"), "SAME", "OPPO")
#final_VcCas_dat$ON_TARG <- ifelse((final_VcCas_dat$TNjunc >= 335514 & final_VcCas_dat$TNjunc <= 335713) | (final_VcCas_dat$TNjunc >= 750068 & final_VcCas_dat$TNjunc <= 750267), TRUE, FALSE)
#final_VcCas_dat$ORI <- ifelse((final_VcCas_dat$TNjunc >= 335514 & final_VcCas_dat$TNjunc <= 335713 & final_VcCas_dat$R1_DIR == "-") | (final_VcCas_dat$TNjunc >= 750068 & final_VcCas_dat$TNjunc <= 750267 & final_VcCas_dat$R1_DIR == "+"), "SAME", "OPPO")


## 200bp Target Window
# VC CAS WINDOW 1 --> 335,514 - 335,713
# VC CAS WINDOW 2 --> 750,068 - 750,267

# VcCas bulk data testing  Test
#foo <- subset(TNjunc_Final, (PRIME_POS >= 335514 & PRIME_POS <= 335713) | (PRIME_POS >= 750068 & PRIME_POS <= 750267))
#foo <- subset(tmp_hit_filt, (TNjunc >= 335514 & TNjunc <= 335713) | (TNjunc >= 750068 & TNjunc <= 750267))








#### ShCas Analysis 



### Create Final ShCas Container
final_ShCas_dat <- data.table()

### Read In Data

## Barcode References 
#all_bc_ref <- fread("~/Dropbox/Banfield_Lab_Files/Projects/mCAFE/mCAFE_Project/Studies/19_10_21_JD_BR_VCas/ET_Seq/hc/all_bc_cluster_master_assignment_filt.txt")
tmp_hit_table <- fread("~/Dropbox/Banfield_Lab_Files/Projects/mCAFE/mCAFE_Project/Studies/19_10_21_JD_BR_VCas/ET_Seq/hc/ShCas_rep3.hits")


### Perform Basic MapQ Filtering
#tmp_hit_table <- subset(tmp_hit_table, MAPQ1 >= 3 & MAPQ2 >= 3)


### Aggregate number of reads by genome/model/BCcluster + include swap rate
tmp_hit_summary <- data.table(tmp_hit_table %>% 
                                group_by(GENOME,model,clstID,SWAP_RATE) %>%
                                summarise(counts = n()))

# Identify barcode clusters that pass filtering criteria 
tmp_hit_summary <- data.frame(tmp_hit_summary,
                              PASS_FILTER = ifelse(tmp_hit_summary$counts >= 2 & tmp_hit_summary$counts >= tmp_hit_summary$SWAP_RATE, TRUE, FALSE))

#tmp_hit_summary <- data.frame(tmp_hit_summary,
#                              PASS_FILTER = ifelse(tmp_hit_summary$counts > 5, TRUE, FALSE))

good_clust <- subset(tmp_hit_summary, PASS_FILTER == TRUE)$clstID


# Filter Reads for only those associated with good clusters 
tmp_hit_filt <- subset(tmp_hit_table, clstID %in% good_clust)




### Identify Confidenct Unique Hit Locations (Clusters) and Strandedness 

### Summarize Reads comming from Identical Positions and Directions 
TNjunc_summary <- data.table(tmp_hit_filt %>%
                               group_by(clstID, TNjunc, STRAND1) %>%
                               summarise(reads = n()))

### Identify good mapping positions based on clusters that overwelmingly go to a certian position
TNjunc_stats <- data.table(TNjunc_summary %>% 
                             group_by(clstID) %>%
                             summarise(PRIME_POS = TNjunc[which.max(reads)],
                                       PRIME_POS_RDS = sum(reads[which(TNjunc == PRIME_POS)]),
                                       PRIME_POS_FRC = PRIME_POS_RDS / sum(reads),
                                       PRIME_DIR = STRAND1[which.max(reads)],
                                       PRIME_DIR_RDS = sum(reads[which(STRAND1 == PRIME_DIR)]),
                                       PRIME_DIR_FRC = PRIME_DIR_RDS / sum(reads),
                                       PRIME_POS3_RDS = sum(reads[which(TNjunc %in% (PRIME_POS-3):(PRIME_POS+3))]),
                                       PRIME_POS3_FRC = PRIME_POS3_RDS / sum(reads)))


### Filter Out Unique Hits that dont have <= 75% of positions (within 3 bp) or directionality matching
TNjunc_Final <- subset(TNjunc_stats, PRIME_DIR_FRC >= 0.75 & PRIME_POS3_FRC >= 0.75)
#TNjunc_Final <- TNjunc_stats

### Create Final Hit Table
sample <- "ShCas_3"
final_tmp <- data.table(SAMPLE = sample,
                        clstID = TNjunc_Final$clstID,
                        TNjunc = TNjunc_Final$PRIME_POS,
                        R1_DIR = TNjunc_Final$PRIME_DIR)


### Aggregate Samples 
final_ShCas_dat <- rbind(final_ShCas_dat, final_tmp)






### Add In metadata
final_ShCas_dat$ON_TARG <- ifelse((final_ShCas_dat$TNjunc >= 335585 & final_ShCas_dat$TNjunc <= 335784), TRUE, FALSE)
final_ShCas_dat$ORI <- ifelse((final_ShCas_dat$TNjunc >= 335585 & final_ShCas_dat$TNjunc <= 335784 & final_ShCas_dat$R1_DIR == "-"), "SAME", "OPPO")
#final_ShCas_dat$ON_TARG <- ifelse((final_ShCas_dat$TNjunc >= 335585 & final_ShCas_dat$TNjunc <= 335784) | (final_ShCas_dat$TNjunc >= 749997 & final_ShCas_dat$TNjunc <= 750196), TRUE, FALSE)
#final_ShCas_dat$ORI <- ifelse((final_ShCas_dat$TNjunc >= 335585 & final_ShCas_dat$TNjunc <= 335784 & final_ShCas_dat$R1_DIR == "-") | (final_ShCas_dat$TNjunc >= 749997 & final_ShCas_dat$TNjunc <= 750196 & final_ShCas_dat$R1_DIR == "+"), "SAME", "OPPO")





## 200bp Target Window
# SH CAS WINDOW 1 --> 335,585 - 335,784
# SH CAS WINDOW 2 --> 749,997 - 750,169
# ShCas Test
#foo <- subset(TNjunc_Final, (PRIME_POS >= 335585 & PRIME_POS <= 335784) | (PRIME_POS >= 749997 & PRIME_POS <= 750196))
#foo <- subset(tmp_hit_table, (TNjunc >= 335585 & TNjunc <= 335784) | (TNjunc >= 749997 & TNjunc <= 750196))



### Plotting Genome Wide Inserts - VcCas

## Aggregate hit position data
Vc_hit_pos_agg <- plyr::count(final_VcCas_dat$TNjunc)
colnames(Vc_hit_pos_agg) <- c("bp", "UBC")
Vc_hit_pos_agg$Mb <- Vc_hit_pos_agg$bp/1000000
Vc_hit_pos_agg$dist <- abs(335714 - Vc_hit_pos_agg$bp)

## Full Genome Insertion Plot 
ggplot(Vc_hit_pos_agg, aes(x = Mb, y = log10(UBC + 1))) +
  geom_bar(stat = "identity", fill = "steelblue", color = "steelblue") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 4.559)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 4)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black"),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black")) +
  ylab("Unique Insertions (log10)") +
  xlab("Genomic Position (Mb)")

#ggsave("~/Desktop/ET_Seq_Figs/Pass3/Fig4B1_v3.pdf",device = "pdf", width = 178, height = 60, units = "mm", dpi = 300)


## Inset Plot - with coords 
ggplot(Vc_hit_pos_agg, aes(x = bp, y = log10(UBC + 1))) +
  geom_bar(stat = "identity", fill = "steelblue", color = "steelblue") +
  geom_rect(aes(xmin = 335714, xmax = 335745, ymin = 0, ymax = 0.2)) +
  scale_x_continuous(expand = c(0, 0), limits = c(335514, 335760)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 4)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black", size = 5),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black")) +
  ylab(NULL) +
  xlab(NULL)


## Inset Plot - with dist 
ggplot(Vc_hit_pos_agg, aes(x = dist, y = log10(UBC + 1))) +
  geom_bar(stat = "identity", fill = "steelblue", color = "steelblue") +
  #geom_rect(aes(xmin = 335714, xmax = 335745, ymin = 0, ymax = 0.2)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 60)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 4)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black", size = 5),
        axis.title = element_text(color = "black", size = 7),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black")) +
  ylab("Unique Insertions (log10)") +
  xlab("Distance (bp)")

#ggsave("~/Desktop/ET_Seq_Figs/Pass3/Fig4B1_inset_v3.pdf",device = "pdf", width = 60, height = 40, units = "mm", dpi = 300)






### Plotting Genome Wide Inserts - ShCas

## Aggregate hit position data
Sh_hit_pos_agg <- plyr::count(final_ShCas_dat$TNjunc)
colnames(Sh_hit_pos_agg) <- c("bp", "UBC")
Sh_hit_pos_agg$Mb <- Sh_hit_pos_agg$bp/1000000
Sh_hit_pos_agg$dist <- abs(335785 - Sh_hit_pos_agg$bp)


ggplot(Sh_hit_pos_agg, aes(x = Mb, y = log10(UBC + 1))) +
  geom_bar(stat = "identity", fill = "steelblue", color = "steelblue") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 4.559)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 4)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black"),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black")) +
  ylab("Unique Insertions (log10)") +
  xlab("Genomic Position (Mb)")

#ggsave("~/Desktop/ET_Seq_Figs/Pass3/Fig4B2_v4.pdf",device = "pdf", width = 178, height = 60, units = "mm", dpi = 300)
  

## Inset Plot - with coords 
ggplot(Sh_hit_pos_agg, aes(x = bp, log10(UBC + 1))) +
  geom_bar(stat = "identity", fill = "steelblue", color = "steelblue") +
  geom_rect(aes(xmin = 335785, xmax = 335807, ymin = 0, ymax = 0.2)) +
  scale_x_continuous(expand = c(0, 0), limits = c(335585, 335812)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 4)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black", size = 5),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black")) +
  ylab(NULL) +
  xlab(NULL)


## Inset Plot - with dist 
ggplot(Sh_hit_pos_agg, aes(x = dist, y = log10(UBC + 1))) +
  geom_bar(stat = "identity", fill = "steelblue", color = "steelblue") +
  #geom_rect(aes(xmin = 335785, xmax = 335807, ymin = 0, ymax = 0.2)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 60)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 4)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black", size = 5),
        axis.title = element_text(color = "black", size = 7),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black")) +
  ylab("Unique Insertions (log10)") +
  xlab("Distance (bp)")


#ggsave("~/Desktop/ET_Seq_Figs/Pass3/Fig4B2_inset_v3.pdf",device = "pdf", width = 60, height = 40, units = "mm", dpi = 300)



### Summary Statistics

## VcCas
VcCas_summary_stats <- data.table(final_VcCas_dat %>%
                                    group_by(SAMPLE) %>%
                                    summarise(TYPE = "Vc_Cas",
                                              TOTAL = n(),
                                              ON = sum(ON_TARG),
                                              OFF = TOTAL - ON,
                                              RATE = ON/TOTAL))

## ShCas
ShCas_summary_stats <- data.table(final_ShCas_dat %>%
                                    group_by(SAMPLE) %>%
                                    summarise(TYPE = "Sh_Cas",
                                              TOTAL = n(),
                                              ON = sum(ON_TARG),
                                              OFF = TOTAL - ON,
                                              RATE = ON/TOTAL))



## Aggregate Summaries 
All_summary_stats <- rbind(VcCas_summary_stats, ShCas_summary_stats)

All_summary_stats_errorbar <- data.table(All_summary_stats %>%
                                           group_by(TYPE) %>%
                                           summarise(MEAN_RATE = mean(RATE),
                                                     SD_RATE = sd(RATE)))


All_summary_stats$TYPE <- factor(All_summary_stats$TYPE, levels=c("Vc_Cas", "Sh_Cas"))
All_summary_stats_errorbar$TYPE <- factor(All_summary_stats_errorbar$TYPE, levels=c("Vc_Cas", "Sh_Cas"))

## Plot Summary of on target efficency 
ggplot(All_summary_stats, aes(x = reorder(TYPE, RATE), y = RATE, fill = TYPE)) +
  geom_linerange(data = All_summary_stats_errorbar, aes(x = TYPE, y = SD_RATE, ymin = SD_RATE, ymax = SD_RATE), color = "black") +
  geom_crossbar(data = All_summary_stats_errorbar, aes(x = TYPE, y = MEAN_RATE, ymin = MEAN_RATE, ymax = MEAN_RATE), width = 0.5, color = "black") +
  geom_point(position = position_jitter(width = 0.1), size = 3, shape = 21, alpha = 0.8) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black"),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_blank(),
        legend.position = "none") +
  scale_fill_brewer(palette = "Set1") +
  xlab(NULL) +
  ylab("On Target (200 bp window)")


ggsave("~/Desktop/ET_Seq_Figs/Pass3/Fig4C_v1.svg",device = "svg", width = 60, height = 70, units = "mm", dpi = 300)





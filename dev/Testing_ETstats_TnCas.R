library(data.table)
library(plyr)
library(ggplot2)

### VcCas

# Vc Protospacer = 5'-[CC]-ccctttcgccagctggcgtaatagcgaagagg-3' (PAM in brackets)
# Vc Target site 1 (intended) = 335741-335708 (- strand)
# Vc Target site 2 (in DE3 region) = 750028-750061 (+ strand)


hit_dat <- fread("../Studies/19_10_21_JD_BR_VCas/TnCas_Inserts/VcCas/JD_BR_VcSh_6G.hits", header = T, stringsAsFactors = F)



### Subset and filter
e_coli_hits <- subset(hit_dat, GENOME1 == "Escherichia_coli_BL21_CP001509.3")
hits_filt <-  subset(e_coli_hits, GENOME1 == GENOME2 & MAPQ1 > 0 & MAPQ2 > 0 & NM1 < 5 & NM2 < 5)



### Parse mapping locations and account for stranded-ness 
hits_pos_strand <- subset(hits_filt, STRAND1 == "+")
loci_pos_strand <- count(hits_pos_strand$START1)
loci_pos_strand$dir <- "+"


hits_neg_strand <- subset(hits_filt, STRAND1 == "-")
loci_neg_strand <- count(hits_neg_strand$END1)
loci_neg_strand$dir <- "-"

loci <- rbind(loci_pos_strand, loci_neg_strand)

loci <- data.frame(coord = loci$x, reads = loci$freq, strand = loci$dir, uniq_bc = NA)



### Identify On Target Locations - within 100bp in front and 200bp in front
loci$on_targ_100 <- ifelse(between(loci$coord,335607,335707) | between(loci$coord,750062,750162), TRUE, FALSE) 
loci$on_targ_200 <- ifelse(between(loci$coord,335607,335807) | between(loci$coord,750062,750262), TRUE, FALSE)


### Count unique barcodes at correct start site for strand & save hits to df 

# Create container to store all correct hits 
correct_hits <- data.table()

# Loop over loci and store good hits with correct start site and count unique barcodes for summary data
for (i in 1:nrow(loci)){
  
  if(loci$strand[i] == "+"){
    tmp_loci <- subset(hits_pos_strand, START1 == loci[i,1])
    tmp_loci$START <- tmp_loci$START1
    correct_hits <- rbind(correct_hits, tmp_loci)
    tmp_bc_count <- length(unique(tmp_loci$barcodes))
    loci[i,4] <- tmp_bc_count
  }
  
  if(loci$strand[i] =="-"){
    tmp_loci <- subset(hits_neg_strand, END1 == loci[i,1])
    tmp_loci$START <- tmp_loci$END1
    correct_hits <- rbind(correct_hits, tmp_loci)
    tmp_bc_count <- length(unique(tmp_loci$barcodes))
    loci[i,4] <- tmp_bc_count
  }
  
  cat(paste0("Processed ",i," Loci\r"))
}



### Calculate Summary Statistics
summary_100 <- aggregate(uniq_bc ~ on_targ_100, loci, sum)
summary_200 <- aggregate(uniq_bc ~ on_targ_200, loci, sum)

# reads vs unique_bc per position
# ggplot(loci, aes(x = log(reads + 1), y = log(uniq_bc + 1))) +
#   geom_point(alpha = 0.5) +
#   stat_smooth()


# Full plot
ggplot(loci, aes(x = coord, y = uniq_bc, color = strand, fill = strand)) +
  geom_bar(stat = "identity", position = "stack") +
  xlab("Genome Position (bp)") +
  scale_y_sqrt() +
  ylab("Unique Barcodes") 

ggsave("../Studies/19_10_21_JD_BR_VCas/TnCas_Inserts/Output_Figs/6G_hit_loci_global.pdf")


# Target Locus 1 335741-335708 (- strand)
ggplot(loci, aes(x = coord, y = uniq_bc, fill = strand)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_vline(aes(xintercept=335608), color = "black", linetype = 2) +
  geom_vline(aes(xintercept=335707), color = "black", linetype = 2) +
  geom_rect(aes(xmin=335708, xmax=335741, ymin=1, ymax=5), fill="firebrick", color = "black", inherit.aes = FALSE) +
  xlim(c(335500, 335850)) +
  scale_y_sqrt() +
  xlab("Genome Position (bp)") +
  ylab("Unique Barcodes")

ggsave("../Studies/19_10_21_JD_BR_VCas/TnCas_Inserts/Output_Figs/6G_hit_loci_site1.pdf")

# Target Locus 1 Density Plot
foo <- correct_hits[between(correct_hits$START1,335507, 335707),] # plot 200bp in front only
bar <- foo[!duplicated(foo$barcodes),]

ggplot(bar, aes(x = START)) +
  xlim(c(335500, 335850)) +
  geom_histogram(aes(y = ..density..),fill="grey60", binwidth = 1) +
  stat_density(color = "blue", fill = NA, bw = 4, linetype = 2) +
  geom_hline(yintercept=0, colour="white", size=1) +
  geom_rect(aes(xmin=335708, xmax=335741, ymin=0, ymax=0.001), fill="firebrick", color = "black", inherit.aes = FALSE) +
  scale_y_sqrt() +
  xlab("Genome Position (bp)") +
  ylab("Density Unique Barcodes")

ggsave("../Studies/19_10_21_JD_BR_VCas/TnCas_Inserts/Output_Figs/6G_hit_loci_site1_density.pdf")


# Target Locus 2 750028-750061 (+ strand)
ggplot(loci, aes(x = coord, y = uniq_bc, fill = strand)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_vline(aes(xintercept=750062), color = "black", linetype = 2) +
  geom_vline(aes(xintercept=750162), color = "black", linetype = 2) +
  geom_rect(aes(xmin=750028, xmax=750061, ymin=1, ymax=5), fill = "firebrick", color="black", inherit.aes = FALSE) +
  xlim(c(749928, 750261)) +
  scale_y_sqrt() +
  xlab("Genome Position (bp)") +
  ylab("Unique Barcodes")

ggsave("../Studies/19_10_21_JD_BR_VCas/TnCas_Inserts/Output_Figs/6G_hit_loci_site2.pdf")

# Target Locus 2 Density Plot
foo <- correct_hits[between(correct_hits$START1,750062, 750262),] # plot 200bp in front only
bar <- foo[!duplicated(foo$barcodes),]

ggplot(bar, aes(x = START)) +
  xlim(c(749928, 750261)) +
  geom_histogram(aes(y = ..density..),fill="grey60", binwidth = 1) +
  stat_density(color = "blue", fill = NA, bw = 4, linetype = 2) +
  geom_hline(yintercept=0, colour="white", size=1) +
  geom_rect(aes(xmin=750028, xmax=750061, ymin=0, ymax=0.001), fill="firebrick", color = "black", inherit.aes = FALSE) +
  scale_y_sqrt() +
  xlab("Genome Position (bp)") +
  ylab("Density Unique Barcodes")

ggsave("../Studies/19_10_21_JD_BR_VCas/TnCas_Inserts/Output_Figs/6G_hit_loci_site2_density.pdf")











### ShCas

# Sh Protospacer = 5'-[GTT]-ttacaacgtcgtgactgggaaaa-3' (PAM in brackets)
# Sh Target site 1 (intended) = 335801-335779 (- strand)
# Sh Target site 2 (in DE3 region) = 749968-749990 (+ strand) (edited)
# Suposed Insert Distance from Proto ~ 60-66bp


bc_dat <- fread("../Studies/19_10_21_JD_BR_VCas/TnCas_Inserts/ShCas/JD_BR_VcSh_15G.bc", header = T, stringsAsFactors = F)
hit_dat <- fread("../Studies/19_10_21_JD_BR_VCas/TnCas_Inserts/ShCas/JD_BR_VcSh_15G.hits", header = T, stringsAsFactors = F)


bc_dat$read <- sub(" 1:N:0:.*$","", bc_dat$read, perl = T)

merge_dat <- merge(hit_dat, bc_dat, by.x = "Read", by.y ="read", all.x = T)

e_coli_hits <- subset(merge_dat, GENOME == "Escherichia_coli_BL21_CP001509.3")

hits_filt <- subset(e_coli_hits, MAPQ > 0 & NM < 5 & LEN >= 50)



### Parse mapping locations and account for stranded-ness 
hits_pos_strand <- subset(hits_filt, STRAND == "+")
loci_pos_strand <- count(hits_pos_strand$START)
loci_pos_strand$dir <- "+"


hits_neg_strand <- subset(hits_filt, STRAND == "-")
loci_neg_strand <- count(hits_neg_strand$END)
loci_neg_strand$dir <- "-"

loci <- rbind(loci_pos_strand, loci_neg_strand)

loci <- data.frame(coord = loci$x, reads = loci$freq, strand = loci$dir, uniq_bc = NA)



### Identify On Target Locations - within 100bp in front and 200bp in front
loci$on_targ_100 <- ifelse(between(loci$coord,335778,335678) | between(loci$coord,749990,750090), TRUE, FALSE)  # added +-1 to last target bp end
loci$on_targ_200 <- ifelse(between(loci$coord,335778,335578) | between(loci$coord,749990,750190), TRUE, FALSE)  # added +-1 to last target bp end



### Count unique barcodes at correct start site for strand & save hits to df 

# Create container to store all correct hits 
correct_hits <- data.table()

# Loop over loci and store good hits with correct start site and count unique barcodes for summary data
for (i in 1:nrow(loci)){
  
  if(loci$strand[i] == "+"){
    tmp_loci <- subset(hits_pos_strand, START == loci[i,1])
    tmp_loci$STARTa <- tmp_loci$START
    correct_hits <- rbind(correct_hits, tmp_loci)
    tmp_bc_count <- length(unique(tmp_loci$barcodes))
    loci[i,4] <- tmp_bc_count
  }
  
  if(loci$strand[i] =="-"){
    tmp_loci <- subset(hits_neg_strand, END == loci[i,1])
    tmp_loci$STARTa <- tmp_loci$END
    correct_hits <- rbind(correct_hits, tmp_loci)
    tmp_bc_count <- length(unique(tmp_loci$barcodes))
    loci[i,4] <- tmp_bc_count
  }
  
  cat(paste0("Processed ",i," Loci\r"))
}



### Calculate Summary Statistics
summary_100 <- aggregate(uniq_bc ~ on_targ_100, loci, sum)
summary_200 <- aggregate(uniq_bc ~ on_targ_200, loci, sum)



### Produce Output Plots 

# Full plot
ggplot(loci, aes(x = coord, y = uniq_bc, color = strand, fill = strand)) +
  geom_bar(stat = "identity", position = "stack") +
  xlab("Genome Position (bp)") +
  ylab("Unique Barcodes") 

ggsave("../Studies/19_10_21_JD_BR_VCas/TnCas_Inserts/Output_Figs/15G_hit_loci_global.pdf")

# Target Locus 1 335801-335779 (- strand)
ggplot(loci, aes(x = coord, y = uniq_bc, fill = strand)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_vline(aes(xintercept=335778), color = "black", linetype = 2) +
  geom_vline(aes(xintercept=335678), color = "black", linetype = 2) +
  geom_rect(aes(xmin=335779, xmax=335801, ymin=1, ymax=1.5), fill="firebrick", color = "black", inherit.aes = FALSE) +
  xlim(c(335500, 335850)) +
  #scale_y_sqrt() +
  xlab("Genome Position (bp)") +
  ylab("Unique Barcodes")

ggsave("../Studies/19_10_21_JD_BR_VCas/TnCas_Inserts/Output_Figs/15G_hit_loci_site1.pdf")

# Target Locus 1 Density Plot
foo <- correct_hits[between(correct_hits$STARTa,335578, 335778),] # plot 200bp in front only
bar <- foo[!duplicated(foo$barcodes),]

ggplot(bar, aes(x = STARTa)) +
  geom_histogram(aes(y = ..density..),fill="grey60", binwidth = 1) +
  stat_density(color = "blue", fill = NA, bw = 6, linetype = 2) +
  geom_hline(yintercept=0, colour="white", size=1) +
  geom_rect(aes(xmin=335779, xmax=335801, ymin=0, ymax=0.001), fill="firebrick", color = "black", inherit.aes = FALSE) +
  xlim(c(335500, 335850)) +
  #scale_y_sqrt() +
  xlab("Genome Position (bp)") +
  ylab("Density Unique Barcodes")

ggsave("../Studies/19_10_21_JD_BR_VCas/TnCas_Inserts/Output_Figs/15G_hit_loci_site1_density.pdf")


# Target Locus 2 749968-749990 (+ strand)
ggplot(loci, aes(x = coord, y = uniq_bc, fill = strand)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_vline(aes(xintercept=749991), color = "black", linetype = 2) +
  geom_vline(aes(xintercept=750091), color = "black", linetype = 2) +
  geom_rect(aes(xmin=749968, xmax=749990, ymin=1, ymax=1.5), fill="firebrick", color = "black", inherit.aes = FALSE) +
  xlim(c(749928, 750261)) +
  #scale_y_sqrt() +
  xlab("Genome Position (bp)") +
  ylab("Unique Barcodes")

ggsave("../Studies/19_10_21_JD_BR_VCas/TnCas_Inserts/Output_Figs/15G_hit_loci_site2.pdf")

# Target Locus 1 Density Plot
foo <- correct_hits[between(correct_hits$STARTa,749991, 750191),] # plot 200bp in front only
bar <- foo[!duplicated(foo$barcodes),]

ggplot(bar, aes(x = STARTa)) +
  geom_histogram(aes(y = ..density..),fill="grey60", binwidth = 1) +
  stat_density(color = "blue", fill = NA, bw = 6, linetype = 2) +
  geom_hline(yintercept=0, colour="white", size=1) +
  geom_rect(aes(xmin=749968, xmax=749990, ymin=0, ymax=0.001), fill="firebrick", color = "black", inherit.aes = FALSE) +
  xlim(c(749928, 750261)) +
  #scale_y_sqrt() +
  xlab("Genome Position (bp)") +
  ylab("Density Unique Barcodes")

ggsave("../Studies/19_10_21_JD_BR_VCas/TnCas_Inserts/Output_Figs/15G_hit_loci_site2_density.pdf")











### VcCas Fwd Read Only
# bc_dat <- fread("../Studies/19_10_21_JD_BR_VCas/VcCas/JD_BR_VcSh_4G.bc", header = T, stringsAsFactors = F)
# hit_dat <- fread("../Studies/19_10_21_JD_BR_VCas/test_mapping/4G_fwdonly_maptest.hits", header = F, stringsAsFactors = F)
# colnames(hit_dat) <- c("genome", "read", "MAPQ1", "NM1", "length", "START1")
# 
# bc_dat$read <- sub(" 1:N:0:ATCACG","", bc_dat$read)
# 
# merge_dat <- merge(hit_dat, bc_dat, by.x = "read", by.y ="read", all.x = T)
# 
# 
# e_coli_hits <- subset(merge_dat, genome == "Escherichia_coli_BL21_CP001509.3")
# 
# e_coli_hits_filt <- subset(e_coli_hits, MAPQ1 > 0 & NM1 < 5)
# 
# loci <- count(e_coli_hits_filt$START1)
# 
# loci <- data.frame(pos = loci$x, reads = loci$freq, uniq_bc = NA)
# 
# loci$on_targ <- ifelse(between(loci$pos,335441,336008) | between(loci$pos,749728,750361), TRUE, FALSE)
# 
# 
# for (i in 1:nrow(loci)){
#   
#   tmp_loci <- subset(e_coli_hits_filt, START1 == loci[i,1])
#   tmp_bc_count <- length(unique(tmp_loci$barcodes))
#   
#   loci[i,3] <- tmp_bc_count
#   
# }
# 
# # reads vs unique_bc per position
# ggplot(loci, aes(x = log(reads + 1), y = log(uniq_bc + 1))) +
#   geom_point(alpha = 0.5) +
#   stat_smooth()
# 
# 
# # Full plot
# ggplot(loci, aes(x = pos, y = uniq_bc)) +
#   geom_bar(stat = "identity", color = "blue") +
#   xlab("Genome Position (bp)") +
#   ylab("Unique Barcodes")
# 
# 
# # Target Locus 1 (335741-335708 +- 300bp)
# ggplot(loci, aes(x = pos, y = uniq_bc)) +
#   geom_bar(stat = "identity", color = "blue") +
#   xlim(c(335441, 336008))
# 
# 
# # Target Locus 2 (750028-750061 +- 300bp)
# ggplot(loci, aes(x = pos, y = uniq_bc)) +
#   geom_bar(stat = "identity", color = "blue") +
#   xlim(c(749728, 750361))
# 
# 
# 
# 
# ### Analysis of on (within 300bp of target) vs off target hits
# #all_samples <- list()
# #all_samples[[3]] <- loci
# 
# summary_dat <- data.frame()
# sample_names <- c("4G", "5G", "6G")
# 
# for (i in 1:length(all_samples)){
#   tmp <- all_samples[[i]]
#   tmp_ag <- aggregate(uniq_bc ~ on_targ, tmp, sum)
#   
#   tmp_ag$sample <- sample_names[i]
#   
#   summary_dat <- rbind(tmp_ag, summary_dat)
#   
# }
# 
# ggplot(summary_dat, aes(x = on_targ, y = uniq_bc)) +
#   geom_boxplot(fill = "steelblue") +
#   scale_y_log10() +
#   xlab(NULL) +
#   ylab("Unique Insertion Barcodes")

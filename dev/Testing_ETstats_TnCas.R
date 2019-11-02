library(data.table)
library(plyr)
library(ggplot2)

### VcCas

bc_dat <- fread("../Studies/19_10_21_JD_BR_VCas/VcCas/JD_BR_VcSh_4G.bc", header = T, stringsAsFactors = F)
hit_dat <- fread("../Studies/19_10_21_JD_BR_VCas/test_mapping/4G_bam_parser_test2.txt", header = T, stringsAsFactors = F)

bc_dat$read <- sub(" 1:N:0:ATCACG","", bc_dat$read)

merge_dat <- merge(hit_dat, bc_dat, by.x = "Read", by.y ="read", all.x = T)


### Subset and filter
e_coli_hits <- subset(merge_dat, GENOME1 == "Escherichia_coli_BL21_CP001509.3")

hits_filt <- subset(e_coli_hits, MAPQ1 > 0 & MAPQ2 > 0 & NM1 < 5 & NM2 < 5)


### Parse mapping locations and account for stranded-ness 
hits_pos_strand <- subset(hits_filt, STRAND1 == "+")
loci_pos_strand <- count(hits_pos_strand$START1)
loci_pos_strand$dir <- "+"


hits_neg_strand <- subset(hits_filt, STRAND1 == "-")
loci_neg_strand <- count(hits_neg_strand$END1)
loci_neg_strand$dir <- "-"

loci <- rbind(loci_pos_strand, loci_neg_strand)

loci <- data.frame(coord = loci$x, reads = loci$freq, strand = loci$dir, uniq_bc = NA)



### Identify On Target Locations 
loci$on_targ <- ifelse(between(loci$coord,335441,336008) | between(loci$coord,749728,750361), TRUE, FALSE)



### 
for (i in 1:nrow(loci)){
  
  if(loci$strand[i] == "+"){
    tmp_loci <- subset(hits_pos_strand, START1 == loci[i,1])
    tmp_bc_count <- length(unique(tmp_loci$barcodes))
    loci[i,4] <- tmp_bc_count
  }
  
  if(loci$strand[i] =="-"){
    tmp_loci <- subset(hits_neg_strand, END1 == loci[i,1])
    tmp_bc_count <- length(unique(tmp_loci$barcodes))
    loci[i,4] <- tmp_bc_count
  }
  
  cat(paste0("Processed ",i," Loci\r"))
}

# reads vs unique_bc per position
# ggplot(loci, aes(x = log(reads + 1), y = log(uniq_bc + 1))) +
#   geom_point(alpha = 0.5) +
#   stat_smooth()


# Full plot
ggplot(loci, aes(x = coord, y = log(uniq_bc+1))) +
  geom_bar(stat = "identity", color = "steelblue") +
  xlab("Genome Position (bp)") +
  ylab("Unique Barcodes") 
  


# Target Locus 1 (335741-335708 +- 200bp)
ggplot(loci, aes(x = coord, y = uniq_bc, fill = strand)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_rect(aes(xmin=335708, xmax=335741, ymin=-10, ymax=10), fill="black", inherit.aes = FALSE) +
  xlim(c(335500, 335850)) +
  xlab("Genome Position (bp)") +
  ylab("Unique Barcodes")
  


# Target Locus 2 (750028-750061 +- 200bp)
ggplot(loci, aes(x = coord, y = uniq_bc, fill = strand)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_rect(aes(xmin=750028, xmax=750061, ymin=-10, ymax=10), fill="black", inherit.aes = FALSE) +
  xlim(c(749928, 750261)) +
  xlab("Genome Position (bp)") +
  ylab("Unique Barcodes")






### Analysis of on (within 300bp of target) vs off target hits
#all_samples <- list()
#all_samples[[3]] <- loci

summary_dat <- data.frame()
sample_names <- c("4G", "5G", "6G")

for (i in 1:length(all_samples)){
  tmp <- all_samples[[i]]
  tmp_ag <- aggregate(uniq_bc ~ on_targ, tmp, sum)
  
  tmp_ag$sample <- sample_names[i]
  
  summary_dat <- rbind(tmp_ag, summary_dat)
  
}

ggplot(summary_dat, aes(x = on_targ, y = uniq_bc)) +
  geom_boxplot(fill = "steelblue") +
  scale_y_log10() +
  xlab(NULL) +
  ylab("Unique Insertion Barcodes")
  




### VcCas Fwd Read Only
bc_dat <- fread("../Studies/19_10_21_JD_BR_VCas/VcCas/JD_BR_VcSh_4G.bc", header = T, stringsAsFactors = F)
hit_dat <- fread("../Studies/19_10_21_JD_BR_VCas/test_mapping/4G_fwdonly_maptest.hits", header = F, stringsAsFactors = F)
colnames(hit_dat) <- c("genome", "read", "MAPQ1", "NM1", "length", "START1")

bc_dat$read <- sub(" 1:N:0:ATCACG","", bc_dat$read)

merge_dat <- merge(hit_dat, bc_dat, by.x = "read", by.y ="read", all.x = T)


e_coli_hits <- subset(merge_dat, genome == "Escherichia_coli_BL21_CP001509.3")

e_coli_hits_filt <- subset(e_coli_hits, MAPQ1 > 0 & NM1 < 5)

loci <- count(e_coli_hits_filt$START1)

loci <- data.frame(pos = loci$x, reads = loci$freq, uniq_bc = NA)

loci$on_targ <- ifelse(between(loci$pos,335441,336008) | between(loci$pos,749728,750361), TRUE, FALSE)


for (i in 1:nrow(loci)){
  
  tmp_loci <- subset(e_coli_hits_filt, START1 == loci[i,1])
  tmp_bc_count <- length(unique(tmp_loci$barcodes))
  
  loci[i,3] <- tmp_bc_count
  
}

# reads vs unique_bc per position
ggplot(loci, aes(x = log(reads + 1), y = log(uniq_bc + 1))) +
  geom_point(alpha = 0.5) +
  stat_smooth()


# Full plot
ggplot(loci, aes(x = pos, y = uniq_bc)) +
  geom_bar(stat = "identity", color = "blue") +
  xlab("Genome Position (bp)") +
  ylab("Unique Barcodes")


# Target Locus 1 (335741-335708 +- 300bp)
ggplot(loci, aes(x = pos, y = uniq_bc)) +
  geom_bar(stat = "identity", color = "blue") +
  xlim(c(335441, 336008))


# Target Locus 2 (750028-750061 +- 300bp)
ggplot(loci, aes(x = pos, y = uniq_bc)) +
  geom_bar(stat = "identity", color = "blue") +
  xlim(c(749728, 750361))




### Analysis of on (within 300bp of target) vs off target hits
#all_samples <- list()
#all_samples[[3]] <- loci

summary_dat <- data.frame()
sample_names <- c("4G", "5G", "6G")

for (i in 1:length(all_samples)){
  tmp <- all_samples[[i]]
  tmp_ag <- aggregate(uniq_bc ~ on_targ, tmp, sum)
  
  tmp_ag$sample <- sample_names[i]
  
  summary_dat <- rbind(tmp_ag, summary_dat)
  
}

ggplot(summary_dat, aes(x = on_targ, y = uniq_bc)) +
  geom_boxplot(fill = "steelblue") +
  scale_y_log10() +
  xlab(NULL) +
  ylab("Unique Insertion Barcodes")








### ShCas

bc_dat <- fread("../Studies/19_10_21_JD_BR_VCas/ShCas/JD_BR_VcSh_14G.bc", header = T, stringsAsFactors = F)
hit_dat <- fread("../Studies/19_10_21_JD_BR_VCas/ShCas/JD_BR_VcSh_14G.hits", header = F, stringsAsFactors = F)
colnames(hit_dat) <- c("genome", "read", "MAPQ2", "NM2", "length", "START2")

bc_dat$read <- sub(" 1:N:0:ACAGTG","", bc_dat$read)

merge_dat <- merge(hit_dat, bc_dat, by.x = "read", by.y ="read", all.x = T)


e_coli_hits <- subset(merge_dat, genome == "Escherichia_coli_BL21_CP001509.3")

e_coli_hits_filt <- subset(e_coli_hits, MAPQ2 > 0 & NM2 < 5 & length >= 50)

loci <- count(e_coli_hits_filt$START2)

loci <- data.frame(pos = loci$x, reads = loci$freq, uniq_bc = NA)

loci$on_targ <- ifelse(between(loci$pos,335441,336008) | between(loci$pos,749728,750361), TRUE, FALSE)


for (i in 1:nrow(loci)){
  
  tmp_loci <- subset(e_coli_hits_filt, START2 == loci[i,1])
  tmp_bc_count <- length(unique(tmp_loci$barcodes))
  
  loci[i,3] <- tmp_bc_count
  
}

# Full plot
ggplot(loci, aes(x = pos, y = uniq_bc)) +
  geom_bar(stat = "identity", color = "blue") +
  xlab("Genome Position (bp)") +
  ylab("Unique Barcodes")


# Target Locus 1 (335741-335708 +- 300bp)
ggplot(loci, aes(x = pos, y = uniq_bc)) +
  geom_bar(stat = "identity", color = "blue") +
  xlim(c(335441, 336008))


# Target Locus 2 (750028-750061 +- 300bp)
ggplot(loci, aes(x = pos, y = uniq_bc)) +
  geom_bar(stat = "identity", color = "blue") +
  xlim(c(749728, 750361))



### Analysis of on (within 300bp of target) vs off target hits
#all_samples <- list()
#all_samples[[2]] <- loci

summary_dat <- data.frame()
sample_names <- c("13G", "14G", "15G")

for (i in 1:length(all_samples)){
  tmp <- all_samples[[i]]
  tmp_ag <- aggregate(uniq_bc ~ on_targ, tmp, sum)
  
  tmp_ag$sample <- sample_names[i]
  
  summary_dat <- rbind(tmp_ag, summary_dat)
  
}

ggplot(summary_dat, aes(x = on_targ, y = uniq_bc)) +
  geom_boxplot(fill = "steelblue") +
  scale_y_log10() +
  xlab(NULL) +
  ylab("Unique Insertion Barcodes")

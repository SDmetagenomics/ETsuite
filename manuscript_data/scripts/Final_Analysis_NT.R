library(data.table)
library(ggplot2)
library(MASS)
library(dplyr)
library(rcompanion)

'%notin%' <- Negate('%in%')


### LOAD ALL DATA FOR FIG 1C + 2A (NOVASEQ + METAGENOMICS)

## ET-seq Data 
dat <- fread("~/Dropbox/Banfield_Lab_Files/Projects/mCAFE/mCAFE_Project/Studies/20_2_27_BR_30_ET_KlebCurve_NT/NT_FINAL_ANALYSIS/ETstats_hc_output.txt")
dat <- dat[,-c(4:7,10,11)]

## Metagenomics Data - load
mg_dat <- fread("~/Dropbox/Banfield_Lab_Files/Projects/mCAFE/mCAFE_Project/Studies/20_2_27_BR_30_ET_KlebCurve_NT/NT_FINAL_ANALYSIS/genome_coverage.txt")


## Save Data for github upload
fwrite(dat, "~/Dropbox/Banfield_Lab_Files/Projects/mCAFE/mCAFE_Project/ETsuite/manuscript_data/data/Nat_Trans_et.txt", sep = "\t", quote = F)
fwrite(mg_dat, "~/Dropbox/Banfield_Lab_Files/Projects/mCAFE/mCAFE_Project/ETsuite/manuscript_data/data/Nat_Trans_mg.txt", sep = "\t", quote = F)




### APPLY MIN ABUNDANCE FILTERING (MG == 10 reads ; ET == 2 BC)

## Make Seperate data.frames for filtering
dat_filt <- dat
mg_filt <- mg_dat

## ET-filter - Must have >= 2 uique barcodes
dat_filt$UBCdC <- ifelse(dat_filt$UBCdC < 2, 0, dat_filt$UBCdC)
dat_filt$RDSdC <- ifelse(dat_filt$UBCdC == 0, 0, dat_filt$RDSdC)

## MG-filter - Must have >= 10 reads
mg_filt$TOT_RDS <- ifelse(mg_filt$TOT_RDS < 10, 0, mg_filt$TOT_RDS)
mg_filt$MEAN_COV <- ifelse(mg_filt$TOT_RDS == 0, 0, mg_filt$MEAN_COV)




### CALCULATE AND INTEGRATE SIZE FACTORS - ET-SEQ Hits

## Manual Size Factor Calulation raito of b.theta to geometric mean of b.theta
ET_sio_dat <- subset(dat_filt, GENOME == "Bacteroides_thetaiotaomicron_VPI_5482_GCF_000011065.1")

## Size factor calculated for SIO reads and barcodes (Sample_i/Geometric mean)
ET_sio_dat$ET_sf_rds <- ET_sio_dat$RDSdC/exp(mean(log(ET_sio_dat$RDSdC)))
ET_sio_dat$ET_sf_ubc <- ET_sio_dat$UBCdC/exp(mean(log(ET_sio_dat$UBCdC)))

## Create Size Factor Normalized Data
dat_filt <- merge(dat_filt, ET_sio_dat[,c("SAMPLE", "ET_sf_rds", "ET_sf_ubc")], by = "SAMPLE", all.x = T)
dat_filt <- data.frame(dat_filt,
                       RDS_NRM = dat_filt$RDSdC/dat_filt$ET_sf_rds,
                       UBC_NRM = dat_filt$UBCdC/dat_filt$ET_sf_ubc)




### CALCULATE AND INTEGRATE SIZE FACTORS - Metagenomics 

## Manual Size Factor Calulation raito of b.theta to geometric mean of b.theta
MG_sio_dat <- subset(mg_filt, GENOME == "Bacteroides_thetaiotaomicron_VPI_5482_GCF_000011065.1")

## Size factor calculated for SIO reads (Sample_i/Geometric mean)
MG_sio_dat$MG_sf_cov <- MG_sio_dat$MEAN_COV/exp(mean(log(MG_sio_dat$MEAN_COV)))

## Create Size Factor Normalized Data
mg_filt <- merge(mg_filt, MG_sio_dat[,c("SAMPLE", "MG_sf_cov")], by = "SAMPLE", all.x = T)
mg_filt <- data.frame(mg_filt,
                      COV_NRM = mg_filt$MEAN_COV/mg_filt$MG_sf_cov)




## GENERATE ET COUNTS NORMALIZED BY METAGENOME DATA

## Create Master Data Sheet
master_dat <- merge(dat_filt, mg_filt[,c("SAMPLE", "GENOME", "COV_NRM", "COM_FRAC")], by = c("SAMPLE", "GENOME"), all.x = T)


## Filter out things we dont want to quantify
master_dat_filt <- subset(master_dat, GENOME %notin% c("pHLL250", "pHLL249", "Escherichia_coli_BW25113_GCF_000750555.1", "E.coli_WM3604_scaffold_min1000", "pBFC0882_Barcoded_VcCas", "Sinorhizobium_meliloti_1021_GCF_000006965.1"))


## Create Basic Normalized Values

# Direct Division of normalized reads/barcodes by Metagenomic RPK (Results in some NA due to 0/0 ; Also not whole numbers so cant be used for glm.nb but could be used for Gamma)
master_dat_filt$ET_RDS_final <- master_dat_filt$RDS_NRM / master_dat_filt$COV_NRM
master_dat_filt$ET_UBC_final <- master_dat_filt$UBC_NRM / master_dat_filt$COV_NRM
master_dat_filt[is.na(master_dat_filt)] <- 0

# Collect B. theta normlaized values (B. theta = 100 % and divide each organisms normalized values by B. theta * 100 to give percent)
sio_norm <- subset(master_dat_filt, GENOME == "Bacteroides_thetaiotaomicron_VPI_5482_GCF_000011065.1")
master_dat_filt <- subset(master_dat_filt, GENOME != "Bacteroides_thetaiotaomicron_VPI_5482_GCF_000011065.1")

# Divide normalized counts by constant of B. theta normalized counts to get transformation efficency estimate for Reads and Barcodes
master_dat_filt$ET_RDS_final_bt <- master_dat_filt$ET_RDS_final / sio_norm$ET_RDS_final
master_dat_filt$ET_UBC_final_bt <- master_dat_filt$ET_UBC_final / sio_norm$ET_UBC_final

# Multiply Organism Transformation Efficnecy (READS) by Relative Abundance - Gives Fraction in total sample
master_dat_filt$Frac_in_Samp <- master_dat_filt$ET_RDS_final_bt * master_dat_filt$COM_FRAC

# Take log of ET_RDS_final_bt for later use
master_dat_filt$log_ET_RDS_final_bt <- log10(master_dat_filt$ET_RDS_final_bt)

### Determine if fractional values vs total culture are below LOD and impute LOD

## Put in Logical values 
## Put in Logical values 
ET_LOD <- 8.408345e-06
ET_LOQ <- 2.179848e-05
master_dat_filt$Above_LOD <- ifelse(master_dat_filt$Frac_in_Samp >= ET_LOD, TRUE, FALSE)
master_dat_filt$Above_LOQ <- ifelse(master_dat_filt$Frac_in_Samp >= ET_LOQ, TRUE, FALSE)


## Subset Kleb Only for Plotting 
master_dat_filt <- subset(master_dat_filt, GENOME == "Klebsiella_michiganensis_M5a1_GCF_002090195.1")


### Extended Data Figure XXX - Plot Stuff - Establish Quantitation of Natural Transformation Experiment (exported as 5inX7in landscape)

ggplot(master_dat_filt, aes(x = reorder(DESC, -Frac_in_Samp), y = Frac_in_Samp, fill = EXP)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1),
             shape = 21, size = 3) +
  geom_hline(aes(yintercept = ET_LOQ), size = .5, linetype = "longdash", color = "black") +
  geom_hline(aes(yintercept = ET_LOD), size = .5, linetype = "dashed", color = "black") +
  scale_fill_manual(values = c("grey60", "#7570B3")) +
  scale_y_log10(limits = c(1e-06, 1e-03), labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black"),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_blank()) +
  coord_cartesian(clip = "off") +
  annotation_logticks(outside = T, sides = "l") +
  xlab(NULL) +
  ylab("Fraction Transformed Kleb in Community")



## Figure 2C - Transformation Efficency from Natural Transformatin Experiment, only Experimentals (exported as 4inX6in portrait)
tmp_plot <- subset(master_dat_filt, EXP == "Experimental" & DESC != "KlebOnly_0h")

ggplot(tmp_plot, aes(x = reorder(DESC, -ET_RDS_final_bt), y = ET_RDS_final_bt, alpha = Above_LOD)) +
  geom_point(position = position_jitter(width = 0.1),
            shape = 21, size = 3, fill = "#7570B3") +
  scale_y_log10(limits = c(1e-06, 1e-03), labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black"),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_blank()) +
  coord_cartesian(clip = "off") +
  annotation_logticks(outside = T, sides = "l") +
  scale_alpha_discrete(range = c(0.5, 1)) +
  xlab(NULL) +
  ylab("Transformation Efficency")




### OLDER PLOTS NOT USED 

## Plot Stuff - Transformation Efficency 
ci_dat <- subset(master_dat_filt, DESC == "Com_3h")
ci_means <- groupwiseMean(log_ET_RDS_final_bt ~ DESC, data = ci_dat, conf = 0.95)
ci_means$mean <- 10^ci_means$Mean
ci_means$lower <- 10^ci_means$Trad.lower
ci_means$upper <- 10^ci_means$Trad.upper

ggplot(ci_means, aes(x = reorder(DESC, -mean), y = mean)) +
  geom_point(data = master_dat_filt, aes(x = DESC, y = ET_RDS_final_bt),
             position = position_jitter(width = .1),
             size = 2) +
  geom_linerange(aes(x = DESC, ymin = lower, ymax = upper)) +
  geom_point(shape = 21, size = 4, fill = "steelblue", alpha = 0.7) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_log10(limits = c(1e-06, 1e-03)) +
  xlab(NULL) +
  ylab("Organism Transformation Efficency")



library(data.table)
library(ggplot2)
library(MASS)
library(dplyr)
#library(fitdistrplus)

'%notin%' <- Negate('%in%')


### LOAD ALL DATA FOR FIG 1B + 2A (NOVASEQ + METAGENOMICS)

## ET-seq Data 
dat <- fread("~/Dropbox/Banfield_Lab_Files/Projects/mCAFE/mCAFE_Project/Studies/20_4_6_Kleb_Std_Curve_Final/Kcurve_FINAL_ANALYSIS/ETstats_hc_output.txt")
dat <- dat[,-c(4:7)]

## Metagenomics Data - load
mg_dat <- fread("~/Dropbox/Banfield_Lab_Files/Projects/mCAFE/mCAFE_Project/Studies/20_4_6_Kleb_Std_Curve_Final/Kcurve_FINAL_ANALYSIS/MG_INPUT_ETMAPPER.txt")

## Save Data for github upload
#fwrite(dat, "~/Dropbox/Banfield_Lab_Files/Projects/mCAFE/mCAFE_Project/ETsuite/manuscript_data/data/Kcurve_et.txt", sep = "\t", quote = F)
#fwrite(mg_dat, "~/Dropbox/Banfield_Lab_Files/Projects/mCAFE/mCAFE_Project/ETsuite/manuscript_data/data/Kcurve_mg.txt", sep = "\t", quote = F)



### APPLY FILTERING (MG == 10 reads ; ET == 2 BC)

## Make Seperate data.frames for filtering
dat_filt <- dat
mg_filt <- mg_dat ## MG Filter not needed but keeping for consistency with other code 

## Remove Genomes from ETseq Data that are not Kleb or B. theta
dat_filt <- subset(dat_filt, GENOME == "Bacteroides_thetaiotaomicron_VPI_5482_GCF_000011065.1" | GENOME == "Klebsiella_michiganensis_M5a1_GCF_002090195.1")

## ET-filter - Must have >= 2 uique barcodes
dat_filt$UBCdC <- ifelse(dat_filt$UBCdC < 2, 0, dat_filt$UBCdC)
dat_filt$RDSdC <- ifelse(dat_filt$UBCdC == 0, 0, dat_filt$RDSdC)




### CALCULATE AND INTEGRATE VARIOUS SIZE FACTORS - ET-SEQ Hits

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




### CALCULATE AND INTEGRATE VARIOUS SIZE FACTORS - Metagenomics 

## Manual Size Factor Calulation raito of b.theta to geometric mean of b.theta
MG_sio_dat <- subset(mg_filt, GENOME == "Bacteroides_thetaiotaomicron_VPI_5482_GCF_000011065.1")

## Size factor calculated for SIO reads (Sample_i/Geometric mean)
MG_sio_dat$MG_sf_cov <- MG_sio_dat$MEAN_COV/exp(mean(log(MG_sio_dat$MEAN_COV)))

## Create Size Factor Normalized Data
mg_filt <- merge(mg_filt, MG_sio_dat[,c("SAMPLE", "MG_sf_cov")], by = "SAMPLE", all.x = T)
mg_filt <- data.frame(mg_filt,
                      COV_NRM = mg_filt$MEAN_COV/mg_filt$MG_sf_cov)




### GENERATE ET COUNTS NORMALIZED BY METAGENOME DATA

## Create Master Data Sheet
master_dat <- merge(dat_filt, mg_filt[,c("SAMPLE", "GENOME", "COV_NRM", "COM_FRAC")], by = c("SAMPLE", "GENOME"), all.x = T)


## Create Basic Normalized Values

# Direct Division of normalized reads/barcodes by Metagenomic COV
master_dat$ET_RDS_final <- master_dat$RDS_NRM / master_dat$COV_NRM
master_dat$ET_UBC_final <- master_dat$UBC_NRM / master_dat$COV_NRM

# Collect B. theta normlaized values (B. theta = 100 % and divide each organisms normalized values by B. theta * 100 to give percent)
sio_norm <- subset(master_dat, GENOME == "Bacteroides_thetaiotaomicron_VPI_5482_GCF_000011065.1")
master_dat <- subset(master_dat, GENOME != "Bacteroides_thetaiotaomicron_VPI_5482_GCF_000011065.1")

# Divide normalized counts by constant of B. theta normalized counts to get transformation efficency estimate for Reads and Barcodes
master_dat$ET_RDS_final_bt <- master_dat$ET_RDS_final / sio_norm$ET_RDS_final
master_dat$ET_UBC_final_bt <- master_dat$ET_UBC_final / sio_norm$ET_UBC_final

# Multiply Normalized Counts and Concentration by Relative Abundance - Gives Fraction in total sample
master_dat$CONC_in_Samp <- master_dat$CONC * master_dat$COM_FRAC
master_dat$Frac_in_Samp <- master_dat$ET_RDS_final_bt * master_dat$COM_FRAC




### AGGREGATE AND TAKE MEANS OF DATA 

## Rename samples for merging
master_dat$SAMPLE <- sub("B_Kcurve", "A_Kcurve", master_dat$SAMPLE)

## Aggregate sequencing technical replicates + Output Simplifed table --> Based only on ET_RDS_bt
tmp_final_1 <- aggregate(ET_RDS_final_bt ~ SAMPLE + CONC, master_dat, mean)

## Aggregate sequencing technical replicates + Output Simplifed table --> Based only on Frac_in_Samp
tmp_final_2 <- aggregate(Frac_in_Samp ~ SAMPLE + CONC_in_Samp + CONC, master_dat, mean)


## Create Means of run technical replicates - organism fraction
final_means <- aggregate(ET_RDS_final_bt ~ CONC, tmp_final_1, mean)


## Create Means of run technical replicates - culture fraction
final_means$CONC_in_Samp <- aggregate(CONC_in_Samp ~ CONC, tmp_final_2, mean)[,2]
final_means$Frac_in_Samp <- aggregate(Frac_in_Samp ~ CONC, tmp_final_2, mean)[,2]


#fwrite(final_means, "~/Desktop/ET_Seq_Figs/Std_curve_means.txt", quote = F, sep = "\t")


### Calculate Linear Regression to identify LOD 
## LOD is calculated based on Ability to Detect Fraction of Edited Cells Out of Total Cell Population

## THIS IS OPTIONAL AND ALLOWS YOU TO SWITCH ON AN OFF ANALYSIS WITH ZEROS
final_means <- subset(final_means, CONC != 0)

## Quick Plot of Data
qplot(final_means$CONC_in_Samp, final_means$Frac_in_Samp, log = "xy") + 
  stat_smooth(method = "lm", se = F, linetype = 2)

## lm of data used for calculating LOD & LOQ
mod <- lm(Frac_in_Samp ~ CONC_in_Samp , data = final_means)
#plot(mod)
mean(residuals(mod))
summary(mod)

## LOD is 3.3 * std err of intercept / slope
LOD <- 3.3 * coef(summary(mod))[1,2] / coef(mod)[2]
LOD

## ET LOD is the min value that ETseq data can return that constitutes a positive detection
ET_LOD <- predict(mod, data.frame(CONC_in_Samp = LOD))
ET_LOD

## LOQ is 10 * std err of intercept / slope
LOQ <- 10 * coef(summary(mod))[1,2] / coef(mod)[2]
LOQ

## ET LOQ is the min value that ETseq data can return that constitutes a quantifiable value
ET_LOQ <- predict(mod, data.frame(CONC_in_Samp = LOQ))
ET_LOQ

### Without Zeros
#LOD = 3.086878e-05 
#LOQ = 9.354175e-05
#ET_LOD = 8.408345e-06
#ET_LOQ = 2.179848e-05

### Plotting

## Add predicted values
final_means$pred <- predict(mod)

## Plot Within Organism Response 
ggplot(final_means, aes(x = CONC, y = ET_RDS_final_bt)) +
  geom_point(shape = 21, fill = "steelblue3", size = 5) +
  geom_abline(aes(slope = 1, intercept = 0), size = .2, linetype = 2) +
  scale_x_log10(limits = c(1e-6, 1e-2)) +
  scale_y_log10(limits = c(1e-6, 1e-2)) +
  xlab("Known Klebsiella Fraction") +
  ylab("Measured Klebsiella Fraction") +
  ggtitle("Known vs Measured Transformed Klebsiella Fraction")


## Plot Within Sample Response with LOD 
ggplot(final_means, aes(x = CONC_in_Samp, y = Frac_in_Samp)) +
  #geom_smooth(aes(y = pred), color = "blue", size = 1, linetype = 2, method = "lm") +
  stat_smooth(method = "lm", se = F, size = 1, linetype = 2, color = "steelblue4") +
  geom_point(shape = 21, fill = "steelblue3", size = 5) +
  geom_point(aes(x = LOD, y = ET_LOD), size = 5, shape = 23, fill = "red") +
  geom_point(aes(x = LOQ, y = ET_LOQ), size = 5, shape = 23, fill = "blue") +
  scale_x_log10(limits = c(1e-7, 5e-03)) +
  scale_y_log10(limits = c(1e-7, 5e-03)) +
  xlab("Known Klebsiella Community Fraction") +
  ylab("Measured Klebsiella Community Fraction") +
  ggtitle("Known vs Measured Transformed Klebsiella Fraction of Community")


## Figure 1B - Plot Within Sample Response with LOD (Saved as 4inx6in landscape pdf)
ggplot(final_means, aes(x = CONC_in_Samp, y = Frac_in_Samp)) +
  #geom_smooth(aes(y = pred), color = "blue", size = 1, linetype = 2, method = "lm") +
  stat_smooth(method = "lm", se = F, size = 1, linetype = 1, color = "black") +
  geom_point(shape = 21, fill = "steelblue3", size = 5) +
  geom_segment(aes(x = LOD, y = ET_LOD, xend = LOD, yend = 0), linetype = "dashed", color = "black") +
  geom_segment(aes(x = 0, y = ET_LOD, xend = LOD, yend = ET_LOD), linetype = "dashed", color = "black") +
  geom_segment(aes(x = LOQ, y = ET_LOQ, xend = LOQ, yend = 0), linetype = "longdash", color = "black") +
  geom_segment(aes(x = 0, y = ET_LOQ, xend = LOQ, yend = ET_LOQ), linetype = "longdash", color = "black") +
  scale_x_log10(limits = c(1e-7, 5e-03), labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(limits = c(1e-7, 5e-03), labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black"),
        axis.line = element_line(color = "black"),
        axis.ticks = element_blank()) +
  coord_cartesian(clip = "off") +
  annotation_logticks(outside = T) +
  xlab("Known Klebsiella Community Fraction") +
  ylab("Measured Klebsiella Community Fraction") +
  ggtitle("Known vs Measured Transformed Klebsiella Fraction of Community")



## Output Data Used
fwrite(final_means, "~/Desktop/Std_Curve_dat.txt", quote = F, sep = "\t")








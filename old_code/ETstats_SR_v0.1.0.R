library(data.table)


foo <- fread(file.choose(), sep = "\t")
#foo <- fread(file = "~/Desktop/Hit_Table_Raw.txt", sep = "\t")
#foo <- fread(file = "~/Dropbox/Banfield_Lab_Files/Projects/mCAFE/mCAFE_Project/Studies/18_9_20_Pilot_ETseq_Testing/ET_seq_Barcode_Data/SR_Results/BR2/Hit_Table_Raw.txt", sep = "\t")

colnames(foo) <- c("genome", "read", "mapq", "mis", "length", "pos")


# capture all mapped genome names
gen <- unique(foo$genome)

# generate raw results dataframe
results <- data.frame()

# collect all summary data for unfiltered reads
for (i in 1:length(gen)){

  tmp_genome <- foo[foo$genome == gen[i],]
  
  # calculates stats for each genome - unique hits by unique positions 
  tmp_dat <- data.frame(genome = gen[i],
                        hits = nrow(tmp_genome),
                        uniq_hit = length(unique(tmp_genome$pos)),
                        mean_mq = mean(tmp_genome$mapq),
                        mean_mm = mean(tmp_genome$mis))
  
  results <- rbind(results, tmp_dat)
  
} 




# generate filt results dataframe
results_filt <- data.frame()

# collect all summary data for filtered reads
foo_filt <- subset(foo, mapq >= 20 & mis <=5)

# collect all summary data for filtered reads
for (i in 1:length(gen)){
  
  tmp_genome <- foo_filt[foo_filt$genome == gen[i],]
  
  tmp_dat <- data.frame(genome = gen[i],
                        hits = nrow(tmp_genome),
                        uniq_hit = length(unique(tmp_genome$pos)),
                        mean_mq = mean(tmp_genome$mapq),
                        mean_mm = mean(tmp_genome$mis))
  
  results_filt <- rbind(results_filt, tmp_dat)
  
} 
results_filt <- data.frame(results_filt, fac_hits = (results_filt$hits/sum(results_filt$hits)) * 100)


#tmp <- foo_filt[foo_filt$genome == gen[12],]

write.table(results, "~/Desktop/Res_Raw.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(results_filt, "~/Desktop/Res_Filt.txt", quote = F, sep = "\t", row.names = F, col.names = T)

### Check and Load Libraries
if("seqinr" %in% installed.packages() == F){
  print("Please install R package seqinr. Program quitting...")
  q(save="no")
}
# libraries 
library(seqinr)
library(tools)



### Set up path variables

# Get relative path of ETdb install directory
initial.options <- commandArgs(trailingOnly = F)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)



### Collect and Parse arguments
args <- commandArgs(trailingOnly = T)

# Display help if no args or -h
if("-h" %in% args | !("-g" %in% args) | length(args) == 0) {
  cat("
  
    Usage: ETdb.R -g [genomes] -o [output_dir] [Additonal_Options]

    Mandatory Arguments:
      
      -g: Genomes in fasta format (No Default)

    Additional Options:
    
      -o: Output directory (Default: ./ETdb)
      -s: Don't Produce Genome Stats (Default: TRUE)
      -p: Run prodigal to predict genes (Default: FALSE)
      -cpu: Number of cores (Default: 1)
      -h: Bring up this help menu\n\n")
      
  q(save="no")

}


## Testing Arguments
#genome_folder <- "test_data/ETdb/genomes/"

## Mandatory Arguments

# Genomes in fasta format (No Default)
genome_folder <- args[which(args == "-g") + 1]
genome_files <- list.files(genome_folder)
genome_paths <- normalizePath(paste0(genome_folder,genome_files))
genome_names <- file_path_sans_ext(genome_files)


## Program Control Options

# Output directory (Default: ./ETdb)
out_dir <- "ETdb"
if("-o" %in% args){
  out_dir <- args[which(args == "-o") + 1]
  dir.create(out_dir, recursive = T)
} else{dir.create(out_dir, recursive = T)}

# Produce Genome Stats (Default: FALSE)
stats_run <- TRUE
if("-s" %in% args){
  stats_run <- FALSE
}

# Run prodigal to predict genes (Default: FALSE)
prodigal_run <- FALSE
if("-p" %in% args){
  prodigal_run <- TRUE
}

# Number of cores to use where applicable
cpu <- 1
if("-cpu" %in% args){
  cpu <- as.numeric(args[which(args == "-cpu") + 1])
}




############
############ BEGIN PROGRAM
############



### Concatenate Genomees to file, Build Scaff2Bin, Genome Stats

# output data frames 
scaff2bin <- data.frame()

if(stats_run == TRUE){
  genome_stats <- data.frame()
}


# Process each genome in for loop
for (i in 1:length(genome_files)){
  
  # Indicate current operation
  cat(paste0("Processing Genome: ",genome_names[i],"\n"))
  
  # Read in genome
  tmp_genome <- read.fasta(genome_paths[i], seqtype = "DNA", forceDNAtolower = F)
  
  # Build scaff2bin
  tmp_s2b <- data.frame(scaff = names(tmp_genome),
                        bin = genome_names[i])
  scaff2bin <- rbind(scaff2bin, tmp_s2b)
  
  # Build concatenated genome scaffold file
  write.fasta(tmp_genome,
              names = names(tmp_genome),
              file.out = paste0(out_dir,"/All_Genomes.fa"),
              open = "a")
  
  # Calculate genome statistics if stats_run == T
  if (stats_run == TRUE){
    
    # get genome name, # scaffs, length, and max scaff size
    tmp_name <- genome_names[i]
    tmp_scaf <- length(tmp_genome)
    tmp_size <- sum(getLength(tmp_genome))
    tmp_long <- max(getLength(tmp_genome))
    
    # create tmp vector to hold gc
    gc_tmp <- vector(length = length(tmp_genome), mode = "numeric")
    
    # calculate gc for each contig
    for (j in 1:length(tmp_genome)){
      gc_tmp[j] <- GC(tmp_genome[[j]])
    
    }
    
    # take mean of gc 
    gc_tmp <- round(mean(gc_tmp),digits = 4) * 100
    
    # build tmp stats df for genome
    tmp_stats <- data.frame(Genome = tmp_name,
                            Num_Scaf = tmp_scaf,
                            Size = tmp_size,
                            Long_Scaf = tmp_long,
                            Avg_GC = gc_tmp)  
    
    # row bind tmp_stats to output stats df 
    genome_stats <- rbind(genome_stats, tmp_stats)
    
  } ## End stats run
  
} ## End genome processing for loop


# Write scaff2bin output
write.table(scaff2bin,
            paste0(out_dir,"/scaff2bin.txt"),
            quote = F,
            sep = "\t",
            row.names = F, 
            col.names = F)



### Build bowtie2 index

# Indicate current operation
cat("\nBuilding bowtie2 index...\n")

# Build bowtie2 index from All_Genomes.txt
dir.create(paste0(out_dir,"/bt2"), recursive = T)
system(paste0("bowtie2-build",
              " --threads ",cpu,
              " ",out_dir,"/All_Genomes.fa ",
              out_dir,"/bt2/All_Genomes"),
       ignore.stdout = T,
       ignore.stderr = T)



### Predict Proteins (if desired) - This will gnerate coordinates for later plotting ect

if(prodigal_run == T){
  # Indicate current operation
  cat("Predicting Proteins with Prodigal...\n")

  # Predict proteins for each genome with prodigal
  dir.create(paste0(out_dir,"/Proteins"), recursive = T)

  for (i in 1:length(genome_files)){
    cat(paste0("Predicting Proteins for: ",genome_names[i],"\n"))
    
    system(paste0("prodigal -i ",genome_paths[i],
                  " -a ",out_dir,"/Proteins/",genome_names[i],".faa",
                  " -m -p single"), ignore.stdout = T, ignore.stderr = T)
  
  }
  
  # Count number of proteins identified per genome
  prot_num <- data.frame()
  
  for (i in 1:length(genome_files)){
    tmp_genome <- genome_names[i]
    tmp_count <- system(paste0("grep -c '^>' ",out_dir,"/Proteins/",genome_names[i],".faa"), intern = T)  
    
    tmp_count_df <- data.frame(Genome = tmp_genome, Proteins = tmp_count)
    
    prot_num <- rbind(prot_num, tmp_count_df)
  
  }

  # merge genome_stats with proteins counts
  genome_stats <- merge(genome_stats, prot_num, by = "Genome")
  
}


### Write stats output last so protein number can be concatenated into sheet
if(stats_run == TRUE){
  write.table(genome_stats,
              paste0(out_dir,"/genome_stats.txt"),
              quote = F,
              sep = "\t",
              row.names = F, 
              col.names = T)
}














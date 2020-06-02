## VERSION ETdb v1.0

### Check and Load Libraries
# Check data.table
if ("data.table" %in% installed.packages() == FALSE){
  print("Please install R package data.table. Program quitting...")
  q(save="no")
}
# Check seqinr
if("seqinr" %in% installed.packages() == F){
  print("Please install R package seqinr. Program quitting...")
  q(save="no")
}
# libraries
library(data.table)
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
      
      -g: Genome Directory (Files in FASTA format; No Default)

    Additional Options:
    
      -o: Output directory (Default: ./ETdb)
      -d: Don't Produce Concatenated Genome File (Default: Runs)
      -b: Don't Produce Scaff2Bin File (Default: Runs)
      -s: Don't Produce Genome Stats (Default: Runs)
      -x: Don't Create Bowtie2 Index (Default: Runs)
      -p: Don't Run Prodigal to Predict Genes (Default: Runs)
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

# Produce Concatenated Genome (Default: TRUE)
concat_run <- TRUE
if("-d" %in% args){
  concat_run <- FALSE
}

# Produce Scaff2Bin File (Default: TRUE)
s2b_run <- TRUE
if("-b" %in% args){
  s2b_run <- FALSE
}

# Produce Genome Stats (Default: TRUE)
stats_run <- TRUE
if("-s" %in% args){
  stats_run <- FALSE
}

# Produce Bowtie2 Index (Default: TRUE)
bt2_run <- TRUE
if("-x" %in% args){
  bt2_run <- FALSE
}

# Run prodigal to predict genes (Default: TRUE)
prodigal_run <- TRUE
if("-p" %in% args){
  prodigal_run <- FALSE
}

# Number of cores to use where applicable
cpu <- 1
if("-cpu" %in% args){
  cpu <- as.numeric(args[which(args == "-cpu") + 1])
}

# Set DT Threads
setDTthreads(cpu)

# Check for incompatable arguments
if(concat_run == F & bt2_run == T){
  cat("\nERROR: Can't Make Bowtie Index Without Concantenated Genome File, Adjust Program Options\n")
  q()
}




############
############ BEGIN PROGRAM
############



### Concatenate Genomees to file, Build Scaff2Bin, Genome Stats

## output data frames and folders 
if(stats_run == TRUE){
  genome_stats <- data.frame()
}

if(prodigal_run == TRUE){
  dir.create(paste0(out_dir,"/Proteins"), recursive = T)
}

if(bt2_run == TRUE){
  dir.create(paste0(out_dir,"/bt2"), recursive = T)
}



### Process each genome in for loop (concat, s2b, stats, prodigal)
for (i in 1:length(genome_files)){
  
  ## Indicate current operation
  cat(paste0("Processing Genome: ",genome_names[i],"\n"))
  
  ## Read in genome
  tmp_genome <- read.fasta(genome_paths[i], seqtype = "DNA", forceDNAtolower = F)
  
  
  ## Create/Append genome to concatenated genome file if concat_run == T
  if(concat_run == T){
    
    # Build concatenated genome scaffold file
    write.fasta(tmp_genome,
                names = names(tmp_genome),
                file.out = paste0(out_dir,"/All_Genomes.fa"),
                open = "a")
  }
  
  
  ## Create/Append scaff2bin file if s2b_run == T
  if(s2b_run == TRUE){
    
    # Build scaff2bin
    tmp_s2b <- data.frame(scaff = names(tmp_genome),
                        bin = genome_names[i])
    
    fwrite(tmp_s2b,
           file = paste0(out_dir,"/scaff2bin.txt"),
           sep = "\t",
           row.names = F,
           col.names = F,
           append = T)
  }
  
  
  ## Calculate genome statistics if stats_run == T
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
                            C_Name = NA,        ####*** This is where the placeholders for Genome Type and 
                            Type = "Target",    ####*** common name are added 
                            Num_Scaf = tmp_scaf,
                            Size = tmp_size,
                            Long_Scaf = tmp_long,
                            Avg_GC = gc_tmp)  
    
    # row bind tmp_stats to output stats df 
    genome_stats <- rbind(genome_stats, tmp_stats)
    
  } ## End stats run... See end of program for final output 
  
  
  ## Predict proteins with prodigal if prodigal_run == T - This will gnerate coordinates for later plotting ect
  if(prodigal_run == T){
  
    cat(paste0("Predicting Proteins for: ",genome_names[i],"\n"))
      
    system(paste0("prodigal -i ",genome_paths[i],
                  " -a ",out_dir,"/Proteins/",genome_names[i],".faa",
                  " -m -p single"), ignore.stdout = T, ignore.stderr = T)
      
    }
}



### Add number of proteins predicted per genome if prodigal_run == T
if(prodigal_run == T){
  
  # Create df to hold # of proteins identified per genome
  prot_num <- data.frame()
  
  # loop across .faa files and count proteins 
  for (i in 1:length(genome_files)){
    tmp_genome <- genome_names[i]
    tmp_count <- system(paste0("grep -c '^>' ",out_dir,"/Proteins/",genome_names[i],".faa"), intern = T)  
    
    tmp_count_df <- data.frame(Genome = tmp_genome, Proteins = tmp_count)
    
    prot_num <- rbind(prot_num, tmp_count_df)
    
  }
  
  # merge genome_stats with proteins counts
  genome_stats <- merge(genome_stats, prot_num, by = "Genome")
  
  # write genome_stats w/ protien count output
  write.table(genome_stats,
              paste0(out_dir,"/genome_stats.txt"),
              quote = F,
              sep = "\t",
              row.names = F, 
              col.names = T)
  
}else{
  
  write.table(genome_stats,
              paste0(out_dir,"/genome_stats.txt"),
              quote = F,
              sep = "\t",
              row.names = F, 
              col.names = T)
}



### Build bowtie2 index if bt2_run == T
if(bt2_run == T){
  
  # Indicate current operation
  cat("\nBuilding bowtie2 index...\n")

  # Build bowtie2 index from All_Genomes.txt
  system(paste0("bowtie2-build",
                " --threads ",cpu,
                " ",out_dir,"/All_Genomes.fa ",
                out_dir,"/bt2/All_Genomes"),
        ignore.stdout = T,
        ignore.stderr = T)
}


### Indicate Program Completed
cat("\nProgram Complete..:-)\n")











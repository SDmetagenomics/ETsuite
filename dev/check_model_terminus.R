library(data.table)

'%notin%' <- Negate('%in%')

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}


md <- "/Users/sdiamond/Dropbox/Banfield_Lab_Files/Projects/mCAFE/mCAFE_Project/ETsuite/db/ETseq_newprimers_allmodels_v2.fa"


## Recover Model Sequences 
models <- system(paste0("grep -B 1 \"^NN.*\" ",md," | sed 's/>//'"), intern = T)

## Drop grep Applied Divider Characters 
drop_char <- "--"
models <- models[models %notin% drop_char]

## Build Selection Vectors and Subset models + names into model_df 
name_coord <- seq(from = 1, to = length(models) - 1, by = 2)
model_coord <-  seq(from = 2, to = length(models), by = 2)   
model_df <- data.table(model_name = models[name_coord],
                       model_seq = models[model_coord])

## Add 5 Terminal Bases of Each Model to df 
model_df$junc_seq <- substrRight(model_df$model_seq, 5)

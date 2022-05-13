library(Seurat)
library(tidyverse)
library(argparse)
library(data.table)

## SET VARIABLES
parser <-  ArgumentParser()
parser$add_argument('--data',
                    type = "character", 
                    # nargs='+',
                    help = 'path to rdata object')
parser$add_argument('-i',                                   # immune.combined.sct
                    '--obj_integ',
                    type = 'character',
                    help = 'integrated seurat_object')
parser$add_argument('-o',
                    '--obj_name',
                    type = 'character',
                    help = 'name of a seurat_object')
parser$add_argument('--out_file',
                    type = "character",
                    help = 'output filename')

argss <- parser$parse_args()

print(argss)

# load("/mnt/tank/scratch/rasmirnov/code/jac_index/data/pbmc3k/pbmc3k_seurat.RData")
#####!! Don't forget to replace "args" on "argss"
load(argss$data)
seurat_obj<- get(argss$obj_integ)
obj_name <- argss$obj_name
output <- argss$out_file

#### real example
# colnames
res_names <- grep('snn_res', colnames(seurat_obj@meta.data), value = T)
# content: cell idents
all_res <- list()
for (res in res_names) {
  all_res[[res]] <- seurat_obj[[res]][[1]]     # only values without colnames
  names(all_res[[res]]) <- colnames(seurat_obj)    
}

all_res_df <- data.frame(resolution = res_names)
all_res_df$original_ident <- all_res
# leave only digit values
all_res_df$resolution <- gsub("[^[:digit:]., ]", "", all_res_df$resolution)
all_res_df$resolution <- gsub("^[\\.*]", "", all_res_df$resolution)
saveRDS(all_res_df, file = output)

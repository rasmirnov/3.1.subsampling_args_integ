library(tidyverse)
library(argparse)
library(Seurat)

## SET VARIABLES

parser <-  ArgumentParser()
parser$add_argument('--data',
                    type = "character", 
                    nargs='+',
                    help = 'path to rdata object')
# parser$add_argument('--out_file',
#                     type = "character",
#                     help = 'output filename')

args <- parser$parse_args()

print(args)
rdss<- args$data
# added output variable
# output <- args$out_file

get_df<- function(rds){
	res<- readRDS(rds)
	return(res)
}

dat.list<- lapply(rdss, get_df)
# why we don't use args$out_file arg here?
# do.call applyes a function to args
# bind_rows - binds many datasets into one
gather_idents<- do.call(bind_rows, dat.list)
saveRDS(gather_idents, file = "gather_subsample.rds")

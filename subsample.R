library(Seurat)
library(tidyverse)
library(argparse)
library(data.table)

## SET VARIABLES
parser <-  ArgumentParser()
parser$add_argument('--data',
                    type = "character",
                    help = 'path to rdata object')
parser$add_argument('-i',                                   # immune.combined.sct
                    '--obj_integ',
                    type = 'character',
                    help = 'integrated seurat_object')
parser$add_argument('--obj_name',
                    type = 'character',
                    help = 'name of a seurat_object')
# parser$add_argument('--res',
#                     nargs = '+',
#                     type = "character",
#                     help = 'resolution')
parser$add_argument('--run_id',
                    type = "character",
                    help = 'run_id')
parser$add_argument('--rate',
                    type = "double",
                    help = 'rate')

argss <- parser$parse_args()

print(argss)

####### Real params
load(argss$data)
seurat_obj<- get(argss$obj_integ)
obj_name <- argss$obj_name
# resolution<- argss$res
run_id<- argss$run_id
rate <- argss$rate


#подавать на вход датасет с одним резолюш
RandomSubsetData<- function(object, rate, random.subset.seed = NULL, ...){
  ncells<- nrow(object@meta.data)                 # pbmc3k$SCT_snn_res.0.05 %>% data.frame() %>% nrow()
  ncells.subsample<- round(ncells * rate)
  
  set.seed(random.subset.seed)
  
  selected.cells<- sample(colnames(object), ncells.subsample)
  object<- subset(object, cells =  selected.cells,
                  ...)
  return(object)
}

subset_seurat_obj<- RandomSubsetData(seurat_obj, rate = rate)            # =0.8

# original_ident<- Idents(subset_seurat_obj) 
##### Way 2: rows - res, cols - list with vectors of cells 
res_names <- grep('snn_res', colnames(subset_seurat_obj@meta.data), value = T)
# content: cell idents
all_res <- list()
for (res in res_names) {
  all_res[[res]] <- subset_seurat_obj[[res]][[1]]     # only values without colnames
  names(all_res[[res]]) <- colnames(subset_seurat_obj)    
}

all_res_df <- data.frame(resolution = res_names)
all_res_df$original_ident <- all_res

PreprocessSubsetData<- function(object,
                                num.pc = 20,
                                pc.use = NULL,
                                useSCTransform = TRUE,
                                ...){
  if(!is.null(pc.use)){
    if(pc.use > num.pc){
      stop("Specify the maximum pc.use number as less than or equal to the total num.pc")
    }
  }
  meta.data.colnames<- object@meta.data %>% colnames()
  # vars.to.regress<- c("percent.mt","nFeature_RNA")
  # vars.to.regress<- vars.to.regress[vars.to.regress %in% meta.data.colnames]
  if(useSCTransform==TRUE){
    object<- SCTransform(object, 
                         method = "glmGamPoi", 
                         vars.to.regress = "percent.mt", 
                         verbose = FALSE)  
  }else{
    stop("The SCTransform method for normalization is the only method currently supported by this function.  If you wish to use the approach that involves NormalizeData, ScaleData, and FindVariableFeatures and enables use of the Jackstraw procedure for determining which PCs to use please use the PreprocessSubsetDataV2 function from the scclusteval R package.")
    
  }
  object <- FindVariableFeatures(object,           ###
                                 selection.method = "vst", 
                                 nfeatures = 2000)
  all.genes <- rownames(object)                   ###
  object <- ScaleData(object, features = all.genes)     ###
  
  object<- RunPCA(object = object, features = VariableFeatures(object = object))
  if(is.null(pc.use)){
    pc.use <- num.pc
    message("SCTransform is being used and the Jackstraw procedure for determining which PCs to use is not compatable with this procedure. Since pc.use was not specified it is being automatically set to num.pc")
  }
  # add significant pc number to metadata, need to have names same as the cells
  pc.use.meta<- rep(pc.use, length(colnames(object)))
  names(pc.use.meta)<- colnames(object)
  object<- AddMetaData(object = object, metadata = pc.use.meta, col.name = "pc.use")
  object<- FindNeighbors(object, dims = 1:20)
  object <- FindClusters(object = object,                  
                         resolution = c(0.1, 0.2, 0.4, 0.6, 0.8, 1, 1.5, 2, 5))          
  return(object)
}

###########? источник ошибки?
# command<- paste("PreprocessSubsetData", "(", "subset_seurat_obj,",
#                 "resolution=", resolution, ")")
command<- paste("PreprocessSubsetData", "(", "subset_seurat_obj,", ")")

subset_seurat_obj<- eval(parse(text=command))

res_names <- grep('snn_res', colnames(subset_seurat_obj@meta.data), value = T)            # 'ed_snn_res' for integrated
# content: cell idents
all_res2 <- list()
for (res in res_names) {
  all_res2[[res]] <- subset_seurat_obj[[res]][[1]]     # only values without colnames
  names(all_res2[[res]]) <- colnames(subset_seurat_obj)    
}

all_res_df2 <- data.frame(resolution = res_names)
all_res_df2$recluster_ident <- all_res2
all_res_df2$round <- run_id                          

res <- cbind(all_res_df, all_res_df2[ , 2:3])                            # and add round

res$resolution <- gsub("[^[:digit:]., ]", "", res$resolution)
res$resolution <- gsub("^[\\.*]", "", res$resolution)

outfile<- paste0("subsample/subsample_", "round_", run_id, ".rds")
# dir.create("/mnt/tank/scratch/rasmirnov/code/jac_index/scripts/argp_subsampling2/subsample", recursive = T)
# outfile<- paste0("/mnt/tank/scratch/rasmirnov/code/jac_index/scripts/argp_subsampling2/",
#                  "subsample/subsample_", "round_", run_id, ".rds")
saveRDS(res, file = outfile)

## make sure it is not empty file
info<- file.info(outfile)
if (info$size == 0) {
  quit(status = 1)
}

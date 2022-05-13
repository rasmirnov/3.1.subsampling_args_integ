configfile: "config.yaml"

# link на исходные скрипты с scclustevel: https://github.com/crazyhottommy/pyflow_seurat_parameter

INPUT_SEURAT = config["input_seurat"]
integ_seurat_obj = config['integ_seurat']
name_seurat_obj = config['name_seurat']
resolutions = config["subsample_resolutions"].strip().split()          # new: split(' ') 
NUM_OF_SUBSAMPLE = config["num_of_subsample"]
rate = config["subsample_rate"]
rscript_path = config['rscript']

#expand equal to the list comprehension in python
SUBSAMPLE_RUN = expand("subsample/subsample_round_{run_id}.rds", \
	run_id = range(NUM_OF_SUBSAMPLE))

TARGETS = []

TARGETS.extend(SUBSAMPLE_RUN)
TARGETS.append("gather_subsample.rds")
TARGETS.append("gather_full_sample.rds")

localrules: all, gather_subsample, gather_full_sample_preprocess
rule all:
    input: TARGETS


## subsample e.g. 80% of the cells and re-do the clustering for n times
rule subsample_cluster:
        input: rds = INPUT_SEURAT
        output: "subsample/subsample_round_{run_id}.rds"
        log: "00log/subsample_round_{run_id}.log"
        shell: """
               {rscript_path}  scripts/subsample.R \
               --data {input.rds} \
               --obj_integ {integ_seurat_obj} \
               --obj_name {name_seurat_obj} \
               --res {resolutions} --run_id {wildcards.run_id} --rate {rate} 2> {log}\   
               """

## gather the subsampled and reclustered cell idents
rule gather_subsample:
        input: rds = SUBSAMPLE_RUN
        output: "gather_subsample.rds"
        log: "00log/gather_subsample.log"
        shell: "{rscript_path}  scripts/gather_subsample.R --data {input.rds}"
        
        
rule gather_full_sample_preprocess: 
        input: rds = INPUT_SEURAT
        output: rds="gather_full_sample.rds"
        log: "00log/gather_full_masha_idents.log"
        shell: "{rscript_path} scripts/gather_full_masha.R --data {input.rds} --obj_integ {integ_seurat_obj} --obj_name {name_seurat_obj} --out_file {output.rds}"   
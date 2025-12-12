import os
import re


configfile: "config/config.yaml"
CFG = config


localrules:
    all,
    build_dlbclone_model,
    dlbclone_predict


def as_r_list(d):
    """Convert Python dict-of-lists to an R list of c(...) vectors."""
    if d is None:
        return "NULL"
    items = []
    for k, v in d.items():
        # v is always a list, so no need for fallback
        r_vec = "c(" + ",".join(f'"{x}"' for x in v) + ")"
        items.append(f"{k}={r_vec}")
    return "list(" + ",".join(items) + ")"

def as_r_c(vec):
    """Convert a Python list to an R c(...) vector"""
    if vec is None or len(vec) == 0:
        return "NULL"
    return "c(" + ",".join(f'"{x}"' for x in vec) + ")"


rule all:
    input:
        #CFG["opt_model_path"] + "/" + CFG["model_name_prefix"] + "_model.rds",
        #CFG["opt_model_path"] + "/" + CFG["model_name_prefix"] + "_umap.uwot"
        CFG["dlbclone_predict"]["pred_dir"] + "/" + CFG["model_name_prefix"] + "_DLBCLone_predictions.tsv"

rule build_dlbclone_model: # Uses only GAMBL metadata and binary mutation data as training for DLBCLone models 
    output:
        opt_model_rds_path = CFG["opt_model_path"] + "/" + CFG["model_name_prefix"] + "_model.rds",
        opt_model_uwot_path = CFG["opt_model_path"] + "/" + CFG["model_name_prefix"] + "_umap.uwot"
    params:
        core_features = as_r_list(CFG["build_dlbclone_model"]["core_features"]),
        truth_classes = as_r_c(CFG["build_dlbclone_model"]["truth_classes"]),
        hidden_features = as_r_c(CFG["build_dlbclone_model"]["hidden_features"])
    container:
        "dlbclone:latest"
    shell:
        """
        Rscript DLBCLone_model.R \
            --opt_model_path {CFG[opt_model_path]} \
            --model_name_prefix {CFG[model_name_prefix]} \
            --core_features '{params.core_features}' \
            --core_feature_multiplier {CFG[build_dlbclone_model][core_feature_multiplier]} \
            --hidden_features '{params.hidden_features}' \
            --truth_classes '{params.truth_classes}' \
            --min_k {CFG[build_dlbclone_model][min_k]} \
            --max_k {CFG[build_dlbclone_model][max_k]}
        """

rule dlbclone_predict:
    input:
        mutation_matrix = CFG["dlbclone_predict"]["test_data_dir"] # first column is sample ID and all other columns are features
    output:
        predictions = CFG["dlbclone_predict"]["pred_dir"] + "/" + CFG["model_name_prefix"] + "_DLBCLone_predictions.tsv"
    container:
        "dlbclone:latest"
    shell:
        """
        Rscript DLBCLone_predict.R \
            --test_features {input.mutation_matrix} \
            --output_dir {output.predictions} \
            --model_path {CFG[opt_model_path]} \
            --model_prefix {CFG[model_name_prefix]} \
            --fill_missing {CFG[dlbclone_predict][fill_missing]} \
            --drop_extra {CFG[dlbclone_predict][drop_extra]}
        """

# RUN IT
# cd C:\Users\lukeg\Documents\Morin_Lab\DLBCLone_Pipeline\Modules\DLBCLone
# conda activate snakemake
# docker run -it -v C:/Users/lukeg/Documents/Morin_Lab:/mnt/work -w /mnt/work/DLBCLone_Pipeline/Modules/DLBCLone dlbclone:latest snakemake -s DLBCLone.smk -j 8

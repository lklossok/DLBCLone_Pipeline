import os
import re


configfile: "config/config.yaml"
CFG = config()


localrules:
    all,
    build_dlbclone_model,
    test_meta_maf_formatter,
    assemble_genetic_features,
    dlbclone_predict,


def as_r_list(d): # handle list objects from configfile
    if d is None:
        return "NULL"
    parts = [f'{k}="{v}"' for k, v in d.items()]
    return "list(" + ",".join(parts) + ")"


rule all:
    input:
        #CFG["opt_model_path"] + "/" + CFG["model_name_prefix"] + "_model.rds",
        #CFG["opt_model_path"] + "/" + CFG["model_name_prefix"] + "_umap.uwot"
        #"data/test_metadata_formatted.tsv",
        #"data/test_maf_formatted.maf"
        CFG["dlbclone_predict"]["test_data_dir"]
        #CFG["dlbclone_predict"]["pred_dir"] + "/" + CFG["model_name_prefix"] + "_DLBCLone_predictions.tsv"

rule test_meta_maf_formatter:
    input:
        test_metadata = CFG["test_meta_maf_formatter"]["test_metadata_dir"],
        test_maf = CFG["test_meta_maf_formatter"]["test_maf_dir"]
    output:
        form_metadata = temp("data/test_metadata_formatted.tsv"),
        form_maf = temp("data/test_maf_formatted.maf")
    params:
        sv_from_metadata = ",".join(CFG["test_meta_maf_formatter"]["sv_from_metadata"]),
        translocation_status = as_r_list(CFG["test_meta_maf_formatter"]["translocation_status"]),
    container:
        "dlbclone:latest"
    shell:
        """
        Rscript Rscript/test_meta_maf_formatter.R \
            --metadata {input.test_metadata} \
            --maf {input.test_maf} \
            --output_meta_dir {output.form_metadata} \
            --output_maf_dir {output.form_maf} \
            --maf_sample_id_colname {CFG[test_meta_maf_formatter][maf_sample_id_colname]} \
            --metadata_sample_id_colname {CFG[test_meta_maf_formatter][metadata_sample_id_colname]} \
            --sv_from_metadata '{params.sv_from_metadata}' \
            --translocation_status '{params.translocation_status}' \
            --truth_column {CFG[test_meta_maf_formatter][truth_column]} \
            --truth_column_colname {CFG[test_meta_maf_formatter][truth_column_colname]}
        """

rule assemble_genetic_features:
    input:
        form_metadata = rules.test_meta_maf_formatter.output.metadata_formatted,
        form_maf = rules.test_meta_maf_formatter.output.maf_formatted
    output:
        test_mutation_matrix = CFG["dlbclone_predict"]["test_data_dir"]
    params:
        sv_from_metadata = ",".join(CFG["test_meta_maf_formatter"]["sv_from_metadata"])
    container: 
        "dlbclone:latest"
    shell:
        """
        Rscript Rscript/assemble_genetic_features.R \
            --form_metadata {input.form_metadata} \
            --form_maf {input.form_maf} \
            --output_matrix_dir {output.test_mutation_matrix} \
            --sv_from_metadata '{params.sv_from_metadata}'  
        """

rule build_dlbclone_model:
    input:
        training_matrix = CFG["build_dlbclone_model"]["training_matrix_path"], # first column is sample ID and all other columns are features
        training_metadata = CFG["build_dlbclone_model"]["metadata_path"]
    output:
        opt_model_rds_path = CFG["opt_model_path"] + "/" + CFG["model_name_prefix"] + "_model.rds",
        opt_model_uwot_path = CFG["opt_model_path"] + "/" + CFG["model_name_prefix"] + "_umap.uwot"
    params:
        core_features = as_r_list(CFG["build_dlbclone_model"]["core_features"])
    container:
        "dlbclone:latest"
    shell:
        """
        Rscript Rscript/DLBCLone_model.R \
            --metadata {input.training_metadata} \
            --training_matrix {input.training_matrix} \
            --opt_model_path {CFG[opt_model_path]} \
            --model_name_prefix {CFG[model_name_prefix]}
        """

rule dlbclone_predict:
    input:
        mutation_matrix = CFG["dlbclone_predict"]["test_data_dir"], # first column is sample ID and all other columns are features
        model_path = CFG["opt_model_path"] 
        model_prefix = CFG["model_name_prefix"]
    output:
        predictions = CFG["dlbclone_predict"]["pred_dir"] + "/" + CFG["model_name_prefix"] + "_DLBCLone_predictions.tsv"
    container:
        "dlbclone:latest"3
    shell:
        """
        Rscript Rscript/DLBCLone_predict.R \
            --test_features {input.mutation_matrix} \
            --output_dir {output.predictions} \
            --model_path {input.model_path} \
            --model_prefix {input.model_prefix} \
            --fill_missing {CFG[dlbclone_predict][fill_missing]} 
        """

# RUN IT
# cd C:\Users\lukeg\Documents\Morin_Lab\DLBCLone_Pipeline\Modules\DLBCLone
# conda activate snakemake
# docker run -it -v C:/Users/lukeg/Documents/Morin_Lab:/mnt/work -w /mnt/work/DLBCLone_Pipeline/Modules/DLBCLone dlbclone:latest snakemake -s DLBCLone.smk -j 8

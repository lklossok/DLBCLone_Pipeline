# Load necessary libraries
suppressWarnings(suppressMessages(library(optparse)))
suppressWarnings(suppressMessages(library(GAMBLR.predict)))
suppressWarnings(suppressMessages(library(uwot)))
suppressWarnings(suppressMessages(library(tidyverse)))


# Define command line options
option_list <- list(
    make_option(c("-m", "--metadata"), type = "character", help = "path to metadata file with truth labels"),
    make_option(c("-t", "--training_matrix"), type = "character", help = "path to training mutation matrix file"),
    make_option(c("-p", "--opt_model_path"), type = "character", help = "Output path for optimized model"),
    make_option(c("-n", "--model_name_prefix"), type = "character", help = "Name prefix for optimized model"),
    
    make_option(c("-a", "--truth_column"), type = "character", default="lymphgen", help = "Column name in metadata file with truth labels"),
    make_option(c("-b", "--core_features"), type = "character", default="list(ST2=c('SGK1','DUSP2','TET2','SOCS1'),N1='NOTCH1', EZB=c('EZH2','BCL2_SV'), MCD=c('MYD88HOTSPOT','CD79B','PIM1'), BN2=c('BCL6_SV','NOTCH2','SPEN','CD70'))", help = "Core features as a list in R syntax"),
    make_option(c("-c", "--core_feature_multiplier"), type = "double", default=1.5, help = "Multiplier for core features during UMAP construction"),
    make_option(c("-d", "--hidden_features"), type = "character", default=NULL, help = "Hidden features as a vector in R syntax"),
    make_option(c("-e", "--truth_classes"), type = "character", default="c('EZB','MCD','ST2','N1','BN2','Other')", help = "Truth classes as a vector in R syntax"),
    make_option(c("-f", "--other_class"), type = "character", default="Other", help = "Name of the 'Other' class"),
    make_option(c("-g", "--min_k"), type = "integer", default=3, help = "Minimum k to test during optimization"),
    make_option(c("-i", "--max_k"), type = "integer", default=23, help = "Maximum k to test during optimization"),
)

opt <- parse_args(OptionParser(option_list=option_list))

metadata_path <- opt$metadata
training_matrix <- opt$training_matrix
opt_model_path <- opt$opt_model_path
model_name_prefix <- opt$model_name_prefix
truth_column <- opt$truth_column
target_weight <- opt$target_weight
core_features <- eval(parse(text=opt$core_features))
core_feature_multiplier <- opt$core_feature_multiplier
hidden_features <- if (!is.null(opt$hidden_features)) eval(parse(text = opt$hidden_features)) else NULL
truth_classes <- eval(parse(text=opt$truth_classes))
other_class <- opt$other_class
min_k <- opt$min_k
max_k <- opt$max_k

# Load data
metadata <- readr::read_tsv(metadata_path) %>%
    dplyr::filter(.data[[truth_column]] %in% truth_classes)

mutation_matrix <- readr::read_tsv(training_matrix) %>%
    tibble::column_to_rownames("sample_id")


umap_model <- make_and_annotate_umap(
    df = mutation_matrix,
    metadata = metadata,
    truth_column = truth_column,
    target_weight = target_weight,
    core_features = core_features,
    core_feature_multiplier = core_feature_multiplier,
    hidden_features = hidden_features
)

opt_model <- DLBCLone_optimize_params(  
    combined_mutation_status_df = umap_model$features, 
    metadata_df = metadata,
    umap_out = umap_model,
    truth_classes = truth_classes,
    other_class = other_class,
    truth_column = truth_column,
    min_k = min_k,
    max_k = max_k
)

DLBCLone_save_optimized( # saving optimized model
    DLBCLone_model = opt_model,
    base_path = opt_model_path,
    name_prefix = model_name_prefix
)

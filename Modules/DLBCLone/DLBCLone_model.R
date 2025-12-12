# Load necessary libraries
suppressWarnings(suppressMessages(library(optparse)))
suppressWarnings(suppressMessages(library(GAMBLR.predict)))
suppressWarnings(suppressMessages(library(uwot)))


# Define command line options
option_list <- list(
    make_option(c("-p", "--opt_model_path"), type = "character", help = "Output path for optimized model"),
    make_option(c("-n", "--model_name_prefix"), type = "character", help = "Name prefix for optimized model"),
    make_option(c("-b", "--core_features"), type = "character", help = "Core features as a list in R syntax"),
    make_option(c("-c", "--core_feature_multiplier"), type = "double", help = "Multiplier for core features during UMAP construction"),
    make_option(c("-d", "--hidden_features"), type = "character", help = "Hidden features as a vector in R syntax"),
    make_option(c("-e", "--truth_classes"), type = "character", help = "Truth classes as a vector in R syntax"),
    make_option(c("-g", "--min_k"), type = "integer", help = "Minimum k to test during optimization"),
    make_option(c("-i", "--max_k"), type = "integer", help = "Maximum k to test during optimization")
)

opt <- parse_args(OptionParser(option_list=option_list))

opt_model_path <- opt$opt_model_path
model_name_prefix <- opt$model_name_prefix
core_features <- eval(parse(text=opt$core_features))
core_feature_multiplier <- opt$core_feature_multiplier
hidden_features <- if (!is.null(opt$hidden_features)) eval(parse(text = opt$hidden_features)) else NULL
truth_classes <- eval(parse(text=opt$truth_classes))
min_k <- opt$min_k
max_k <- opt$max_k

# Load GAMBL data
metadata <- readr::read_tsv(system.file("extdata/dlbcl_meta_with_dlbclass.tsv",package = "GAMBLR.predict")) %>%
    dplyr::filter(lymphgen %in% truth_classes)

mutation_matrix <- readr::read_tsv(system.file("extdata/all_full_status.tsv",package = "GAMBLR.predict")) %>%
  tibble::column_to_rownames("sample_id") 


umap_model <- make_and_annotate_umap(
    df = mutation_matrix,
    metadata = metadata,
    core_features = core_features,
    core_feature_multiplier = core_feature_multiplier,
    hidden_features = hidden_features
)

opt_model <- DLBCLone_optimize_params(  
    combined_mutation_status_df = umap_model$features, 
    metadata_df = metadata,
    umap_out = umap_model,
    truth_classes = truth_classes,
    min_k = min_k,
    max_k = max_k
)

DLBCLone_save_optimized( # saving optimized model
    DLBCLone_model = opt_model,
    base_path = opt_model_path,
    name_prefix = model_name_prefix
)

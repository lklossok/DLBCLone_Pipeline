# Load necessary libraries
suppressWarnings(suppressMessages(library(optparse)))
suppressWarnings(suppressMessages(library(GAMBLR.predict)))
suppressWarnings(suppressMessages(library(tidyverse)))


# Define command line options
option_list <- list(
    make_option(c("-f", "--test_features"), type = "character", help = "Input test mutation matrix file path"),
    make_option(c("-o", "--output_dir"), type = "character", help = "Output file path for predictions"),
    make_option(c("-p", "--model_path"), type = "character", help = "Pre-trained model file path"),
    make_option(c("-n", "--model_prefix"), type = "character", help = "Pre-trained model name prefix"),
    make_option(c("-m", "--fill_missing"), action="store_true", default=FALSE, help="if TRUE, any features present in the model but not seen in test_features are added and set to zero. If FALSE, missing model features cause an error"),
    make_option(c("-e", "--drop_extra"), action="store_true", default=FALSE, help="if TRUE, any features present in test_features but not seen during model training are dropped")
)

opt <- parse_args(OptionParser(option_list=option_list))

test_features_path <- opt$test_features
output_dir <- opt$output_dir
model_path <- opt$model_path
model_prefix <- opt$model_prefix
fill_missing <- opt$fill_missing
drop_extra <- opt$drop_extra

test_features <- read_tsv(test_features_path) %>%
    column_to_rownames("sample_id")

# loading saved DLBCLone model
loaded_model <- DLBCLone_load_optimized( 
    path = model_path,
    name_prefix = model_prefix
)

# activating loaded DLBCLone model
active_model <- DLBCLone_activate( 
    loaded_model, 
    force = TRUE
)

# predicting DLBCL genetic subgroups for test samples
predictions <- DLBCLone_predict(
    mutation_status = test_features, 
    optimized_model = active_model,
    fill_missing = fill_missing,
    drop_extra = drop_extra
)

write_tsv(predictions$prediction, output_dir)

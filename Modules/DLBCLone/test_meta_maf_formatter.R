# Load necessary libraries
suppressWarnings(suppressMessages(library(optparse)))
suppressWarnings(suppressMessages(library(tidyverse)))

option_list <- list(
    make_option(c("-m", "--metadata"), type = "character", help = "path to test metadata file"),
    make_option(c("-f", "--maf"), type = "character", help = "path to test maf file"),
    make_option(c("-o", "--output_meta_dir"), type = "character", help = "path to meta output"),
    make_option(c("-y", "--output_maf_dir"), type = "character", help = "path to maf output"),
    
    make_option(c("-a", "--maf_sample_id_colname"), type = "character", help = "Column name in MAF file for sample IDs"),
    make_option(c("-b", "--metadata_sample_id_colname"), type = "character", help = "Column name in metadata file for sample IDs"),
    make_option(c("-s", "--sv_from_metadata"), type = "character", help = "A named vector that specifies the columns containing the oncogene translocation status for any SV that is annotated in the metadata"),
    make_option(c("-z", "--translocation_status"), type = "character", help = "A named vector that specifies the values corresponding to positive translocation, negative translocation, and NA for any SV that is annotated in the metadata"),
    make_option(c("-c", "--truth_column_colname"), type = "character", help = "Column name in metadata file with truth labels")
)

opt <- parse_args(OptionParser(option_list=option_list))

metadata <- opt$metadata
metadata_sample_id_colname <- opt$metadata_sample_id_colname

maf <- opt$maf
maf_sample_id_colname <- opt$maf_sample_id_colname

sv_from_metadata <- strsplit(opt$sv_from_metadata, ",")[[1]] 

trans_raw <- eval(parse(text = opt$translocation_status))
# Convert the special placeholder to actual NA
if ("NA_STATUS" %in% names(trans_raw)) {
    names(trans_raw)[names(trans_raw) == "NA_STATUS"] <- NA
}
# Build reverse map: e.g. "yes" → "POS", "no" → "NEG", "Not Available" → NA
value_map <- setNames(names(trans_raw), unlist(trans_raw))


truth_column_colname <- opt$truth_column_colname

output_meta_dir <- opt$output_meta_dir
output_maf_dir <- opt$output_maf_dir

test_metadata <- readr::read_tsv(
        metadata
    ) %>%
    rename(
        sample_id = metadata_sample_id_colname,
        lymphgen = truth_column_colname
    ) %>%
    mutate(across( # Assign POS, NEG, NA values
        all_of(sv_from_metadata),   
        ~ dplyr::recode(.x, !!!value_map, .default = .x)  
    ))

test_maf <- readr::read_tsv(
        maf
    ) %>%  
    mutate(
        sample_id = .data[[maf_sample_id_colname]]
    )

write_tsv(test_metadata, output_meta_dir)
write_tsv(test_maf, output_maf_dir)

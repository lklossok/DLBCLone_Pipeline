# Load necessary libraries
suppressWarnings(suppressMessages(library(optparse)))
suppressWarnings(suppressMessages(library(GAMBLR.open)))
suppressWarnings(suppressMessages(library(tidyverse)))

# Define command line options
option_list <- list(
    make_option(c("-m", "--form_metadata"), type = "character", help = "path to formatted test metadata file"),
    make_option(c("-f", "--form_maf"), type = "character", help = "path to formatted test maf file"),
    make_option(c("-o", "--output_matrix_dir"), type = "character", help = "output mutation matrix path"),
    make_option(c("-s", "--sv_from_metadata"), type = "character", help = "A named vector that specifies the columns containing the oncogene translocation status for any SV that is annotated in the metadata"),

    make_option(c("-g", "--genes"), type = "character", default="c('ACTB','ACTG1','BCL10','BCL2','BCL2L1','BCL6','BIRC3','BRAF','BTG1','BTG2','BTK','CD19','CD70','CD79B','CD83','CDKN2A','CREBBP','DDX3X','DTX1','DUSP2','EDRF1','EIF4A2','EP300','ETS1','ETV6','EZH2','FAS','FCGR2B','FOXC1','FOXO1','GNA13','GRHPR','HLA-A','HLA-B','HNRNPD','IRF4','IRF8','ITPKB','JUNB','KLF2','KLHL14','KLHL6','KMT2D','MEF2B','MPEG1','MYD88','NFKBIA','NFKBIE','NFKBIZ','NOL9','NOTCH1','NOTCH2','OSBPL10','PIM1','PIM2','PRDM1','PRDM15','PRKDC','PRRC2C','PTPN1','RFTN1','S1PR2','SETD1B','SGK1','SOCS1','SPEN','STAT3','STAT6','TBCC','TBL1XR1','TET2','TMEM30A','TMSB4X','TNFAIP3','TNFRSF14','TOX','TP53','TP73','UBE2A','WEE1','XBP1','ZFP36L1','MYD88HOTSPOT','BCL2_SV','BCL6_SV')",help = "Vector of gene symbols to include"),
    make_option(c("-y", "--synon_genes"), type = "character", default=NULL, help = "Vector of gene symbols for synonymous mutations"),
    make_option(c("-z", "--hotspot_genes"), type="character", default="MYD88", help = "Vector specifying genes for which hotspot mutations should be separately annotated"),
    make_option(c("-b", "--genome_build"), type = "character", default="grch37", help = "genome build default: grch37"),
    make_option(c("-j", "--sv_value"), type = "integer", default=2, help = "structural variant value default: 2"),
    make_option(c("-k", "--synon_value"), type = "integer", default=1, help = "synonymous mutation value default: 1"),
    make_option(c("-l", "--coding_value"), type = "integer", default=2, help = "coding mutation value default: 2"),
    make_option(c("-i", "--include_ashm"), type = "logical", default=FALSE, help = "include ASHM features default: FALSE"),
    make_option(c("-a", "--annotated_sv"), type = "character", default=NULL, help = "path to annotated sv file"),
    make_option(c("-u", "--include_GAMBL_sv"), type = "logical", default=FALSE, help = "include GAMBLR sv features default: FALSE")
)

opt <- parse_args(OptionParser(option_list=option_list))

metadata <- opt$form_metadata
maf <- opt$form_maf
output_matrix_dir <- opt$output_matrix_dir
sv_from_metadata <- strsplit(opt$sv_from_metadata, ",")[[1]]

genes <- eval(parse(text=opt$genes))

if (is.null(opt$synon_genes)) {
    synon_genes <- NULL
} else {
    synon_genes <- strsplit(opt$synon_genes, ",")[[1]]
    synon_genes <- trimws(synon_genes)
}

if (is.null(opt$hotspot_genes)) {
    hotspot_genes <- NULL
} else {
    hotspot_genes <- strsplit(opt$hotspot_genes, ",")[[1]]
    hotspot_genes <- trimws(hotspot_genes)
}

genome_build <- opt$genome_build
sv_value <- opt$sv_value
synon_value <- opt$synon_value
coding_value <- opt$coding_value
include_ashm <- opt$include_ashm
annotated_sv <- opt$annotated_sv
include_GAMBL_sv <- opt$include_GAMBL_sv


test_metadata <- readr::read_tsv(metadata) 

maf <- readr::read_tsv(maf)

#convert to S3 object for compatability
test_maf <- create_maf_data(
        maf_df = maf,
        genome_build = genome_build
    ) %>%
    mutate(Chromosome = as.character(Chromosome))

test_matrix = assemble_genetic_features(
    these_samples_metadata = test_metadata,
    maf_with_synon = test_maf,
    sv_from_metadata = sv_from_metadata,
    genes = genes,
    synon_genes = synon_genes,
    hotspot_genes = hotspot_genes,
    genome_build = genome_build,
    sv_value = sv_value,
    synon_value = synon_value,
    coding_value = coding_value,
    include_ashm = include_ashm,
    annotated_sv = annotated_sv,
    include_GAMBL_sv = include_GAMBL_sv
    ) %>%
    rownames_to_column(var = "sample_id")

write_tsv(test_matrix, output_matrix_dir)

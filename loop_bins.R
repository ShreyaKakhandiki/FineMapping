
require(reticulate)
library(coloc)
library(stringr) 

target_chr <- commandArgs(trailingOnly = T)[1] # instead of manually specified "chr3"

# alternative: target_chr <- 9

chr_str <- paste0("chr",target_chr,"_") 


# Load the Python function
py_run_string("
import numpy as np
import pandas as pd
import scipy.sparse as sparse

def load_ld_npz(ld_prefix):

    passed = 'true'

    # Base directory for the files
    
    # DON'T USE SHORTCUT
    base_dir = '/u/project/cluo_scratch/shreyaka/broad-alkesgroup-ukbb-ld/UKBB_LD/'
    gz_file = base_dir + '%s.gz' % ld_prefix

    # Load the gz file into a dataframe
    try:
        df_ld_snps = pd.read_csv(gz_file, sep='\\s+')
    except:
        print('Error with' , gz_file) 
        passed = 'false'
        df_ld_snps = pd.read_csv(gz_file, sep='\\s+', error_bad_lines='skip')
  
    # Rename columns for consistency
    rename_dict = {
        'rsid': 'SNP',
        'chromosome': 'CHR',
        'position': 'BP',
        'allele1': 'A1',
        'allele2': 'A2'
    }
    df_ld_snps.rename(columns=rename_dict, inplace=True, errors='ignore')

    # Ensure necessary columns exist
    required_columns = ['SNP', 'CHR', 'BP', 'A1', 'A2']
    for col in required_columns:
        if col not in df_ld_snps.columns:
            raise ValueError(f'Missing column: {col}')

    # Create a unique index for the dataframe
    df_ld_snps.index = df_ld_snps['CHR'].astype(str) + '.' + df_ld_snps['BP'].astype(str)

    # Load the LD matrix
    npz_file = base_dir + f'/{ld_prefix}.npz'
    try:
        R = sparse.load_npz(npz_file).toarray()
        R += R.T
    except ValueError:
        raise IOError(f'Corrupt file: {npz_file}')

    # Construct and return the resultant dataframe
    df_R = pd.DataFrame(R, index=df_ld_snps.index, columns=df_ld_snps.index)

    return df_R, df_ld_snps, passed
")



# extract a list of files that we want to open 

base_dir <- "/u/project/cluo_scratch/shreyaka/broad-alkesgroup-ukbb-ld/UKBB_LD/"

gz_files <- list.files(path = base_dir , pattern = "\\.gz$", full.names = FALSE)

ld_prefix_list <- sub("\\.gz$", "", gz_files) 

# Remove the first two files

skip_files <- 2

ld_prefix_list <- ld_prefix_list[-seq_len(skip_files)]


ld_prefix_list <- ld_prefix_list[grep(chr_str, ld_prefix_list)]

ld_prefix_list_sub <- ld_prefix_list[48:length(ld_prefix_list)]

# ld_prefix_list_subset <- ld_prefix_list[1:2]


# Read in sumstats file
sumstats <- read.table(file = "/u/home/s/shreyaka/project-cluo/bed_files/PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.tsv", sep = '\t', header = TRUE)


sumstats <- sumstats[complete.cases(sumstats), ] # Remove rows with NA values

## test the specific one 

summary_list <- list()

for (ld_prefix in ld_prefix_list_sub) {
    
    print(paste0("currently running on ", ld_prefix))
    
     
    # Load LD data for the current prefix
    result <- py$load_ld_npz(ld_prefix)
    df_R <- result[[1]]
    df_ld_snps <- result[[2]]
    
    
     snp_names <- df_ld_snps$SNP
  
     # Make SNP names unique by adding a unique identifier
      unique_snp_names <- make.unique(snp_names)
  
     # Set row and column names of df_R using unique SNP names
      rownames(df_R) <- unique_snp_names
      colnames(df_R) <- unique_snp_names

    
    # Extract chromosome number, start BP, and end BP from ld_prefix
    ld_parts <- unlist(strsplit(ld_prefix, "_"))
    chromosome <- as.numeric(str_extract(ld_prefix, "\\d+"))
    start_bp <- as.numeric(ld_parts[2])
    end_bp <- as.numeric(ld_parts[3])
    
    # Extract rows from sumstats file within the LD prefix range
    matching_sumstats <- subset(sumstats, CHROM == chromosome & POS >= start_bp & POS <= end_bp)
    
    # Update LD matrix and sumstats data to include only common SNPs
    common_snps <- intersect(matching_sumstats$ID, df_ld_snps$SNP)
    subset_df_R <- df_R[common_snps, common_snps]
    sumstats_ordered <- matching_sumstats[matching_sumstats$ID %in% common_snps, ]
    
    
    # reorganize the sumstats file to be in the same order as the matrix 

    ordered_snps <- rownames(subset_df_R)

    # Order the sumstats dataframe based on the SNP names
    sumstats_ordered <- matching_sumstats[matching_sumstats$ID %in% ordered_snps, ]

    
        # Format sumstats data for coloc.SuSIE
        sumstats_ordered$MAF <- pmin(ifelse(sumstats_ordered$A1 == 'A', sumstats_ordered$FCAS, sumstats_ordered$FCON), 1 - ifelse(sumstats_ordered$A1 == 'A', sumstats_ordered$FCAS, sumstats_ordered$FCON))
    
    
    coloc_data <- list(
        beta = sumstats_ordered$BETA,
        varbeta = sumstats_ordered$SE^2,
        type = "quant",
        position = sumstats_ordered$POS, 
        snp = sumstats_ordered$ID, 
        se = sumstats_ordered$SE, 
        Neff = sumstats_ordered$NEFFDIV2, 
        MAF = sumstats_ordered$MAF, 
        N = nrow(subset_df_R),
        LD = as.matrix(subset_df_R)
    )
    
    
# Check if the dataset is valid
tryCatch({
  dataset_check <- check_dataset(coloc_data, req = "LD")
  dataset_bool <- is.null(dataset_check)
}, error = function(e) {
  # Handle the error (print a message or take appropriate action)
  print("Error checking dataset:", e$message)
  dataset_bool <- FALSE  # Assume dataset is invalid
})

print(dataset_check) 

n <- nrow(subset_df_R)

# If the dataset is valid (returns NULL) and there are more than 1 rows in subset_df_R, proceed with runsusie
    
tryCatch({
  if (dataset_bool && n > 1) {
    # Use tryCatch to handle errors during Susie execution
    tryCatch({
      S <- runsusie(coloc_data)
      summary_list[[ld_prefix]] <- S  
    }, error = function(e) {
      print(paste("Error running Susie for", ld_prefix, ":"))
      print(e)
    })
  } else {
    # If the dataset is not valid or there are not enough rows in subset_df_R
    print("Error: Dataset is not valid or does not have enough rows.")
  }
}, error = function(e) {
  print("Error in if statement:")
  print(e)
})

}


#save the list of summary objects

saveRDS(summary_list, paste0("/u/home/s/shreyaka/project-cluo/Luo_Lab_Initial2/reticulate/summary/", chr_str, "susie_objects.rds"))



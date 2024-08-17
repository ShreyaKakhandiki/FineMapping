# now that we have basically fully processed one matrix, write a loop to process all matrices and have the results stored in a big R object (at the end)

# keep a counter of the matching SNPs 

# Load reticulate for Python integration
require(reticulate)

print("the code has started 0") 

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

print("function has been defined 1") 

# CODE BLOCK TO EXTRACT THE FILES THAT WE WANT TO OPEN


base_dir <- "/u/project/cluo_scratch/shreyaka/broad-alkesgroup-ukbb-ld/UKBB_LD/"
gz_files <- list.files(path = base_dir , pattern = "\\.gz$", full.names = FALSE)
ld_prefix_list <- sub("\\.gz$", "", gz_files) # remove .gz to get the prefix

skip_files <- 2

# Remove the first two files
ld_prefix_list <- ld_prefix_list[-seq_len(skip_files)]

print("LD prefixes list created 2")


# Initialize an empty list to store temporary matrices
temp_matrices <- list()

skip_count <- 0
success_count <- 0 

# Loop through each LD prefix
for (ld_prefix in ld_prefix_list) {
  # Call the Python function to load LD matrix and SNP info
  result <- py$load_ld_npz(ld_prefix)
  
  # Extract LD matrix and SNP info
  df_R <- result[[1]]
  df_ld_snps <- result[[2]]
  passed <- result[[3]]
  
  if (passed == 'false') {
    skip_count <- skip_count + 1
    next
}
  success_count <- success_count + 1
    
  # Assuming the SNP names are in a column named "SNP" in df_ld_snps dataframe
  snp_names <- df_ld_snps$SNP
  
  # Make SNP names unique by adding a unique identifier
  unique_snp_names <- make.unique(snp_names)
  
  # Set row and column names of df_R using unique SNP names
  rownames(df_R) <- unique_snp_names
  colnames(df_R) <- unique_snp_names
  
  # Read the sumstats file
  sumstats <- read.table(file = "/u/home/s/shreyaka/project-cluo/bed_files/PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.tsv",    sep = '\t', header = TRUE)
  
  # Extract SNP columns as vectors
  sumstats_snps <- sumstats$ID
  
  # Find rows in df_R where SNP names match those in sumstats_snps
  common_snps <- intersect(rownames(df_R), sumstats_snps)
    
 # Extract rows and columns based on common SNP names
  subset_df_R <- df_R[common_snps, common_snps]
  
  # Store the matching rows in a temporary matrix 
  temp_matrices[[ld_prefix]] <- as.matrix(subset_df_R)
    
  print("matrix processing 3")   
    
}

reticulate::py_last_error()

print(paste0("skip_count: ", skip_count))

print(paste0("success_count: ", success_count))


# Combine all temporary matrices into one big matrix
all_results <- do.call(rbind, temp_matrices)

# Write the combined matrix to an RDS object
save(all_results, file = "/u/home/s/shreyaka/project-cluo/Luo_Lab_Initial2/reticulate/all_results.RData")
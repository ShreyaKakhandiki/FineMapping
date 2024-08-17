# Read in sumstats file
    sumstats <- read.table(file = "/u/home/s/shreyaka/project-cluo/bed_files/PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.tsv", sep = '\t', header = TRUE)


    sumstats <- sumstats[complete.cases(sumstats), ] # Remove rows with NA values
    
    colnames(sumstats)[colnames(sumstats) == "ID"] ="SNP"
    

# directory
directory <- "/u/home/s/shreyaka/project-cluo/Luo_Lab_Initial2/reticulate/summary"




# Initialize an empty list to store results from each file
all_results <- list()

# List all RDS files in the directory
rds_files <- list.files(path = directory, pattern = "*.rds", full.names = TRUE)

rds_files_test <- rds_files[1:3]

# Loop through each RDS file
for (rds_file in rds_files_test) {
    
  susie_list <- readRDS(rds_file) 
    
    
  # Load the RDS file
  print(rds_file)
    
   
   # Get column names from sumstats
    column_names <- colnames(sumstats)
    
    
    # Create an empty dataframe with the same column names as sumstats
    all_causal <- data.frame(matrix(ncol = length(column_names), nrow = 0))
    colnames(all_causal) <- column_names
     
    
    # Loop through each Susie object in the list
    for (susie_obj in susie_list) {
        # Extract causal SNPs and add them to the list
        # Initialize an empty vector to store rsids
        rsids <- c()

        # Extract causal SNPs list from susie_obj
        causal_snps_list <- susie_obj$sets$cs

        
        # Iterate over all elements in causal_snps_list
        for (set_list in causal_snps_list) {
    
        rsids <- c(rsids, names(set_list)) 
    
        }

        # Combine all causal SNPs into a single vector
        causal_snps <- rsids  

        # Create an empty dataframe with the same column names as sumstats
        causal_sumstats <- data.frame(matrix(ncol = length(column_names), nrow = 0))
        colnames(causal_sumstats) <- column_names
    
        
        causal_sumstats <- sumstats[sumstats$SNP %in% causal_snps, ]
        
        library(dplyr)

    # If all_causal is not already initialized as a data frame, you can do so
    all_causal <- data.frame()

    # Combine data frames while preserving column names
    all_causal <- bind_rows(all_causal, causal_sumstats)
    
    head(all_causal)
    }
      

   all_results[[basename(rds_file)]] <- all_causal
  
}





# Initialize an empty list to store combined data frames
combined_dfs <- list()

# Iterate over each susie object and combine its data frame into the list
for (susie_obj_file in all_results) {
  susie_obj <- readRDS(susie_obj_file)
  combined_dfs[[susie_obj_file]] <- susie_obj
}

# Combine all data frames into a single data frame
giant_df <- do.call(rbind, combined_dfs)

                                  
head(giant_df)
                                  
                                  
# Save the giant dataframe to an RDS file
saveRDS(giant_df, "/u/home/s/shreyaka/project-cluo/Luo_Lab_Initial2/reticulate/giant_susie_df.rds")
                                  
print("RDS object saved to file")
                                  
                                  
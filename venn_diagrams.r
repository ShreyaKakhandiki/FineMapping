credible_snps <- readRDS("/u/home/s/shreyaka/project-cluo/Luo_Lab_Initial2/reticulate/giant_susie_df.rds")
#credible_snps <- t(credible_snps)

credible_snps

# credible_snps


write.csv(credible_snps, file = "/u/home/s/shreyaka/project-cluo/Luo_Lab_Initial2/reticulate/scz_causal_snps")

 

paper_indeces <- read.table("/u/home/s/shreyaka/project-cluo/Luo_Lab_Initial2/reticulate/sup_table_1")
# Create a new data frame with unique IDs
# unique_index_df <- data.frame(unique_ids = unique(paper_indeces$ids))

paper_indeces <- as.data.frame(paper_indeces)

col_names <- as.character(paper_indeces[1, ])

# Remove the first row (which contains the column names)
paper_indeces <- paper_indeces[-1, ]

# Set the column names
colnames(paper_indeces) <- col_names

colnames(paper_indeces)[1] <- "SNP"

nrow(paper_indeces)

# Keep only distinct rows based on the "ids" column
# unique_ids_df <- distinct(paper_indeces, ids)

# colnames(paper_indeces) <- c("Index SNP")

unique_ids_df <- paper_indeces[!duplicated(paper_indeces$SNP, fromLast = TRUE), , drop = FALSE]


nrow(unique_ids_df)

#INDEX SNPs VENN DIAGRAM

# Install and load the VennDiagram package
library(VennDiagram)

# Create vectors of the credible SNPs from the paper results and your results
paper_snps <- paper_indeces$SNP # Replace with your actual rsids from the paper results
your_snps <-  credible_snps$SNP # Replace with your actual rsids from your results

# Create a list of sets for the Venn diagram
venn_list <- list(Paper_Results = paper_snps, Your_Results = your_snps)

# Create the Venn diagram
venn_result <- venn.diagram(
  x = venn_list,
  category.names = c("Paper Result", "Our Results"),
  filename = NULL
)

# Plot the Venn diagram
grid.draw(venn_result)

# check table 11 to see if inex snps match table 1


table_11 <- read.csv("/u/project/cluo/shreyaka/Luo_Lab_Initial2/reticulate/supp_table_11.csv", header = TRUE)

index_snps <- table_11$index_snp

unique_ids <- index_snps[!duplicated(index_snps, fromLast = TRUE), drop = FALSE]

table_11

# make a venn diagram that shows overlap between tables 1 and 11

# Install and load the VennDiagram package
library(VennDiagram)

# Create vectors of the credible SNPs from the paper results and your results
paper_snps <- paper_indeces$SNP # Replace with your actual rsids from the paper results


# Create a list of sets for the Venn diagram
venn_list <- list(Table_1 = paper_snps, Table_11 = unique_ids)

# Create the Venn diagram
venn_result <- venn.diagram(
  x = venn_list,
  category.names = c("Table 1", "Table 11"),
  filename = NULL
)

# Plot the Venn diagram
grid.draw(venn_result)

# create venn diagram of our causal snps vs all of their snps 


# Install and load the VennDiagram package
library(VennDiagram)

# Create vectors of the credible SNPs from the paper results and your results
paper_snps <- table_11$rsid # Replace with your actual rsids from the paper results
your_snps <-  credible_snps$SNP # Replace with your actual rsids from your results

# Create a list of sets for the Venn diagram
venn_list <- list(Paper_Results = paper_snps, Your_Results = your_snps)

# Create the Venn diagram
venn_result <- venn.diagram(
  x = venn_list,
  category.names = c("Paper Results", "Our Results"),
  filename = NULL
)

# Plot the Venn diagram
grid.draw(venn_result)



# create venn diagram of our causal snps vs their index snps


# Install and load the VennDiagram package
library(VennDiagram)

# Create vectors of the credible SNPs from the paper results and your results
paper_snps <- unique_ids # Replace with your actual rsids from the paper results
your_snps <-  credible_snps$SNP # Replace with your actual rsids from your results

# Create a list of sets for the Venn diagram
venn_list <- list(Paper_Results = paper_snps, Your_Results = your_snps)

# Create the Venn diagram
venn_result <- venn.diagram(
  x = venn_list,
  category.names = c("Paper Results", "Our Results"),
  filename = NULL
)

# Plot the Venn diagram
grid.draw(venn_result)


# create a plot that shows a bar for each index snp with the height of the bar being % correspondence with our results

# Count the occurrences of each index SNP
index_snp_counts <- table(index_snps)

# Convert the table to a data frame
index_snp_counts_df <- as.data.frame(index_snp_counts)

# Rename the columns for clarity
colnames(index_snp_counts_df) <- c("index_snp", "count")

# Print the data frame
print(index_snp_counts_df)






# Calculate the percentage of SNPs in credible_snps that match each index SNP
index_snp_percentages <- sapply(unique(table_11$index_snp), function(idx_snp) {
  idx_snp_snps <- table_11[table_11$index_snp == idx_snp, "rsid"]
  match_snps <- sum(credible_snps$SNP %in% idx_snp_snps)
  total_snps <- index_snp_counts_df[index_snp_counts_df$index_snp == idx_snp, "count"]
  percentage <- (match_snps/total_snps) * 100
    
  if(percentage != 0)
      {
      print(index_snp_counts[idx_snp])
  }
  return(percentage)
})

# Create a data frame from the calculated percentages
index_snp_percentages_df <- data.frame(index_snp = names(index_snp_percentages),
                                       percentage = as.numeric(index_snp_percentages))

# Create a bar plot using ggplot2
library(ggplot2)
ggplot(index_snp_percentages_df, aes(x = index_snp, y = percentage)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = "Percentage of Matching SNPs",
       x = "Index SNP",
       y = "Percentage") +
  theme(axis.text.x = element_blank())

# check the index_snps out in the sumstats file to see if they are very significant 

sumstats <- read.table(file = "/u/home/s/shreyaka/project-cluo/bed_files/PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.tsv", sep = '\t', header = TRUE)

    sumstats <- sumstats[complete.cases(sumstats), ] # Remove rows with NA values
    
    colnames(sumstats)[colnames(sumstats) == "ID"] ="SNP"

# Look at the GWAS catalog and see how many of our causal snps overlap 

credible_snps <- readRDS("/u/home/s/shreyaka/project-cluo/Luo_Lab_Initial2/reticulate/giant_susie_df.rds")


gwas_catalog <- read.table(file = "/u/home/s/shreyaka/project-cluo/Luo_Lab_Initial2/gwascatalog.txt", sep = '\t', header = TRUE) 

# head(gwas_catalog)



catalog_overlap <- gwas_catalog[gwas_catalog$IndSigSNP %in% credible_snps$SNP,]

unique(catalog_overlap$Trait)


# bar plot the top 10 categories and # of times they occur in this list 

# Step 1: Calculate frequency of each unique trait
trait_counts <- table(catalog_overlap$Trait)



# Step 2: Select top 10 most frequent traits
top_10_traits <- names(sort(trait_counts, decreasing = TRUE))[1:10]

top_10_traits

# Step 3: Create bar plot
barplot(trait_counts[top_10_traits], 
        names=top_10_traits,
        main = "Top 10 Traits", 
        xlab = "Trait", 
        ylab = "Frequency",
        col = rainbow(10),
       las = 2)  # You can customize the colors as you like


# make a new bar plot with the filtered/combined values 

# also remove any rows with "Trubetskoy V" as the first author 

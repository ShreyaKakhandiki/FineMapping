# FineMapping
Code to identify potential causal variants from GWAS summary statistics and Linkage Disequilibrium information.


**Overview**
This project focuses on fine-mapping to determine causal variants associated with common neurological disorders. Fine-mapping aims to identify the genetic variants most likely to influence the risk of developing such disorders by refining the association signals detected in genome-wide association studies (GWAS). This repository contains the code and documentation for performing fine-mapping analyses, visualizing results, and interpreting the data.

**Project Structure**
The repository is organized as follows:

data/: Contains example data files and datasets used for analysis.
scripts/: Contains R scripts for performing fine-mapping, visualizing results, and generating plots.
results/: Stores output files, including fine-mapped variant lists and visualization plots.
docs/: Includes documentation and instructions for using the code.


**Installation**
To use the code in this repository, you'll need to have R installed along with several R packages. You can install the required packages by running the following commands in your R environment:

r
Copy code
install.packages(c("data.table", "ggplot2", "cowplot", "dplyr", "tidyr", "ggrepel"))


**Usage
Data Preparation**
Ensure that your data files are in the data/ directory. The expected format for input data is:

GWAS Summary Statistics: A table with columns for SNP, effect size, standard error, p-value, etc.
LD Matrix: A matrix of pairwise linkage disequilibrium (LD) values between SNPs.

Running Fine-Mapping
To perform fine-mapping, use the provided R scripts located in the scripts/ directory. For example:

r
Copy code
source("scripts/run_finemapping.R")
This script will process the GWAS summary statistics, perform fine-mapping analysis, and save the results in the results/ directory.

Visualizing Results
You can generate plots to visualize the fine-mapping results using the following command:

r
Copy code
source("scripts/plot_results.R")
This script will produce plots showing the fine-mapping results, including credible sets and association signals, and save them to the results/ directory.

Example
To run a full analysis pipeline, you can execute the following commands in R:

r
Copy code
# Load the required libraries
library(data.table)
library(ggplot2)

# Run the fine-mapping analysis
source("scripts/run_finemapping.R")

# Plot the results
source("scripts/plot_results.R")


**Documentation**
Additional documentation and user instructions can be found in the docs/ directory. This includes detailed explanations of the analysis methods, parameters, and interpretation of results.

**Contributing**
Contributions to this project are welcome. If you have suggestions, improvements, or bug fixes, please submit a pull request or open an issue.

**License**
This project is licensed under the MIT License. See the LICENSE file for details.

**Contact**
For questions or additional information, please contact Shreya Kakhandiki at shreya.kakhandiki@gmail.com.

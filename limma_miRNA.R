# Set working directory
setwd("D:/project")

install.packages("BiocManager")
install.packages("forcats")
install.packages("stringr")
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("readr")
install.packages("tidyr")
install.packages("survminer")
BiocManager::install("GEOquery")
BiocManager::install("limma")
BiocManager::install("pheatmap")
BiocManager::install("org.Hs.eg.db")

library(GEOquery)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(readr)
library(limma)

# Define a function for analysis
analyze_dataset <- function(gse, my_id, output_file) {
  cat("Processing dataset:", my_id, "\n")
  
  # Shift expression data to make them non-negative
  exprs(gse) <- exprs(gse) + abs(min(exprs(gse)))
  
  # Boxplots
  boxplot(exprs(gse), outline = FALSE)
  
  # Sample information
  sampleInfo <- pData(gse)
  sampleInfo <- sampleInfo %>%
    select(source_name_ch1, characteristics_ch1.1) %>%
    rename(group = source_name_ch1, patient = characteristics_ch1.1)
  
  # Correlation heatmap
  corMatrix <- cor(exprs(gse), use = "c")
  pheatmap(corMatrix)
  
  # PCA
  pca <- prcomp(t(exprs(gse)))
  plot_data <- cbind(sampleInfo, pca$x)
  ggplot(plot_data, aes(x = PC1, y = PC2, col = group, label = paste("Patient", patient))) +
    geom_point() + geom_text_repel()
  
  # Feature selection and output
  features <- fData(gse)
  features <- select(features, mirna_ids)
  full_output <- cbind(features, exprs(gse))
  write_csv(full_output, path = output_file)
  
  # Limma analysis
  design <- model.matrix(~0 + sampleInfo$group)
  colnames(design) <- c("Normal", "Tumour")
  
  is_expressed <- exprs(gse) > median(exprs(gse))
  keep <- rowSums(is_expressed) > 2
  gse <- gse[keep,]
  
  fit <- lmFit(exprs(gse), design)
  contrasts <- makeContrasts(Tumour - Normal, levels = design)
  fit2 <- contrasts.fit(fit, contrasts)
  fit2 <- eBayes(fit2)
  topTable(fit2)
}


# Dataset IDs
dataset_ids <- c("GSE45666", "GSE70754", "GSE38167", "GSE41922")

# Apply the function to each dataset
for (my_id in dataset_ids) {
  gse <- getGEO(my_id)[[1]]
  output_file <- paste0("gse_", my_id, "_output.csv")
  analyze_dataset(gse, my_id, output_file)
}

library(readxl)
library(tidyverse)
library(BiocManager)
library(miRBaseConverter)
library(openxlsx)

setwd("E:/wes data")

# Function to process each dataset
process_dataset <- function(file_path, sheet_number, target_version, output_file) {
  sheet_name <- paste("Mirnas found", sheet_number)
  
  cat("Processing dataset:", sheet_name, "\n")
  
  # Read data from Excel file
  dataset <- read_excel(file_path, sheet = sheet_number)  # Assuming sheet number is used as a name
  
  # Extract miRNA IDs
  miRNA_ids <- dataset %>% pull(miRNA_ID)
  
  # Check miRNA version
  miRNA_version <- checkMiRNAVersion(miRNA_ids, verbose = TRUE)
  # Convert miRNA names to Accession numbers
  result <- miRNA_NameToAccession(miRNA_ids, version = miRNA_version)
  
  # Convert Accession numbers to names in the target version
  result_names <- miRNA_AccessionToName(result$Accession, targetVersion = target_version)
  
  # Combine miRNA names and Accession numbers
  complete_result <- cbind(miRNA_ids, result_names)
  
  # Create a new Excel workbook
  wb <- createWorkbook()
  
  # Add data to a new sheet in the workbook
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet = sheet_name, x = complete_result)
  
  # Save the workbook
  saveWorkbook(wb, output_file, overwrite = TRUE)
  
  return(complete_result)
}

# File paths and sheet numbers for the specified datasets
file_paths <- c(
  "E:/wes data/Unilateral/GSE45666_Annotation.xlsx",
  "E:/wes data/Unilateral/GSE70754_Annotation.xlsx",
  "E:/wes data/Unilateral/GSE38167_Annotation.xlsx",
  "E:/wes data/Unilateral/GSE41922_Annotation.xlsx"
)

sheet_numbers <- c(1, 2, 3, 4)

# Target version for conversion
target_version <- "v15"

# Output file names
output_files <- paste0("C:/Users/vinni/OneDrive/Desktop/DEGS/Complete_results_Mirnas_found_", sheet_numbers, ".xlsx")

# Process each dataset
for (i in seq_along(file_paths)) {
  process_dataset(file_paths[i], sheet_numbers[i], target_version, output_files[i])
}

# Combine miRNA names for comparison
miRNA_names_combined <- process_dataset(file_path = "E:/wes data/Unilateral/GSE45666_Annotation.xlsx",
                                        sheet_number = 1,  # Replace with actual sheet number
                                        target_version = "v15",
                                        output_file = "C:/Users/vinni/OneDrive/Desktop/DEGS/Complete_results_Mirnas_found_1.xlsx")

# Process each dataset
result_datasets <- lapply(seq_along(file_paths), function(i) {
  process_dataset(file_paths[i], sheet_numbers[i], target_version, output_files[i])
})

# Extract miRNA names from the first dataset
miRNA_names_combined <- result_datasets[[1]] %>% select("TargetName")

# Find common miRNAs with the rest of the datasets
common_miRNAs <- Reduce(intersect, lapply(result_datasets[-1], function(dataset) {
  dataset %>% select("TargetName")
}))

# Load required libraries
library(DESeq2)

# Define DESeq2 function
DESeq2 <- function(data_folder, metadata_folder) {
  
  # List all data files in the "data" and "metadata" folder
  data_files <- list.files(data_folder, full.names = TRUE)
  metadata_files <- list.files(metadata_folder, full.names = TRUE)
  
  # Check if number of data files and metadata files match
  if (length(data_files) != length(metadata_files)) {
    stop("Number of data files and metadata files do not match.")
  }
  
  # Initialize progress bar
  pb <- txtProgressBar(min = 0, max = length(data_files), style = 3)
  
  # Loop through indices of data files
  for (i in seq_along(data_files)) {
    # Get data file and metadata file paths
    data_file <- data_files[i]
    metadata_file <- metadata_files[i]
    
    # Read data
    raw_counts <- read.table(data_file, header = TRUE)
    row.names(raw_counts) <- raw_counts[,1]
    raw_counts <- raw_counts[,-1]
    raw_counts <- as.matrix(raw_counts)
    data <- raw_counts[rowSums(raw_counts) >= 0, ]
    
    # Read metadata
    meta <- read.table(metadata_file, header = TRUE)
    rownames(meta) <- meta$labels
    meta$group <- factor(meta$group, levels = unique(meta$group))
    
    # Check compatibility
    if (!all(colnames(data) %in% rownames(meta)) | !all(colnames(data) == rownames(meta))) {
      warning(paste("Data and metadata are not compatible for", basename(data_file), ". Skipping."))
      next
    }
    
    # Perform DESeq2 analysis
    dds <- DESeqDataSetFromMatrix(countData = round(data), colData = meta, design = ~ group)
    dds <- estimateSizeFactors(dds) 
    dds <- DESeq(dds)
    res <- results(dds)
    
    # Write results to CSV
    output_file <- paste("DE_", levels(meta$group)[2], "_vs_", levels(meta$group)[1], ".csv", sep="")
    write.csv(res, file = output_file)
  
    # Update progress bar
    setTxtProgressBar(pb, i)
  }
  # Close progress bar
  close(pb)
}

# Extract data and metadata folder paths from command line arguments
data_folder <- commandArgs(trailingOnly = TRUE)[1]
metadata_folder <- commandArgs(trailingOnly = TRUE)[2]

# Run DESeq2 function
DESeq2(data_folder, metadata_folder)


# Load necessary libraries
library(data.table)
library(ggplot2)
library(dplyr)
library(UpSetR)
library(tidyr)

# Set the path to the featureCounts directory
edgeR_featureCounts_dir <- "edgeR_featureCounts"

# Function to list files and extract DEGs (genes) into columns by group_time
list_files_DEGs <- function(dir_path) {
  degs_list <- list()  # Initialize an empty list
  # List directories starting with 'Results_'
  subdirs <- list.dirs(dir_path, recursive = FALSE, full.names = TRUE)
  subdirs <- subdirs[grepl("/Results_", subdirs)]
  
  for (subdir in subdirs) {
    # List CSV files in the subdirectory
    files <- list.files(subdir, pattern = "\\.csv$", full.names = TRUE)
    
    for (file_path in files) {
      # Print file path to check which files are being processed
      print(paste("Processing file:", file_path))
      
      # Read the file to get DEGs
      data <- fread(file_path)
      
      # Extract the first column (gene names)
      genes <- data[[1]]
      
      # Extract sample name from the file path
      sample_name <- gsub("Results_", "", basename(dirname(file_path)))  # Remove "Results_" from the directory name
      
      # Add the genes to the list for this group_time (sample_name)
      degs_list[[sample_name]] <- genes  
    }
  }
  
  # Convert the list of genes into a wide-format data frame
  max_length <- max(sapply(degs_list, length))  # Get the longest list of genes
  degs_df <- data.frame(lapply(degs_list, function(x) c(x, rep(NA, max_length - length(x)))))  # Fill shorter columns with NAs
  
  # Assign group_time (sample_name) as column names
  colnames(degs_df) <- names(degs_list)
  
  return(degs_df)
}

# Call the function to extract DEGs into columns by group_time
degs_featureCounts <- list_files_DEGs(edgeR_featureCounts_dir)

# View the resulting data frame with DEGs for each group_time in separate columns
print(degs_featureCounts)

# Save the DEGs data frame to a text file
write.table(degs_featureCounts, "DEGs_by_group_time.txt", sep = "\t", quote = FALSE, row.names = FALSE)


# Create a binary matrix for upset plot where 1 indicates presence of DEG in that sample
degs_binary_matrix <- lapply(degs_featureCounts, function(genes) {
  # Remove NA values from each list of genes before creating the binary matrix
  genes <- genes[!is.na(genes)]
  
  # Create a binary vector indicating presence/absence of DEGs
  binary_vector <- as.integer(unique(unlist(degs_featureCounts)) %in% genes)
  
  return(binary_vector)
})

# Convert the list to a data frame
degs_binary_matrix <- as.data.frame(do.call(cbind, degs_binary_matrix))

# Check the dimensions of the rebuilt binary matrix
print(paste("Number of rows in binary matrix after rebuilding:", nrow(degs_binary_matrix)))

# Check if any rows are invalid or contain only zeros
invalid_rows <- apply(degs_binary_matrix, 1, function(x) all(x == 0))

# Remove invalid rows
if (sum(invalid_rows) > 0) {
  degs_binary_matrix <- degs_binary_matrix[!invalid_rows, ]
}

# Check the new dimensions after cleaning
print(paste("Number of rows in binary matrix after removing invalid rows:", nrow(degs_binary_matrix)))



# Remove NA values from unique DEGs before assigning row names
unique_degs <- unique(unlist(degs_featureCounts))
unique_degs <- unique_degs[!is.na(unique_degs)]

# Check if the number of unique DEGs matches the cleaned binary matrix
if (length(unique_degs) == nrow(degs_binary_matrix)) {
  rownames(degs_binary_matrix) <- unique_degs
} else {
  print("Mismatch still exists between unique DEGs and binary matrix rows.")
  print(paste("Number of unique DEGs:", length(unique_degs)))
  print(paste("Number of rows in binary matrix:", nrow(degs_binary_matrix)))
}


library(grid)

# Create the upset plot
upset(degs_binary_matrix, 
      sets = names(degs_featureCounts), 
      main.bar.color = "dodgerblue", 
      set_size.show = F,
      order.by = "freq", 
      keep.order = TRUE,
      sets.x.label = NULL,           # Remove set size labels on the left
      text.scale = c(3, 3, 3, 3, 2, 3)
)

# save svg
ggsave("upset_plot_common_DEGs_timepoints.svg", width = 12, height = 8)

# Function to combine DEGs across time points for each cell type
combine_degs_by_cell_type <- function(degs_featureCounts) {
  # Create an empty list to hold combined DEGs by cell type
  cell_type_degs <- list()
  
  # Loop through each sample and combine DEGs by cell type
  for (sample_name in names(degs_featureCounts)) {
    # Extract the cell type from the sample name (e.g., "CEK", "CELi", "HD11")
    cell_type <- strsplit(sample_name, "_")[[1]][1]
    
    # Combine DEGs for each cell type
    if (!cell_type %in% names(cell_type_degs)) {
      cell_type_degs[[cell_type]] <- degs_featureCounts[[sample_name]]  # Initialize with first time point
    } else {
      # Union DEGs across time points for the same cell type
      cell_type_degs[[cell_type]] <- union(cell_type_degs[[cell_type]], degs_featureCounts[[sample_name]])
    }
  }
  
  return(cell_type_degs)
}

# Combine DEGs by cell type (ignoring time points)
cell_type_degs <- combine_degs_by_cell_type(degs_featureCounts)



# Combine DEGs by cell type (ignoring time points)
cell_type_degs <- combine_degs_by_cell_type(degs_featureCounts)

# Get the unique DEGs across all cell types, removing NA values
unique_genes <- unique(unlist(cell_type_degs))
unique_genes <- unique_genes[!is.na(unique_genes)]  # Remove NA values

# Create a binary matrix where each cell type is a column and rows represent unique genes
binary_matrix_cell_types <- sapply(cell_type_degs, function(genes) {
  # Create a binary vector indicating presence of each gene in this cell type
  as.integer(unique_genes %in% genes)
})

# Convert to data frame
binary_matrix_cell_types <- as.data.frame(binary_matrix_cell_types)

# Assign the unique genes as row names to the binary matrix
rownames(binary_matrix_cell_types) <- unique_genes

# Check if the number of unique genes matches the number of rows in the binary matrix
print(paste("Number of unique genes:", length(unique_genes)))
print(paste("Number of rows in binary matrix:", nrow(binary_matrix_cell_types)))

# Create the UpSet plot for cell types without time points
upset(binary_matrix_cell_types, 
      sets = names(cell_type_degs), 
      main.bar.color = "dodgerblue", 
      order.by = "freq", 
      keep.order = TRUE,
      text.scale = c(3, 3, 3, 3, 3, 3),  # Adjust text sizes
)


# print set sizes for each cell type
cell_type_degs <- lapply(cell_type_degs, function(x) x[!is.na(x)])
set_sizes <- sapply(cell_type_degs, length)
print(set_sizes)

# common genes within cell type for timepoints after removing NAs
cell_type_degs <- lapply(cell_type_degs, function(x) x[!is.na(x)])
common_genes <- Reduce(intersect, cell_type_degs)
common_genes

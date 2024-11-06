# Load necessary libraries
library(data.table)
library(ggplot2)
library(dplyr)
library(svglite)


# Set the path to the featureCounts directory
edgeR_featureCounts_dir <- "edgeR_featureCounts"

# Function to list files and their row counts
list_files_row_counts <- function(dir_path) {
  file_info_list <- list()
  # List directories starting with 'Results_'
  subdirs <- list.dirs(dir_path, recursive = FALSE, full.names = TRUE)
  subdirs <- subdirs[grepl("/Results_", subdirs)]
  for (subdir in subdirs) {
    # List CSV files in the subdirectory
    files <- list.files(subdir, pattern = "\\.csv$", full.names = TRUE)
    for (file_path in files) {
      # Count rows excluding header
      row_count <- nrow(fread(file_path, select = 1)) 
      file_info <- list(file_path = file_path, row_count = row_count)
      file_info_list[[length(file_info_list) + 1]] <- file_info
    }
  }
  return(file_info_list)
}





# Get row counts for the featureCounts directory only
row_counts_featureCounts <- list_files_row_counts(edgeR_featureCounts_dir)

# Convert list to data frame
df_featureCounts <- do.call(rbind, lapply(row_counts_featureCounts, function(x) data.frame(Source = "FeatureCounts", Sample = basename(dirname(x$file_path)), RowCount = x$row_count)))

# Clean and extract sample names if necessary
df_featureCounts$Sample <- gsub("Results_", "", df_featureCounts$Sample)

# Let's re-factor the Sample variable to ensure the desired order
# Define the conditions and time points
conditions <- c("CEK", "CELi", "HD11")
time_points <- c("6hpi", "12hpi", "24hpi")

# Generate the ordered list of samples, with spaces instead of underscores
ordered_samples <- as.vector(sapply(conditions, function(c) paste(c, time_points, sep=" ")))
ordered_samples <- unlist(ordered_samples)

# Ensure that sample names match the expected format
df_featureCounts$Sample <- gsub("_", " ", df_featureCounts$Sample)

# Update the Sample variable to reflect this ordered factor
df_featureCounts$Sample <- factor(df_featureCounts$Sample, levels = ordered_samples)

# Example: Creating a new column 'Group' to distinguish between CEK, CELi, and HD11
df_featureCounts$Group <- ifelse(grepl("CEK", df_featureCounts$Sample), "CEK",
                                 ifelse(grepl("CELi", df_featureCounts$Sample), "CELi", "HD11"))

# Example: Creating a new column 'Time' to distinguish time points
df_featureCounts$Time <- ifelse(grepl("6hpi", df_featureCounts$Sample), "06hpi",
                                ifelse(grepl("12hpi", df_featureCounts$Sample), "12hpi", "24hpi"))

# print the counts of DEGs for each combination
summary_DEGs <- df_featureCounts %>% group_by(Group, Time) %>% summarise(RowCount = sum(RowCount))
summary_DEGs

# # save the table in a text file
# write.table(summary_DEGs, "featureCounts_DEGs_summary.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Plotting only FeatureCounts data
font_size = 24
title_size = 36
p <- ggplot(df_featureCounts, aes(x = Time, y = RowCount, fill = Time)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_linedraw() +
  facet_wrap(~Group, scales = "free_x", ncol = 3) +  # Adding faceting
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = font_size, family = "Times New Roman"),
        axis.text.y = element_text(size = 20, family = "Times New Roman"),
        panel.spacing = unit(2, "lines"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.title = element_text(size = font_size, family = "Times New Roman", face = "bold"),
        legend.text = element_text(size = font_size, family = "Times New Roman"),
        plot.title = element_text(size = font_size, face = "bold", family = "Times New Roman"),
        axis.title = element_text(size = font_size, family = "Times New Roman"),
        strip.text = element_text(size = 24, face = "bold", family = "Times New Roman")) +
  
  labs(y = "Number of DEGs Identified", 
       title = "Number of Differentially Expressed Genes Identified",
       fill = "Time (hpi)")

# Display the plot
print(p)

# save svg
ggsave("featureCounts_DEGs_by_cell_type.svg", plot = p, device = "svg", width = 12, height = 8)



p <- ggplot(df_featureCounts, aes(x = Group, y = RowCount, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_linedraw() +
  facet_wrap(~Time, scales = "free_x", ncol = 3) +  # Adding faceting
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = font_size, family = "Times New Roman"),
        axis.text.y = element_text(size = 20, family = "Times New Roman"),
        panel.spacing = unit(2, "lines"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.title = element_text(size = font_size, family = "Times New Roman", face = "bold"),
        legend.text = element_text(size = font_size, family = "Times New Roman"),
        plot.title = element_text(size = font_size, face = "bold", family = "Times New Roman"),
        axis.title = element_text(size = font_size, family = "Times New Roman"),
        strip.text = element_text(size = 24, face = "bold", family = "Times New Roman")) +
  
  labs(y = "Number of DEGs Identified", 
       title = "Number of Differentially Expressed Genes Identified",
       fill = "Cell Type")

# Display the plot
print(p)

# save svg
ggsave("featureCounts_DEGs_by_time.svg", plot = p, device = "svg", width = 12, height = 8) 

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

# # Save the DEGs data frame to a text file
# write.table(degs_featureCounts, "DEGs_by_group_time.txt", sep = "\t", quote = FALSE, row.names = FALSE)


# print number of non NA rows from each column
summary_DEGs <- sapply(degs_featureCounts, function(x) sum(!is.na(x)))
summary_DEGs


# Load necessary libraries
library(ggplot2)
library(dplyr)

# Load data
counts_matrix <- read.csv("final_counts_matrix.csv", row.names = 1)
metadata <- read.csv("metadata.csv")

# Transpose counts matrix for PCA (samples as rows, genes as columns)
counts_matrix_t <- t(counts_matrix)

# Remove columns with zero or near-zero variance
counts_matrix_t <- counts_matrix_t[, apply(counts_matrix_t, 2, var) > 0]

# Run PCA
pca_result <- prcomp(counts_matrix_t, scale. = TRUE)

# Create a data frame with PCA results and metadata
pca_df <- as.data.frame(pca_result$x)
pca_df$Sample_ID <- rownames(pca_df)
pca_df <- merge(pca_df, metadata, by.x = "Sample_ID", by.y = "Sample_ID")


ggplot(pca_df, aes(x = PC1, y = PC2, color = Treatment, shape = Cell_type, size = Time_hpi)) +
  geom_point() +
  theme_linedraw() +
  # facet_grid(~Cell_type, scale = "free_y") +
  # change font size of the grid text
  theme(strip.text = element_text(size = 28, family="Times New Roman")) +
  # increase panel spacing
  theme(panel.spacing = unit(2, "lines")) +
  labs(x = "PC1", y = "PC2") +
  theme(text = element_text(size = 20, family="Times New Roman")) +
  # adjust legend key size
  theme(legend.key.size = unit(1.5, "cm")) +  # This controls spacing
  # control the size of points in the legend
  guides(
    color = guide_legend(override.aes = list(size = 5)),  # Adjust color legend point size
    shape = guide_legend(override.aes = list(size = 5))  # Adjust shape legend point siz   # Adjust size legend point size
  ) +
  scale_color_discrete(name = "Treatment") +
  scale_size_discrete(name = "Time (hpi)") +
  scale_shape_discrete(name = "Cell Type")


# cell type faceted graphs
ggplot(pca_df, aes(x = PC1, y = PC2, color = Treatment, shape = Cell_type, size = Time_hpi)) +
  geom_point() +
  theme_linedraw() +
  facet_wrap(~Cell_type, scales = "free") +
  # change font size of the grid text
  theme(strip.text = element_text(size = 28, family="Times New Roman")) +
  # increase panel spacing
  theme(panel.spacing = unit(2, "lines")) +
  labs(x = "PC1", y = "PC2") +
  theme(text = element_text(size = 20, family="Times New Roman")) +
  # adjust legend key size
  theme(legend.key.size = unit(1.5, "cm")) +  # This controls spacing
  # control the size of points in the legend
  guides(
    color = guide_legend(override.aes = list(size = 5)),  # Adjust color legend point size
    shape = guide_legend(override.aes = list(size = 5))  # Adjust shape legend point siz   # Adjust size legend point size
  ) +
  scale_color_discrete(name = "Treatment") +
  scale_size_discrete(name = "Time (hpi)") +
  scale_shape_discrete(name = "Cell Type")

# permanova
# Calculate Euclidean distance matrix
distance_matrix <- dist(counts_matrix_t, method = "euclidean")

# PERMANOVA for Cell_type
permanova_cell_type <- adonis2(distance_matrix ~ Cell_type, data = metadata, permutations = 999)
print(permanova_cell_type)

# PERMANOVA for Time_hpi
permanova_time_hpi <- adonis2(distance_matrix ~ Time_hpi, data = metadata, permutations = 999)
print(permanova_time_hpi)

# PERMANOVA for Treatment
permanova_treatment <- adonis2(distance_matrix ~ Treatment, data = metadata, permutations = 999)
print(permanova_treatment)


# PERMANOVA for Treatment
permanova_treatment <- adonis2(distance_matrix ~ Treatment + Time_hpi + Cell_type, data = metadata, permutations = 999)
print(permanova_treatment)


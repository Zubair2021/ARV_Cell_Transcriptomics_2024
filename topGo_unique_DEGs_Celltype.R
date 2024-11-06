# Load necessary libraries
library(readxl)

# Read the Excel file
deg_data <- read_excel("combined_list_DEGs.xlsx")



# Extract the gene lists for each cell type
CEK_genes <- deg_data$CEK
# remove nas 
CEK_genes <- CEK_genes[!is.na(CEK_genes)]
length(CEK_genes)

CELi_genes <- deg_data$CELi
CELi_genes <- CELi_genes[!is.na(CELi_genes)]
length(CELi_genes)

HD11_genes <- deg_data$HD11
HD11_genes <- HD11_genes[!is.na(HD11_genes)]
length(HD11_genes)


# Load necessary libraries
library(readr)
library(topGO)

# Load necessary library
library(readr)

# Load the CSV file and set the first column as rownames
counts_data <- read_csv("final_counts_matrix.csv")
colnames(counts_data)

# Now, the rownames will be the gene IDs
gene_universe <- counts_data$Gene_ID


# Create named vectors indicating DEGs for each cell type (DEGs as 1, others as 0)
CEK_list <- factor(as.integer(gene_universe %in% CEK_genes))
names(CEK_list) <- gene_universe

CELi_list <- factor(as.integer(gene_universe %in% CELi_genes))
names(CELi_list) <- gene_universe

HD11_list <- factor(as.integer(gene_universe %in% HD11_genes))
names(HD11_list) <- gene_universe

# Load the appropriate annotation package for your organism
# For human, use 'org.Hs.eg.db'. For mouse, use 'org.Mm.eg.db', etc.
library(org.Gg.eg.db)  # Replace with appropriate organism if not human

# Update the run_topGO function to use annFUN.org
run_topGO <- function(gene_list, gene_selection_function) {
  GO_data <- new("topGOdata",
                 ontology = "BP",  # You can change to "MF" or "CC"
                 allGenes = gene_list,
                 geneSel = gene_selection_function,
                 annot = annFUN.org,  # Use annFUN.org instead of annFUN.db
                 mapping = "org.Gg.eg.db",  # Replace with your organism's database
                 ID = "symbol")  # Assuming gene IDs are gene symbols
  
  # Run the enrichment test
  resultFisher <- runTest(GO_data, algorithm = "classic", statistic = "fisher")
  
  # Get results
  allGO <- GenTable(GO_data, classicFisher = resultFisher, orderBy = "classicFisher", topNodes = 10)
  return(allGO)
}

# Function to select genes (significant genes)
gene_selection <- function(x) x == 1


# Run for CEK
CEK_GO_results <- run_topGO(CEK_list, gene_selection)
write.csv(CEK_GO_results, "CEK_GO_results.csv")

# Run for CELi
CELi_GO_results <- run_topGO(CELi_list, gene_selection)
CELi_GO_results
write.csv(CELi_GO_results, "CELi_GO_results.csv")

# Run for HD11
HD11_GO_results <- run_topGO(HD11_list, gene_selection)
HD11_GO_results
write.csv(HD11_GO_results, "HD11_GO_results.csv")



# Update the run_topGO function to use annFUN.org
run_topGO <- function(gene_list, gene_selection_function) {
  GO_data <- new("topGOdata",
                 ontology = "MF",  # You can change to "MF" or "CC"
                 allGenes = gene_list,
                 geneSel = gene_selection_function,
                 annot = annFUN.org,  # Use annFUN.org instead of annFUN.db
                 mapping = "org.Gg.eg.db",  # Replace with your organism's database
                 ID = "symbol")  # Assuming gene IDs are gene symbols
  
  # Run the enrichment test
  resultFisher <- runTest(GO_data, algorithm = "classic", statistic = "fisher")
  
  # Get results
  allGO <- GenTable(GO_data, classicFisher = resultFisher, orderBy = "classicFisher", topNodes = 10)
  return(allGO)
}

# Function to select genes (significant genes)
gene_selection <- function(x) x == 1

# Run for CEK
CEK_GO_MF_results <- run_topGO(CEK_list, gene_selection)
write.csv(CEK_GO_MF_results, "CEK_GO_MF_results.csv")

# Run for CELi
CELi_GO_MF_results <- run_topGO(CELi_list, gene_selection)
CELi_GO_MF_results
write.csv(CELi_GO_MF_results, "CELi_GO_MF_results.csv")

# Run for HD11
HD11_GO_MF_results <- run_topGO(HD11_list, gene_selection)
HD11_GO_MF_results
write.csv(HD11_GO_MF_results, "HD11_GO_MF_results.csv")



# Load necessary libraries
library(ggplot2)
library(dplyr)

# Function to calculate the expected number of significant genes for each GO term
calculate_observed_expected_ratio <- function(GO_results, gene_universe_size, significant_genes_count) {
  GO_results$Annotated <- as.numeric(GO_results$Annotated)  # Ensure 'Annotated' is numeric
  GO_results$Significant <- as.numeric(GO_results$Significant)  # Ensure 'Significant' is numeric
  
  # # Calculate the expected number of significant genes based on the proportion of annotated genes
  # GO_results$Expected <- (GO_results$Annotated / gene_universe_size) * significant_genes_count
  
  # Calculate the observed vs expected ratio
  GO_results$GeneRatio <- GO_results$Significant / GO_results$Annotated
  
  return(GO_results)
}

# Load your GO results
CEK_GO_BP_results <- read.csv("CEK_GO_results.csv")
CELi_GO_BP_results <- read.csv("CELi_GO_results.csv")
HD11_GO_BP_results <- read.csv("HD11_GO_results.csv")

# Define total size of the gene universe and the number of significant genes for each cell type
gene_universe_size <- length(gene_universe)
CEK_significant_genes_count <- length(CEK_genes)
CELi_significant_genes_count <- length(CELi_genes)
HD11_significant_genes_count <- length(HD11_genes)

# Calculate observed vs expected ratio for each cell type
CEK_GO_BP_results <- calculate_observed_expected_ratio(CEK_GO_BP_results, gene_universe_size, CEK_significant_genes_count)
CELi_GO_BP_results <- calculate_observed_expected_ratio(CELi_GO_BP_results, gene_universe_size, CELi_significant_genes_count)
HD11_GO_BP_results <- calculate_observed_expected_ratio(HD11_GO_BP_results, gene_universe_size, HD11_significant_genes_count)

# Add cell type labels to each dataset
CEK_GO_BP_results$cell_type <- "CEK"
CELi_GO_BP_results$cell_type <- "CELi"
HD11_GO_BP_results$cell_type <- "HD-11"

# Combine all the data into one data frame
combined_GO_results <- bind_rows(CEK_GO_BP_results, CELi_GO_BP_results, HD11_GO_BP_results)


# modify the term column to sentence case 
combined_GO_results$Term <- str_to_sentence(combined_GO_results$Term)

# use short names for the terms
combined_GO_results$Term <- str_replace_all(combined_GO_results$Term, " regulation ", " reg. ")
combined_GO_results$Term <- str_replace_all(combined_GO_results$Term, " response ", " resp. ")
combined_GO_results$Term <- str_replace_all(combined_GO_results$Term, "Response ", "Resp. ")
combined_GO_results$Term <- str_replace_all(combined_GO_results$Term, " signaling ", " signal. ")
combined_GO_results$Term <- str_replace_all(combined_GO_results$Term, "Negative", "Neg.")
combined_GO_results$Term <- str_replace_all(combined_GO_results$Term, "Biological process involved in interspec...", "Interspecies Interaction")


# Create the dot plot with GeneRatio as size and color
ggplot(combined_GO_results, aes(str_wrap(Term, width = 25), y = GeneRatio, color = classicFisher , size = GeneRatio)) +
  geom_point() +
  facet_wrap(~ cell_type, scales = "free_y") +
  coord_flip() +
  scale_color_gradient(low = "red", high = "blue") +
  labs(x = "GO Term", y = "Gene Ratio (Significant / Annotated)", color = "P-value (Fisher)", size = "Gene Ratio") +
  theme_linedraw() +
  theme(text = element_text(size = 12, family = "Times New Roman"),
        axis.text.y = element_text(size = 14), 
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        # x axis title
        axis.title.x = element_text(size = 16, vjust =0.3),
        axis.title.y = element_text(size = 16),
        strip.text = element_text(size = 18))

# save svg
ggsave("GO_BP_dotplot.svg", width = 12, height = 8, units = "in")


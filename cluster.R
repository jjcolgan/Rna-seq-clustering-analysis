# Load necessary libraries
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(DESeq2)
library(tidyverse)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(DOSE)
library(enrichplot)
library(ggupset)
library(ashr)
library(cowplot)
library(enrichR)

# Load gene count data and set 'Geneid' as row names
fulldata <- read.csv('gene_matrix_count.csv')
fulldata <- column_to_rownames(fulldata, 'Geneid')

# Remove sample identifiers and extra characters from column names
colnames(fulldata) <- sub("_S\\d+.*$", "", colnames(fulldata))

# Make a copy of fulldata for processing
counts <- fulldata
counts <- counts[, order(names(counts))]  # Sort columns alphabetically

# Load metadata and set 'library' as row names
metadata <- read.csv('metadata.csv')
metadata <- column_to_rownames(metadata, 'library')
metadata <- metadata[order(rownames(metadata)), ]  # Sort metadata rows

# Transpose fulldata for merging
fulldata <- t(fulldata)

# Remove specific samples from metadata
metadata <- metadata[row.names(metadata) != 'EC_CG_30', ]
metadata <- metadata[row.names(metadata) != 'EC_CG_44', ]
metadata <- metadata[row.names(metadata) != 'EC_CG_06', ]

# Create a new 'sampleDay' column in metadata by combining row names and 'Day' column
metadata$sampleDay <- paste(rownames(metadata), metadata$Day, sep = "_")

# Merge transposed counts and metadata by row names
full <- merge(metadata, t(counts), by = 'row.names')

# Subset the data for specific cohorts and tissue type "Colon"
justBarseqSI <- full[(full$Cohort == "Barseq" | full$Cohort == "Conventionalized" | full$Cohort == "2 week GF") & full$Tissue == "Colon", ]

# Split data into counts and metadata for selected samples
SICounts <- justBarseqSI[c(19:56960)]  # Gene expression counts
SIMetadata <- justBarseqSI[c(1:19)]  # Sample metadata

# Clean up row names for SICounts and SIMetadata
SICounts <- as.data.frame(SICounts)
rownames(SICounts) <- NULL
SICounts <- column_to_rownames(SICounts, 'sampleDay')

SIMetadata <- as.data.frame(SIMetadata)
rownames(SIMetadata) <- NULL
SIMetadata <- column_to_rownames(SIMetadata, 'sampleDay')

# Transpose SICounts for further filtering
SICounts <- t(SICounts)

# Record the number of features (rows) before filtering
preFilterFeatures <- nrow(SICounts)

# Set Group in SIMetadata as a factor and display the levels
SIMetadata$Group <- as.factor(SIMetadata$Group)
levels(SIMetadata$Group)

# Abundance filter to keep genes expressed in at least 'smallestGroupSize' samples with counts >= 15
smallestGroupSize <- 4
keep <- rowSums(SICounts >= 15) >= smallestGroupSize
SICounts <- SICounts[keep, ]

# Record the number of features after filtering for abundance
abundanceFilterfeatures <- nrow(SICounts)

# Load additional coding gene data and filter for 'protein_coding' transcripts
coding <- read_tsv("mart_export (2).txt")
coding <- coding[coding$`Transcript type` == 'protein_coding', ]
coding <- column_to_rownames(coding, 'Gene stable ID')

# Keep only genes present in both coding and SICounts data
keep <- intersect(rownames(coding), rownames(SICounts))
length(keep)  # Number of genes kept
SICounts <- SICounts[rownames(SICounts) %in% keep, ]

# Display summary statistics for row means of SICounts
summary(rowMeans(SICounts))

# Create DESeq2 dataset using SICounts and SIMetadata, with design ~Group
dds <- DESeqDataSetFromMatrix(countData = as.matrix(SICounts), colData = SIMetadata, design = ~Group)

#Run Deseq2
dds <- DESeq(dds)
dim(dds)

# Extract normalized counts as a data frame from the DESeq2 dataset
resdata <- as.data.frame(counts(dds, normalized = TRUE))

# Initialize an empty list to store selected genes across comparisons
total_selec <- list()
x <- 1  # Counter to track levels in Group factor

# Loop through each level of Group to compare against other levels
for(i in levels(SIMetadata$Group)) {
  x <- x + 1
  
  # Continue only if the counter is within the range of Group levels
  if (x <= length(levels(SIMetadata$Group))) {
    
    # Loop to perform pairwise comparisons of Group levels
    for(j in x:length(levels(SIMetadata$Group))) {
      
      # Run DESeq2's results function for each pairwise contrast in Group levels
      res <- results(dds, contrast = c("Group", i, levels(SIMetadata$Group)[j]))
      
      # Create a unique identifier for each comparison pair
      d <- paste(i, levels(SIMetadata$Group)[j], sep = "&")
      
      # Add gene names to results and convert to a data frame
      res$genenames <- rownames(res)
      resul <- as.data.frame(res)
      
      # Filter for significant differentially expressed genes (log2 fold change > 0.5 or < -0.5)
      # and adjusted p-value < 0.05
      significantDE <- resul[resul$log2FoldChange > 0.5 | resul$log2FoldChange < -0.5, ]
      significantDE <- significantDE[significantDE$padj < 0.05, ]
      significantDE <- na.omit(significantDE)  # Remove any NA values
      
      # Extract gene names of significant DE genes and add to total_selec list
      selec <- as.list(significantDE$genenames)
      total_selec <- append(total_selec, selec)
    }
  }
}

# Remove duplicates from total_selec list and convert to a data frame
total_selec <- c(unique(total_selec))
total_selec <- t(as.data.frame(total_selec))

# Subset resdata to include only selected genes based on total_selec list
selection <- resdata[total_selec[, 1], ]

# Display dimensions of the selection data frame (number of genes selected)
dim(selection)

# Convert the selected data to a matrix and scale the data (row-wise scaling)
DEgenes <- as.matrix(selection)
scaledata <- t(scale(t(DEgenes)))

# Cluster by samples
hc <- hclust(as.dist(1-cor(scaledata, method="spearman")), method="complete")

#View sample clusters
sampleTree = as.dendrogram(hc, method="average")
plot(sampleTree,
     main = "Sample Clustering",
     ylab = "Height")

# Cluster genes
hr <- hclust(as.dist(1-cor(t(scaledata), method="spearman")), method="complete") 

# Calculate the initial total within-cluster sum of squares (WSS) for the entire dataset as a single cluster
wss <- (nrow(scaledata) - 1) * sum(apply(scaledata, 2, var))

# Loop through cluster numbers from 2 to 20 to calculate WSS for each k-means solution
for (i in 2:20) {
  # Run k-means clustering with 'i' clusters and store the WSS for the solution
  wss[i] <- sum(kmeans(scaledata, centers = i)$withinss)
}

# Plot the WSS values against the number of clusters to identify the optimal cluster count using the "elbow method"
plot(1:20, wss, type = "b", xlab = "Number of Clusters", ylab = "Within groups sum of squares")

# Define a function to calculate the mean of each cluster's expression levels
clust.core <- function(i, dat, clusters) {
  ind <- (clusters == i)  # Logical index for samples in the current cluster
  colMeans(dat[ind, ])    # Calculate column means for selected samples
}

# Assign clusters and calculate core expression values for each unique cluster
clusters <- hclustk6
cores <- sapply(unique(clusters), clust.core, scaledata, clusters)

# Display column names of scaledata to verify structure
colnames(scaledata)

# Convert core expression data to a long format for easier plotting
moltenCores <- melt(cores)
colnames(moltenCores) <- c('sampleDay', 'cluster', 'value')

# Merge the molten cores data with metadata on 'sampleDay'
moltenCores <- merge(moltenCores, SIMetadata, by = 'sampleDay')

# Define the desired order of Group levels and apply it as a factor
desired_order <- c("14dayGF", "1dayBar-seq", "3dayBar-seq", "7dayBar-seq", "14dayBar-seq", "Conventionalized")
moltenCores$Group <- factor(moltenCores$Group, levels = desired_order)

# Subset moltenCores by each cluster (1 through 4)
cluster1 <- moltenCores[moltenCores$cluster == 1, ]
cluster2 <- moltenCores[moltenCores$cluster == 2, ]
cluster3 <- moltenCores[moltenCores$cluster == 3, ]
cluster4 <- moltenCores[moltenCores$cluster == 4, ]

# Extract genes and their expressions for each cluster and store dimensions
clust1 <- t(scaledata[hclustk6 == 1, ])
nGenesClust1 <- dim(clust1)
clust2 <- t(scaledata[hclustk6 == 2, ])
nGenesClust2 <- dim(clust2)
clust3 <- t(scaledata[hclustk6 == 3, ])
nGenesClust3 <- dim(clust3)
clust4 <- t(scaledata[hclustk6 == 4, ])
nGenesClust4 <- dim(clust4)

# Transform cluster data to long format and add a 'cluster' column
clust1Longer <- clust1 %>%
  as.data.frame() %>%
  rownames_to_column('sampleDay') %>%
  pivot_longer(cols = 2:(nGenesClust1[2] + 1), names_to = 'gene') %>%
  dplyr::rename('scaledExpression' = value) %>%
  mutate('cluster' = 1)

clust2Longer <- clust2 %>%
  as.data.frame() %>%
  rownames_to_column('sampleDay') %>%
  pivot_longer(cols = 2:(nGenesClust2[2] + 1), names_to = 'gene') %>%
  dplyr::rename('scaledExpression' = value) %>%
  mutate('cluster' = 2)

clust3Longer <- clust3 %>%
  as.data.frame() %>%
  rownames_to_column('sampleDay') %>%
  pivot_longer(cols = 2:(nGenesClust3[2] + 1), names_to = 'gene') %>%
  dplyr::rename('scaledExpression' = value) %>%
  mutate('cluster' = 3)

# Rename Group levels for readability
levels(moltenCores$Group)[levels(moltenCores$Group) == '1dayBar-seq'] <- 'Day 1'
levels(moltenCores$Group)[levels(moltenCores$Group) == '3dayBar-seq'] <- 'Day 3'
levels(moltenCores$Group)[levels(moltenCores$Group) == '7dayBar-seq'] <- 'Day 7'
levels(moltenCores$Group)[levels(moltenCores$Group) == '14dayBar-seq'] <- 'Day 14'
levels(moltenCores$Group)[levels(moltenCores$Group) == '14dayGF'] <- 'GF'

# Combine clusters, merge with metadata, and calculate mean expression levels
bp1 <- rbind(clust1Longer, clust2Longer) %>%
  merge(moltenCores, by = c('sampleDay', 'cluster')) %>%
  group_by(gene, Group, cluster) %>%
  summarise(
    'meanGeneExpression' = mean(scaledExpression),
    'meanClusterExpression' = mean(value)
  ) %>%
  
  # Plot with ggplot2
  ggplot() +
  geom_line(aes(x = Group, y = meanGeneExpression, group = gene), col = 'lightblue') +
  
  # Facet wrap with customized labels
  facet_wrap(~cluster, nrow = 1, labeller = as_labeller(
    c(
      `1` = paste0('Cluster 1 (', nGenesClust1[2], ' genes)'),
      `2` = paste0('Cluster 2 (', nGenesClust2[2], ' genes)'),
      `3` = paste0('Cluster 3 (', nGenesClust3[2], ' genes)')
    )
  )) +
  
  # Add lines and points for mean cluster expression
  geom_line(aes(x = Group, y = meanClusterExpression, group = 1), col = 'blue') +
  geom_point(aes(x = Group, y = meanClusterExpression, col = Group)) +
  
  # Apply theme and labels
  theme_bw() +
  theme(
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 15)
  ) +
  
  # Adjust x-axis breaks and labels
  scale_x_discrete(breaks = NULL) +
  labs(y = 'Scaled expression', x = 'Time')

# Display the plot
bp1

# Example of cluster identification
# Perform Gene Ontology (GO) enrichment analysis on genes in cluster 2
clust2Go <- enrichGO(
  gene = colnames(clust2),                # Specify genes for enrichment analysis
  OrgDb = org.Mm.eg.db,                   # Use the mouse genome database
  keyType = 'ENSEMBL',                    # Specify that genes are in ENSEMBL format
  readable = TRUE,                        # Convert ENSEMBL IDs to gene symbols
  ont = 'bp'                              # Focus on Biological Processes
)

# Simplify the GO results to reduce redundancy
clust2go <- clusterProfiler::simplify(clust2Go)

# Generate a bar plot for the top 10 enriched biological processes in cluster 2
p6 <- barplot(
  clust2Go,                               # Data for plotting
  showCategory = 10,                      # Show top 10 categories
  font.size = 8,                          # Set font size
  title = 'Enriched biological\nprocesses cluster 1'  # Title for the plot
)
p6  # Display the plot

# Convert clust2Go results to a data frame
clust2Go %>%
  as.data.frame()

# Map ENSEMBL gene IDs in cluster 2 to gene symbols
queryGenes <- bitr(
  geneID = colnames(clust2),              # Genes to convert
  fromType = 'ENSEMBL',                   # Original gene ID format
  toType = 'SYMBOL',                      # Target gene ID format
  OrgDb = org.Mm.eg.db                    # Use mouse genome database
) %>%
  .$SYMBOL                                # Extract SYMBOL column for use in enrichment

# Perform additional pathway enrichment using Enrichr for multiple databases
out <- enrichr(
  databases = c('GO_Biological_Process_2023', 'KEGG_2019_Mouse', 'WikiPathways_2024_Mouse'),
  genes = queryGenes                       # Genes for enrichment analysis
)

# Filter results for GO Biological Processes with an adjusted p-value < 0.05
out[1] %>%
  as.data.frame() %>%
  filter(`GO_Biological_Process_2023.Adjusted.P.value` < 0.05)

# Filter results for KEGG pathways with an adjusted p-value < 0.05
out[2] %>%
  as.data.frame() %>%
  filter(`KEGG_2019_Mouse.Adjusted.P.value` < 0.05)

# Filter results for WikiPathways with an adjusted p-value < 0.05
out[3] %>%
  as.data.frame() %>%
  filter(`WikiPathways_2024_Mouse.Adjusted.P.value` < 0.05)


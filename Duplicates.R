# This script performs duplicate analysis, using PLINK identity-by-descent (PI_HAT) metrics.
# Genotyping results are exported from Axiom Analysis Suite in PLINK .ped/.map file format.

# Only a list of 5194 quality SNPs common across the 20K/50K/480K arrays are exported (based off a list provided by Nicholas Howard).
# Aaron's MSc samples are combined with genotype data from other available data sets (New Zealand and international datasets).

# PLINK version used: PLINK v1.9.0-b.7.7 64-bit (22 Oct 2024)

# Load packages -----------------------------------------------------------
library(tibble)
library(igraph)
library(dplyr)
library(tidyr)


# Set working directory ---------------------------------------------------

setwd("C:/Users/curly/Desktop/Apple Genotyping/Methods/LowQC_ID/Inputs/Duplicates")



# Run PLINK duplicate analysis --------------------------------------------

# Create a filtered subset of SNPs using PLINK
system("plink --file All_Geno  --missing-genotype 0 --geno 0.02 --maf 0.05 --make-bed --out Filtered_Geno")

# Convert PLINK binary files (.bed/.bim/.fam) to PLINK text files (.ped/.map)
system("plink --bfile Filtered_Geno --recode --out Check_Dupe")

# Run duplicate analysis using the filtered SNP subset
system("plink --file Check_Dupe --missing-genotype 0 --genome full")


# Save PLINK .genome file as .txt -----------------------------------------

# Read .genome file
genome <- read.table("plink.genome", header = TRUE, sep = "", stringsAsFactors = FALSE)
write.table(genome, "C:/Users/curly/Desktop/Apple Genotyping/Results/LowQC_ID/Duplicates/PLINK_genome.txt", sep = "\t", row.names = FALSE, quote = FALSE)


# Group duplicate samples together ----------------------------------------

# Filter for PI_HAT > 0.90 (duplicate threshold; manually assess any between 0.90 and 0.97)
dupes <- genome[!(genome$PI_HAT < 0.90),]

# Save a .txt file of duplicate pairs
write.table(dupes,"C:/Users/curly/Desktop/Apple Genotyping/Results/LowQC_ID/Duplicates/Duplicate_pairs.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Keep only IDs
dupes_ID <- subset(dupes, select = c("IID1","IID2"))

# Group duplicate IDs with igraph
graph <- graph_from_data_frame(dupes_ID, directed = FALSE)
components <- components(graph)

# Sort groupings by number of duplicates
group_sizes <- table(components$membership)
sorted_group_ids <- order(group_sizes)
new_ids <- match(components$membership, sorted_group_ids)
V(graph)$group <- new_ids
grouped_samples <- split(names(components$membership), new_ids)

# Pad group with length less than max length with NA values
max_len <- max(sapply(grouped_samples, length))
padded_list <- lapply(grouped_samples, function(x) {c(x, rep(" ", max_len - length(x)))})

# Write groupings to a dataframe, add number of samples per group, format
dupe_groups <- as.data.frame(do.call(rbind, padded_list))
dupe_groups <- cbind(Group = seq_len(nrow(dupe_groups)), dupe_groups)
sample_counts <- rowSums(dupe_groups[, -1] != " ")
dupe_groups <- add_column(dupe_groups, SampleCount = sample_counts, .after = "Group")
colnames(dupe_groups) <- c("Group", "SampleCount", "ID1","ID2","ID3","ID4","ID5","ID6")

# Save .csv of duplicate ID groupings
write.csv(dupe_groups,"C:/Users/curly/Desktop/Apple Genotyping/Results/LowQC_ID/Duplicates/Duplicate_groups.csv", row.names = FALSE)

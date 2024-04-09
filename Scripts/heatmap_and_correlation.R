#install.packages(c("tidyverse", "viridis", "pheatmap"))
#install.packages("Hmisc")
library(tidyverse)
library(viridis)
library(pheatmap)
library(ggplot2)
library(Hmisc)
library(readr)
library(corrplot)

# Replace 'data.txt' with the actual path to your file
data <- read.table("cocultivations_v2_corrected.txt", header = TRUE, sep = "\t", row.names = 1)
# remove taxa appearing only once
filtered_data <- data[rowSums(data > 0) >= 2, ]
# Create a more visually appealing heatmap using pheatmap
library(pheatmap)

# Use the 'RdBu' color scale for better contrast
pheatmap(filtered_data, 
         color = colorRampPalette(c("white", "orange", "red"))(200), 
         cluster_rows = TRUE, 
         cluster_cols = TRUE)

####spearman correlation

# Read the CSV file
data <- read.table("cocultivations_v2_corrected.txt", header = TRUE, sep = "\t", row.names = 1)

# Transpose the data frame
transposed_data <- t(data)

#correlation of all taxa 
res <- cor(transposed_data, method = "spearman")
res <- round(res, 10)
res
write.table(res, file = "spearman_correlation.txt", row.names = TRUE, col.names = TRUE,sep = "\t")

#according p-values
res2 <- rcorr(res, type = c("spearman"))
p_value <- round (res2[["P"]], 10)
p_value
write.table(p_value, file = "spearman_correlation_pvalue.txt", row.names = TRUE, col.names = TRUE,sep = "\t")

#heatmap for correlation matrix
corrplot (res, tl.cex = 0.5)




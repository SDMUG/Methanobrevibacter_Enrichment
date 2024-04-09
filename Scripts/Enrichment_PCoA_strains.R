
devtools::install_github(c("david-barnett/microViz","jrnold/ggthemes","jbisanz/qiime2R","gmteunisse/Fantaxtic"))
BiocManager::install("microbiomeMarker")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("lefser")


library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("plyr"); packageVersion("plyr")
library("cluster"); packageVersion("cluster")
library("igraph"); packageVersion("igraph")
library("markovchain"); packageVersion("markovchain")
library("RColorBrewer")
library("phytools")
library("gridExtra")
library("grid")
library("ggthemes")
library("ggpubr")
library("dplyr")
library("cowplot")
library("ape")
library("devtools")
library("vegan")
library(cluster)
library(tidyverse)
library(microbiomeMarker)
library("lefser")





#################################################################################

##### Import .txt files #####

#################################################################################


#load data and prepare it for phyloseq

#set working directory
setwd("~/Enrichment")

#read in otu table
otu_table = read.table("eggnog_featuretable.txt", sep = "\t", row.names = 1, header = TRUE, check.names = FALSE)
otu_table = as.matrix(otu_table)
print(otu_table)

#read in taxonomy
#separated by kingdom phylum class order family genus species
taxonomy = read.table("tax_cog.txt", row.names = 1, header=TRUE, fill=TRUE)
taxonomy = as.matrix(taxonomy)


#read in metadata
metadata = read.table("Metadata.txt",row.names = 1, header = TRUE, sep = "\t")
# Create a data frame with the metadata
sample_data <- data.frame(metadata)

#unload microbiome maker again as it interferes with other commands later on
detach("package:microbiomeMarker", unload = TRUE)



#import as phyloseq objects
OTU = otu_table(otu_table, taxa_are_rows = TRUE)
TAX = tax_table(taxonomy)
META = sample_data(metadata)

#check that your OTU names are consistent across objects
taxa_names(TAX)
taxa_names(OTU)


# make sure files have the same sample names
sample_names(OTU)
sample_names(META)

#merge into one phyloseq object
ps = phyloseq(OTU, TAX, META)
ps


#################################################################################

##### make PCoA #####
## unconstrained ordinations ##

#################################################################################


# Ordinate using Principal Coordinate analysis, distance can be changed to "bray", "unifrac", "wunifrac"
ps_pcoa <- ordinate(
  physeq = ps, 
  method = "PCoA", 
  distance = "bray"
)

# Or write in short format
ps_pcoa <- ordinate(ps, "PCoA", "bray")

ps_pcoa

# checking the ordination by making a screen plot
plot_ordination(
  physeq = ps,
  ordination = ps_pcoa,
  type="scree")

# Plot 
PCoA_enrichment <- plot_ordination(
  physeq = ps,
  ordination = ps_pcoa,
  axes = c(1, 2),
  color = "Strain",
  shape = "Enrichment",
) +
  scale_color_manual(values = c("orange", "darkgreen")) +
  geom_point(aes(color = Strain), size = 3, show.legend = FALSE) +
  geom_point(colour = "grey90", size = 1.5) +
  geom_point(aes(shape = Enrichment), fill = "black") +
  stat_ellipse(aes(color = Strain, group = Strain)) +
  labs(title = "PCoA (bray)", subtitle = "COG") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )


PCoA_enrichment + scale_shape_manual(values = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16))

PCoA_enrichment

# Calculate bray curtis distance matrix
ps_bray <- phyloseq::distance(ps, method = "bray")

# make a data frame from the scaled sample_data
sampledf <- data.frame(sample_data(ps))


# Perform PERMANOVA test
permanova_result <- adonis2(ps_bray ~ Strain, data = sampledf)
permanova_result


# test of Homegeneity of dispersion
beta <- betadisper(ps_bray, sampledf$Strain)

# run a permutation test to get a statistic and a significance score
permutest(beta)

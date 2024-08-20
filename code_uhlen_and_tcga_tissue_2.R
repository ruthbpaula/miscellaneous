setwd("C:/Users/ruthb/Documents/UTH-MDA/Master's files/2020-2_and_2021-1/uhlen_tissue_expression/any_gene_locus")
setwd("C:/Users/ruthb/Documents/UTH-MDA/Master's files/2020-2_and_2021-1/new tcga analysis")
options(digits=22)
rm(list=ls())

path_uhlen <- "C:/Users/ruthb/Documents/UTH-MDA/Master's files/2020-2_and_2021-1/uhlen_tissue_expression/any_gene_locus/"
path_tcga <- "C:/Users/ruthb/Documents/UTH-MDA/Master's files/2020-2_and_2021-1/new tcga analysis/"

# Edit files in excel - average of tissues ----

# Reading CSV file

plate_raw <- read.csv('all_tumor_g4_avg_ok - Copy.txt', header = TRUE, sep = "\t")
plate <- as.data.frame(plate_raw)

# log2(fpkm + 1)
plate_subset <- subset(plate, select=kidney:rectum)  # for uhlen
plate_subset <- subset(plate, select=BLCA:UCEC)  # for TCGA

log2 <- log2(plate_subset + 1)

write.table(log2, file = "results_tumor_g4 - Copy.txt", sep = "\t", row.names = FALSE)

# heatmap
library(pheatmap)
library(circlize)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

data <- data.matrix(log2)

# colors <- c(min(data),seq(-2,2,by=0.01),max(data))
# my_palette <- c("green",colorRampPalette(colors = c("red", "yellow", "green"))
#                (n = length(colors)-3), "green")
# pheatmap(data, show_rownames = FALSE, color = my_palette)

#heatmap(data)

#heatmap.2(data)

my_palette_2 <- c("darkslategray2", "darkslategray3", "darkslategray4", "darkslategray", "black")
pheatmap(data, show_rownames = FALSE, color = my_palette_2)

#____

# Define breaks and corresponding colors
breaks0 <- c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20)
breaks1 <- c(0, 4, 8, 12, 16, 20)  # for TCGA
breaks2 <- c(0, 10, 20)  # to create continuous scale - for TCGA
breaks3 <- c(0, 5, 10, 15, 20, 25)
breaks4 <- c(0, 3, 6, 9, 12, 15)  # for Uhlen breaks
breaks5 <- c(0, 7.5, 15)  # to create continuous scale - for Uhlen

# Discrete scale
# colors <- c("darkslategray2", "darkslategray3", "darkslategray4", "darkslategray", "black")
colors <- c("pink","hotpink", "deeppink2", "deeppink4", "black")
# colors <- brewer.pal(n = 6, name = "PuRd")

# Continuous scale
colors = colorRamp2(breaks5, c("pink","deeppink2", "deeppink4"))
# Generate the colors from the color scale
n_colors <- 100  # Choose an appropriate number of colors
colors_vector <- colors(seq(0, 15, length.out = n_colors))

# Plotting - Discrete
pheatmap(data, col = colors, show_rownames = FALSE,
         breaks = breaks1,
         legend_breaks = breaks1)
# 455 width, 420 height

# Plotting - Continuous
pheatmap(data, col = colors_vector, show_rownames = FALSE,
         breaks = seq(0, 15, length.out = 100),
         legend_breaks = breaks4)
# 455 width, 420 height

#--

uhlen_no_g4 <- read.table(paste0(path_uhlen, "results_no_g4.txt"), header=T)
uhlen_g4 <- read.table(paste0(path_uhlen, "results_g4.txt"), header=T)
tcga_normal_no_g4 <- read.table(paste0(path_tcga, "results_normal_nog4.txt"), header=T)
tcga_normal_g4 <- read.table(paste0(path_tcga, "results_normal_g4.txt"), header=T)
tcga_tumor_no_g4 <- read.table(paste0(path_tcga, "results_tumor_nog4.txt"), header=T)
tcga_tumor_g4 <- read.table(paste0(path_tcga, "results_tumor_g4.txt"), header=T)

# setwd("C:/Users/ruthb/Documents/UTH-MDA/Master's files/2020-2_and_2021-1/uhlen_tissue_expression/any_gene_locus")  # uhlen
setwd("C:/Users/ruthb/Documents/UTH-MDA/Master's files/2020-2_and_2021-1/new tcga analysis")  # tcga
data <- uhlen_g4  # change accordingly
# plate_raw <- read.csv('all_tumor_g4_avg_ok.txt', header = TRUE, sep = "\t")
plate_raw <- read.csv('uhlen_g4_ok.txt', header = TRUE, sep = "\t")
plate <- as.data.frame(plate_raw)

sd_res <- data.frame(matrix(ncol = 0, nrow = 0))
max_res <- data.frame(matrix(ncol = 0, nrow = 0))
min_res <- data.frame(matrix(ncol = 0, nrow = 0))
mahal_res <- data.frame(matrix(ncol = 0, nrow = 0))
mean_res <- data.frame(matrix(ncol = 0, nrow = 0))
coef_res <- data.frame(matrix(ncol = 0, nrow = 0))

for (i in seq(1, nrow(data))){
  res1 <- sd(data[i,])
  sd_res <- rbind(sd_res, res1)
  res2 <- max(data[i,])
  max_res <- rbind(max_res, res2)
  res3 <- min(data[i,])
  min_res <- rbind(min_res, res3)
  res5 <- mean(as.numeric(data[i,]))
  mean_res <- rbind(mean_res, res5)
  res6 <- (sd(as.numeric(data[i,])) / mean(as.numeric(data[i,]))) *100
  coef_res <- rbind(coef_res, res6)
}
res4 <- as.data.frame(mahalanobis(data, colMeans(data), cov(data)))
mahal_res <- rbind(mahal_res, res4)

genes_sd <- as.data.frame(plate_raw[,2])  # for tcga, ideally use plate_raw[,1] (that's where the gene names are)
genes_sd <- cbind(genes_sd, sd_res)
genes_max <- as.data.frame(plate_raw[,2])
genes_max <- cbind(genes_max, max_res)
genes_min <- as.data.frame(plate_raw[,2])
genes_min <- cbind(genes_min, min_res)
genes_mahal <- as.data.frame(plate_raw[,2])
genes_mahal <- cbind(genes_mahal, mahal_res)
genes_mean <- as.data.frame(plate_raw[,2])
genes_mean <- cbind(genes_mean, mean_res)
genes_coef <- as.data.frame(plate_raw[,2])
genes_coef <- cbind(genes_coef, coef_res)

write.table(genes_sd, file = "sd_tumor_g4 - Copy.txt", sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(genes_max, file = "max_tumor_g4 - Copy.txt", sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(genes_min, file = "min_tumor_g4 - Copy.txt", sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(genes_mahal, file = "mahal_tumor_g4.txt", sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(genes_mean, file = "mean_g4.txt", sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(genes_coef, file = "coef_g4.txt", sep = "\t", row.names = FALSE, col.names = FALSE)

#--

# sd
uhlen_no_g4 <- read.table(paste0(path_uhlen, "sd_no_g4.txt"))
uhlen_g4 <- read.table(paste0(path_uhlen, "sd_g4.txt"))
tcga_normal_no_g4 <- read.table(paste0(path_tcga, "sd_normal_nog4.txt"))
tcga_normal_g4 <- read.table(paste0(path_tcga, "sd_normal_g4.txt"))
tcga_tumor_no_g4 <- read.table(paste0(path_tcga, "sd_tumor_nog4.txt"))
tcga_tumor_g4 <- read.table(paste0(path_tcga, "sd_tumor_g4.txt"))
# tcga_tumor_no_g4 <- read.table(paste0(path_tcga, "sd_tumor_nog4 - Copy.txt"))
# tcga_tumor_g4 <- read.table(paste0(path_tcga, "sd_tumor_g4 - Copy.txt"))

# max
uhlen_no_g4 <- read.table(paste0(path_uhlen, "max_no_g4.txt"))
uhlen_g4 <- read.table(paste0(path_uhlen, "max_g4.txt"))
tcga_normal_no_g4 <- read.table(paste0(path_tcga, "max_normal_nog4.txt"))
tcga_normal_g4 <- read.table(paste0(path_tcga, "max_normal_g4.txt"))
tcga_tumor_no_g4 <- read.table(paste0(path_tcga, "max_tumor_nog4.txt"))
tcga_tumor_g4 <- read.table(paste0(path_tcga, "max_tumor_g4.txt"))

# min
uhlen_no_g4 <- read.table(paste0(path_uhlen, "min_no_g4.txt"))
uhlen_g4 <- read.table(paste0(path_uhlen, "min_g4.txt"))
tcga_normal_no_g4 <- read.table(paste0(path_tcga, "min_normal_nog4.txt"))
tcga_normal_g4 <- read.table(paste0(path_tcga, "min_normal_g4.txt"))
tcga_tumor_no_g4 <- read.table(paste0(path_tcga, "min_tumor_nog4.txt"))
tcga_tumor_g4 <- read.table(paste0(path_tcga, "min_tumor_g4.txt"))

# mahalanobis distance
uhlen_no_g4 <- read.table(paste0(path_uhlen, "mahal_nog4.txt"))
uhlen_g4 <- read.table(paste0(path_uhlen, "mahal_g4.txt"))
tcga_normal_no_g4 <- read.table(paste0(path_tcga, "mahal_normal_nog4.txt"))
tcga_normal_g4 <- read.table(paste0(path_tcga, "mahal_normal_g4.txt"))
tcga_tumor_no_g4 <- read.table(paste0(path_tcga, "mahal_tumor_nog4.txt"))
tcga_tumor_g4 <- read.table(paste0(path_tcga, "mahal_tumor_g4.txt"))
# tcga_tumor_no_g4 <- read.table(paste0(path_tcga, "mahal_tumor_nog4 - Copy.txt"))
# tcga_tumor_g4 <- read.table(paste0(path_tcga, "mahal_tumor_g4 - Copy.txt"))

# mean
uhlen_no_g4 <- read.table(paste0(path_uhlen, "mean_no_g4.txt"))
uhlen_g4 <- read.table(paste0(path_uhlen, "mean_g4.txt"))
tcga_normal_no_g4 <- read.table(paste0(path_tcga, "mean_normal_nog4.txt"))
tcga_normal_g4 <- read.table(paste0(path_tcga, "mean_normal_g4.txt"))
tcga_tumor_no_g4 <- read.table(paste0(path_tcga, "mean_tumor_nog4.txt"))
tcga_tumor_g4 <- read.table(paste0(path_tcga, "mean_tumor_g4.txt"))

# coef
uhlen_no_g4 <- read.table(paste0(path_uhlen, "coef_no_g4.txt"))
uhlen_g4 <- read.table(paste0(path_uhlen, "coef_g4.txt"))
tcga_normal_no_g4 <- read.table(paste0(path_tcga, "coef_normal_nog4.txt"))
tcga_normal_g4 <- read.table(paste0(path_tcga, "coef_normal_g4.txt"))
tcga_tumor_no_g4 <- read.table(paste0(path_tcga, "coef_tumor_nog4.txt"))
tcga_tumor_g4 <- read.table(paste0(path_tcga, "coef_tumor_g4.txt"))

# OPTIONAL - Subsetting
# uhlen_g4 <- as.data.frame(uhlen_g4)
test <- uhlen_g4[1:4000,]
boxplot(uhlen_no_g4[,2], test[,2])

boxplot(uhlen_no_g4[,2], uhlen_g4[,2])
boxplot(tcga_normal_no_g4[,2], tcga_normal_g4[,2])
boxplot(tcga_tumor_no_g4[,2], tcga_tumor_g4[,2])

## Fancier boxplot ----
library(ggplot2)

uhlen_no_g4$gene_type <- "no_g4"
uhlen_g4$gene_type <- "g4"
myData <- rbind(uhlen_no_g4, uhlen_g4)

plot <- ggplot(myData, aes(x=gene_type, y=V2, fill = gene_type)) + 
  geom_boxplot(width = 0.45) +
  theme_classic() +
  #     geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.1) +
  #     ggtitle("Expression of genes that start at G4 DNA") +
  ylab("(SD/mean) *100") +
  # scale_x_discrete(labels = xlabs) +
  theme(plot.title = element_text(family = "Times New Roman", face = "bold.italic", color = "black", size = 8, hjust = 0.5),
        axis.title.x = element_text(family = "Helvetica", face = "plain", color = "black", size = 8),
        axis.title.y = element_text(family = "Helvetica", face = "plain", color = "black", size = 8),
        axis.text.x = element_text(family = "Helvetica", face = "italic", color = "black", size = 7, angle = 45),
        axis.text.y = element_text(family = "Helvetica", face = "plain", color = "black", size = 8),
        legend.text = element_text(family = "Helvetica", face = "bold", color = "black", size = 7),
        legend.title = element_text(family = "Helvetica", face = "italic", color = "black", size = 7)) #+
  # stat_summary(fun.y = mean, geom = "point", shape = 23, size = 2, col = "black", fill = "black") #+ 
#stat_compare_means(label.x = 1.1, label.y = 15.0, size = 2.7, color = "red")

plot


wilcox.test(uhlen_no_g4[,2], uhlen_g4[,2])  # W = 22277541, p-value < 2.2e-16
wilcox.test(tcga_normal_no_g4[,2], tcga_normal_g4[,2])  # W = 23605062, p-value < 2.2e-16
wilcox.test(tcga_tumor_no_g4[,2], tcga_tumor_g4[,2])  # W = 23775802, p-value < 2.2e-16 / W = 26752633, p-value < 2.2e-16 (after fixing no g4 list)
wilcox.test(tcga_normal_g4[,2], tcga_tumor_g4[,2])  # W = 140387153, p-value = 0.001479
wilcox.test(tcga_normal_no_g4[,2], tcga_tumor_no_g4[,2])  # W = 7572270, p-value = 0.4255
t.test(uhlen_no_g4[,2], uhlen_g4[,2])
t.test(tcga_normal_no_g4[,2], tcga_normal_g4[,2])
t.test(tcga_tumor_no_g4[,2], tcga_tumor_g4[,2])
t.test(tcga_normal_g4[,2], tcga_tumor_g4[,2])
t.test(tcga_normal_no_g4[,2], tcga_tumor_no_g4[,2])

#--

tcga_normal_no_g4 <- read.table(paste0(path_tcga, "all_normal_nog4_avg_ok.txt"), header=T)
tcga_normal_no_g4$gene_type <- "no_g4"
tcga_normal_g4 <- read.table(paste0(path_tcga, "all_normal_g4_avg_ok.txt"), header=T)
tcga_normal_g4$gene_type <- "g4"
tcga_tumor_no_g4 <- read.table(paste0(path_tcga, "all_tumor_nog4_avg_ok.txt"), header=T)
tcga_tumor_no_g4$gene_type <- "no_g4"
tcga_tumor_g4 <- read.table(paste0(path_tcga, "all_tumor_g4_avg_ok.txt"), header=T)
tcga_tumor_g4$gene_type <- "g4"

normal_genes <- rbind(tcga_normal_no_g4, tcga_normal_g4)
tumor_genes <- rbind(tcga_tumor_no_g4, tcga_tumor_g4)
# remove duplicated genes (there are 32 duplicated)
normal_genes <- subset(normal_genes, !duplicated(normal_genes$gene))
tumor_genes <- subset(tumor_genes, !duplicated(tumor_genes$gene))

sd_res <- data.frame(matrix(ncol = 4, nrow = nrow(normal_genes)))
rownames(sd_res) <- normal_genes$gene
colnames(sd_res) <- c("var_normal","var_tumor","var_tumor_minus_normal","gene_type")

for (i in seq(1, nrow(sd_res))){
  sd_res[i,"var_normal"] <- sd(normal_genes[i,2:16])
  sd_res[i,"var_tumor"] <- sd(tumor_genes[i,2:16])
  sd_res[i,"var_tumor_minus_normal"] <- sd_res[i,"var_tumor"] - sd_res[i,"var_normal"]
  sd_res[i,"gene_type"] <- normal_genes[i,17]
}

write.table(sd_res, "sd_res_normal_vs_tumor_all_genes.txt")

# Subsetting for only g4 genes or only no g4 genes
backup <- sd_res
# sd_res <- backup
sd_res <- as.data.frame(sd_res)
sd_res <- sd_res[sd_res$gene_type %in% "g4",]

# % of negative var_tumor_minus_normal
(nrow(sd_res[sd_res$var_tumor_minus_normal <0,]) / nrow(sd_res))*100

# % of positive var_tumor_minus_normal
(nrow(sd_res[sd_res$var_tumor_minus_normal >0,]) / nrow(sd_res))*100

# % of 0 var_tumor_minus_normal
(nrow(sd_res[sd_res$var_tumor_minus_normal %in% 0,]) / nrow(sd_res))*100

#---

# PCA plots ----
setwd("C:/Users/ruthb/Documents/UTH-MDA/Master's files/2020-1/RAW DATA/g4_known-transc_tcga_expression")

library(tidyr)

cancers_of_interest <- c("LIHC","HNSC","ESCA","STAD","LUAD","LUSC","THCA","UCEC","BRCA","BLCA","PRAD","COAD","KICH","KIRC","KIRP")

# Rescue metadata for normal tissue samples
metadata_normal <- as.data.frame(matrix(nrow=0, ncol=2))  # after loop, should have 681 samples
colnames(metadata_normal) <- c("sample_id", "tumor_type")

for (i in cancers_of_interest) {
  
  sample_labels <- as.data.frame(t(read.table(paste0("script_normal_samples/",i,"__geneExpN.txt"))[1,-c(1,2)]))
  colnames(sample_labels) <- "sample_id"
  sample_labels$tumor_type <- i
  metadata_normal <- rbind(metadata_normal, sample_labels)
  
}

# Rescue metadata for tumor samples
metadata_tumor <- as.data.frame(matrix(nrow=0, ncol=2))  # after loop, should have 6366 samples
colnames(metadata_tumor) <- c("sample_id", "tumor_type")

for (i in cancers_of_interest) {
  
  sample_labels <- as.data.frame(t(read.table(paste0("script_tumor_samples/",i,"__geneExpT.txt"))[1,-c(1,2)]))
  colnames(sample_labels) <- "sample_id"
  sample_labels$tumor_type <- i
  metadata_tumor <- rbind(metadata_tumor, sample_labels)
  
}

# Check if there are sample_ids in common among normal and tumor samples
test <- merge(metadata_normal, metadata_tumor, by="sample_id")  # no ids in common

# Only use samples that have correspondence for both normal and tumor
metadata_normal_sep <- separate(metadata_normal, sample_id, into = c("column1", "column2", "subject_id", "column4", "column5", "column6", "column7"), sep = "-")
metadata_normal$subject_id <- metadata_normal_sep$subject_id
colnames(metadata_normal)[1] <- "sample_id_normal"
metadata_tumor_sep <- separate(metadata_tumor, sample_id, into = c("column1", "column2", "subject_id", "column4", "column5", "column6", "column7"), sep = "-")
metadata_tumor$subject_id <- metadata_tumor_sep$subject_id
colnames(metadata_tumor)[1] <- "sample_id_tumor"

metadata <- merge(metadata_normal, metadata_tumor, by="subject_id")  # 644 samples overlap
metadata$sample_id_normal_overlap <- gsub("-", "\\.", metadata$sample_id_normal)
metadata$sample_id_tumor_overlap <- gsub("-", "\\.", metadata$sample_id_tumor)

write.table(metadata, "organized_metadata_paired_samples.txt")

# Concatenate files from different tissues, and only maintain data for samples of interest
normal_matrix <- as.data.frame(matrix(nrow=20530, ncol=2))  # after loop, should have 681 samples
normal_matrix[,1] <- read.table("script_normal_samples/LIHC__geneExpN.txt")[-1,1]  # gene name
normal_matrix[,2] <- read.table("script_normal_samples/LIHC__geneExpN.txt")[-1,2]  # entrez id
colnames(normal_matrix) <- c("gene_name", "entrez_id")

for (i in cancers_of_interest) {
  
  temp_matrix <- as.data.frame(read.table(paste0("script_normal_samples/",i,"__geneExpN.txt"), header=T))[,-c(1,2)]
  selected_columns <- metadata$sample_id_normal_overlap
  valid_columns <- intersect(selected_columns, colnames(temp_matrix))
  temp_matrix <- temp_matrix[, valid_columns]
  normal_matrix <- cbind(normal_matrix, temp_matrix)
  
}

write.table(normal_matrix, "normal_matrix_subset.txt")

tumor_matrix <- as.data.frame(matrix(nrow=20530, ncol=1))  # after loop, should have 681 samples
tumor_matrix[,1] <- read.table("script_tumor_samples/LIHC__geneExpT.txt")[-1,1]  # gene name
tumor_matrix[,2] <- read.table("script_tumor_samples/LIHC__geneExpT.txt")[-1,2]  # entrez id
colnames(tumor_matrix) <- c("gene_name", "entrez_id")

for (i in cancers_of_interest) {
  
  temp_matrix <- as.data.frame(read.table(paste0("script_tumor_samples/",i,"__geneExpT.txt"), header=T))[,-c(1,2)]
  selected_columns <- metadata$sample_id_tumor_overlap
  valid_columns <- intersect(selected_columns, colnames(temp_matrix))
  temp_matrix <- temp_matrix[, valid_columns]
  tumor_matrix <- cbind(tumor_matrix, temp_matrix)
  
}

write.table(tumor_matrix, "tumor_matrix_subset.txt")

# Separate concatenated matrices by pG4 and no pG4 genes
genes_without_g4 <- read.table("C:/Users/ruthb/Documents/UTH-MDA/Re-running g4 at genes and Uhlen graphs/ULTIMATE_GRAPH/genes_without_g4_genes_list_converted_ids.csv", 
                               sep=",", header=T)
genes_without_g4 <- unique(genes_without_g4$converted_alias)
normal_matrix <- read.table("normal_matrix_subset.txt")
tumor_matrix <- read.table("tumor_matrix_subset.txt")
rownames(normal_matrix) <- normal_matrix$entrez_id
rownames(tumor_matrix) <- tumor_matrix$entrez_id

# temp <- normal_matrix
normal_g4 <- t(normal_matrix[!(normal_matrix$entrez_id %in% genes_without_g4), ][,-c(1,2)])
normal_noG4 <- t(normal_matrix[normal_matrix$entrez_id %in% genes_without_g4, ][,-c(1,2)])
tumor_g4 <- t(tumor_matrix[!(tumor_matrix$entrez_id %in% genes_without_g4), ][,-c(1,2)])
tumor_noG4 <- t(tumor_matrix[tumor_matrix$entrez_id %in% genes_without_g4, ][,-c(1,2)])
write.table(normal_g4, "normal_g4")
write.table(normal_noG4, "normal_noG4")
write.table(tumor_g4, "tumor_g4")
write.table(tumor_noG4, "tumor_noG4")

### PCA plots with mean of each subject (inter-individual variability) ----
# Load required libraries
library(ggplot2)
library(dplyr)
library(stats)

# Identify constant/zero-variance columns
# constant_cols <- apply(normal_noG4, 2, function(col) var(col) == 0)
constant_cols <- apply(tumor_noG4, 2, function(col) var(col) == 0)

# Remove constant/zero-variance columns
# mat_filtered <- normal_noG4[, !constant_cols]
mat_filtered <- tumor_noG4[, !constant_cols]

# Perform PCA on the filtered matrix
pca_result <- prcomp(mat_filtered, scale. = TRUE)

# Extract the scores of the first two principal components
pca_scores <- as.data.frame(pca_result$x[, 1:2])
pca_scores$sample_id_normal_overlap <- rownames(pca_scores)  # for normal
pca_scores$sample_id_tumor_overlap <- rownames(pca_scores)  # for tumor

# Combine PCA scores with metadata
pca_data <- merge(pca_scores, metadata[,c(3,6)], by="sample_id_normal_overlap")  # for normal
pca_data <- merge(pca_scores, metadata[,c(5,7)], by="sample_id_tumor_overlap")  # for tumor

# OPTIONAL - Select specific cancers for display
# pca_data <- pca_data[pca_data[,4] %in% c("ESCA","HNSC","STAD"),]
pca_data <- pca_data[pca_data[,4] %in% c("HNSC"),]

# Create a scatterplot with labeled points
# ggplot(pca_data, aes(x = PC1, y = PC2, color = pca_data$tumor_type.x)) +
# ggplot(pca_data, aes(x = PC1, y = PC2, color = pca_data$tumor_type.y)) +
ggplot(pca_data, aes(x = PC1, y = PC2)) +
  geom_point(color = "mediumblue") +
  labs(x = paste("Principal Component 1 (", round(pca_result$sdev[1]^2 / sum(pca_result$sdev^2) * 100), "% variance)", sep = ""),
       y = paste("Principal Component 2 (", round(pca_result$sdev[2]^2 / sum(pca_result$sdev^2) * 100), "% variance)", sep = ""),
       title = "PCA Scatterplot") +
  theme_minimal() +
  scale_x_continuous(limits = c(-125, 125)) +
  scale_y_continuous(limits = c(-125, 125))


### Comparing variabilities using var.test ----
# Load required library
library(stats)

View(normal_g4)
View(normal_noG4)
View(tumor_g4)
View(tumor_noG4)

# Simulated PC1-PC2 data for illustration
# data1 <- matrix(rnorm(100, mean = 0, sd = 1), ncol = 2)
# data2 <- matrix(rnorm(100, mean = 0, sd = 2), ncol = 2)

# Identify constant/zero-variance columns
# constant_cols_1 <- apply(normal_g4, 2, function(col) var(col) == 0)
# constant_cols_2 <- apply(normal_noG4, 2, function(col) var(col) == 0)
constant_cols_1 <- apply(tumor_g4, 2, function(col) var(col) == 0)
constant_cols_2 <- apply(tumor_noG4, 2, function(col) var(col) == 0)

# Remove constant/zero-variance columns
# mat_filtered_1 <- normal_g4[, !constant_cols_1]
# mat_filtered_2 <- normal_noG4[, !constant_cols_2]
mat_filtered_1 <- tumor_g4[, !constant_cols_1]
mat_filtered_2 <- tumor_noG4[, !constant_cols_2]

# Perform PCA on the filtered matrix
pca_result_1 <- prcomp(mat_filtered_1, scale. = TRUE)
pca_result_2 <- prcomp(mat_filtered_2, scale. = TRUE)

# OPTIONAL - Find % of variance explained
pca_result <- pca_result_1
pc_std_dev <- pca_result$sdev

variance_explained <- (pc_std_dev^2) / sum(pc_std_dev^2)

percentage_variance_explained <- variance_explained * 100

View(as.data.frame(percentage_variance_explained))

# Extract the scores of the first two principal components
# pca_scores_1 <- as.data.frame(pca_result_1$x[, 1:2])
pca_scores_1 <- as.data.frame(pca_result_1$x)
pca_scores_1$sample_id_normal_overlap <- rownames(pca_scores_1)  # for normal
# pca_scores_1$sample_id_tumor_overlap <- rownames(pca_scores_1)  # for tumor
# pca_scores_2 <- as.data.frame(pca_result_2$x[, 1:2])
pca_scores_2 <- as.data.frame(pca_result_2$x)
pca_scores_2$sample_id_normal_overlap <- rownames(pca_scores_2)  # for normal
# pca_scores_2$sample_id_tumor_overlap <- rownames(pca_scores_2)  # for tumor

# Combine PCA scores with metadata
pca_data_1 <- merge(pca_scores_1, metadata[,c(3,6)], by="sample_id_normal_overlap")  # for normal
# pca_data_1 <- merge(pca_scores_1, metadata[,c(5,7)], by="sample_id_tumor_overlap")  # for tumor
pca_data_2 <- merge(pca_scores_2, metadata[,c(3,6)], by="sample_id_normal_overlap")  # for normal
# pca_data_2 <- merge(pca_scores_2, metadata[,c(5,7)], by="sample_id_tumor_overlap")  # for tumor

# OPTIONAL - Select specific cancers for display
# pca_data <- pca_data[pca_data[,4] %in% c("ESCA","HNSC","STAD"),]
pca_data_1 <- pca_data_1[pca_data_1[,4] %in% c("STAD"),]
pca_data_2 <- pca_data_2[pca_data_2[,4] %in% c("STAD"),]

# OPTIONAL - If using all data instead of PCs
pca_data_1 <- tumor_g4
pca_data_2 <- tumor_noG4
# View(normal_g4)
# View(normal_noG4)
# View(tumor_g4)
# View(tumor_noG4)

# Calculate variances
# variances_1 <- apply(pca_data_1[,2:3], 2, var)
# variances_2 <- apply(pca_data_2[,2:3], 2, var)
variances_1 <- apply(pca_data_1, 2, var)
variances_2 <- apply(pca_data_2, 2, var)

# Perform F-test
f_test_result <- var.test(variances_1, variances_2)
# f_test_results <- var.test(pca_data_1, pca_data_2)

# Access p-value
p_value <- f_test_result$p.value
p_value

# Interpret results
# if (p_value < 0.05) {
#   cat("The variability between the datasets is significantly different.\n")
# } else {
#   cat("The variability between the datasets is not significantly different.\n")
# }


### Comparing variabilities using sd + t-test ----
# Identify constant/zero-variance columns
constant_cols_1 <- apply(normal_g4, 2, function(col) var(col) == 0)
constant_cols_2 <- apply(normal_noG4, 2, function(col) var(col) == 0)
constant_cols_3 <- apply(tumor_g4, 2, function(col) var(col) == 0)
constant_cols_4 <- apply(tumor_noG4, 2, function(col) var(col) == 0)

# Remove constant/zero-variance columns
mat_filtered_1 <- normal_g4[, !constant_cols_1]
mat_filtered_2 <- normal_noG4[, !constant_cols_2]
mat_filtered_3 <- tumor_g4[, !constant_cols_3]
mat_filtered_4 <- tumor_noG4[, !constant_cols_4]
dim(mat_filtered_1)  # 644 17310
dim(mat_filtered_2)  # 644 2921
dim(mat_filtered_3)  # 644 17333
dim(mat_filtered_4)  # 644 2935

# Calculate sd per tissue
sd_res <- data.frame(matrix(ncol = 5, nrow = nrow(mat_filtered_1)))
rownames(sd_res) <- metadata$subject_id
rownames(metadata) <- metadata$subject_id
colnames(sd_res) <- c("var_normal_g4","var_normal_noG4","var_tumor_g4","var_tumor_noG4","cancer_type")

for (i in seq(1, nrow(sd_res))){
  id_normal <- metadata[rownames(sd_res)[i],"sample_id_normal_overlap"]
  id_tumor <- metadata[rownames(sd_res)[i],"sample_id_tumor_overlap"]
  sd_res[i,"var_normal_g4"] <- sd(mat_filtered_1[id_normal,])
  sd_res[i,"var_normal_noG4"] <- sd(mat_filtered_2[id_normal,])
  sd_res[i,"var_tumor_g4"] <- sd(mat_filtered_3[id_tumor,])
  sd_res[i,"var_tumor_noG4"] <- sd(mat_filtered_4[id_tumor,])
  sd_res[i,"cancer_type"] <- metadata[rownames(sd_res)[i],"tumor_type.y"]
}

# Sample sizes per cancer
cancer_counts <- table(sd_res$cancer_type)  # for UCEC, only 7 people

# Calculate Wilcoxon for each cancer
# wilcox_res <- data.frame(matrix(ncol = 1, nrow = nrow(cancer_counts)))
# rownames(wilcox_res) <- as.data.frame(cancer_counts)$Var1
# colnames(wilcox_res) <- "wilcox.p"
cancer_counts <- as.data.frame(cancer_counts)
rownames(cancer_counts) <- cancer_counts$Var1
colnames(cancer_counts) <- c("cancer_type","sample_size")
cancer_counts$t_test.p_NormG4_NormNoG4 <- NA
cancer_counts$t_test.p_TumG4_TumNoG4 <- NA
cancer_counts$t_test.p_NormG4_TumG4 <- NA
cancer_counts$t_test.p_NormNoG4_TumNoG4 <- NA

options(scipen = 999)  # decimals
# options(scipen = 0)  # scientific notation

for (i in seq(1, nrow(cancer_counts))){
  sd_res_subset <- sd_res[sd_res$cancer_type %in% rownames(cancer_counts)[i],]
  cancer_counts[i,"t_test.p_NormG4_NormNoG4"] <- t.test(sd_res_subset[,"var_normal_g4"],sd_res_subset[,"var_normal_noG4"])$p.value
  cancer_counts[i,"t_test.p_TumG4_TumNoG4"] <- t.test(sd_res_subset[,"var_tumor_g4"],sd_res_subset[,"var_tumor_noG4"])$p.value
  cancer_counts[i,"t_test.p_NormG4_TumG4"] <- t.test(sd_res_subset[,"var_normal_g4"],sd_res_subset[,"var_tumor_g4"])$p.value
  cancer_counts[i,"t_test.p_NormNoG4_TumNoG4"] <- t.test(sd_res_subset[,"var_normal_noG4"],sd_res_subset[,"var_tumor_noG4"])$p.value
}

write.table(cancer_counts, "summary_sd_t-tests_TCGA_matched_subjects.txt", row.names = F, col.names = T, sep="\t")
write.table(sd_res, "summary_sd_results_TCGA_matched_subjects.txt")

#---

### PCA plots with mean of each tissue ----
df <- as.data.frame(tumor_noG4)
df$sample_id_normal_overlap <- rownames(df)  # for normal
df$sample_id_tumor_overlap <- rownames(df)  # for tumor
df <- merge(metadata[,c(3,6)], df, by="sample_id_normal_overlap")  # for normal
df <- merge(metadata[,c(3,7)], df, by="sample_id_tumor_overlap")  # for tumor

mean_df <- as.data.frame(matrix(nrow=15, ncol=(ncol(df)-2)))
rownames(mean_df) <- cancers_of_interest
colnames(mean_df) <- colnames(df)[3:ncol(df)]

for (i in cancers_of_interest) {
  
  temp_df <- df[df$tumor_type.x %in% i,][,-c(1,2)]
  ncol(temp_df)
  temp_df <- t(as.data.frame(colMeans(temp_df)))
  ncol(temp_df)
  mean_df[i,] <- temp_df
  
}

mean_df <- as.matrix(mean_df)

# Load required libraries
library(ggplot2)
library(dplyr)
library(stats)

# Identify constant/zero-variance columns
constant_cols <- apply(mean_df, 2, function(col) var(col) == 0)

# Remove constant/zero-variance columns
mat_filtered <- mean_df[, !constant_cols]

# Perform PCA on the filtered matrix
pca_result <- prcomp(mat_filtered, scale. = TRUE)

# Extract the scores of the first two principal components
pca_scores <- as.data.frame(pca_result$x[, 1:2])
# pca_scores$tumor_type <- rownames(pca_scores)

# # Combine PCA scores with metadata
# pca_data <- merge(pca_scores, metadata[,c(3,6)], by="sample_id_normal_overlap")  # for normal
# pca_data <- merge(pca_scores, metadata[,c(5,7)], by="sample_id_tumor_overlap")  # for tumor

# # OPTIONAL - Select specific cancers for display
# # pca_data <- pca_data[pca_data[,4] %in% c("ESCA","HNSC","STAD"),]
# pca_data <- pca_data[pca_data[,4] %in% c("STAD"),]
pca_data <- pca_scores

# Create a scatterplot with labeled points
ggplot(pca_data, aes(x = PC1, y = PC2, color = rownames(pca_data))) +
  # ggplot(pca_data, aes(x = PC1, y = PC2, color = pca_data$tumor_type.y)) +
  geom_point() +
  labs(x = paste("Principal Component 1 (", round(pca_result$sdev[1]^2 / sum(pca_result$sdev^2) * 100), "% variance)", sep = ""),
       y = paste("Principal Component 2 (", round(pca_result$sdev[2]^2 / sum(pca_result$sdev^2) * 100), "% variance)", sep = ""),
       title = "PCA Scatterplot") +
  theme_minimal() +
  scale_x_continuous(limits = c(-125, 125)) +
  scale_y_continuous(limits = c(-125, 125))

#---

# Boxplots comparing variability in G4 genes (tumor vs. normal) ----
setwd("C:/Users/ruthb/Documents/UTH-MDA/Master's files/2020-1/RAW DATA/g4_known-transc_tcga_expression/")

sd_res <- read.table("summary_sd_results_TCGA_matched_subjects.txt")
# cancer_counts <- read.table("summary_sd_t-tests_TCGA_matched_subjects.txt", header=T)

# Suggested by Albino: Calculate (sd / mean) *100 (because standard deviation is very sensitive to big means)
normal_g4 <- read.table("normal_g4")
normal_noG4 <- read.table("normal_noG4")
tumor_g4 <- read.table("tumor_g4")
tumor_noG4 <- read.table("tumor_noG4")
metadata <- read.table("organized_metadata_paired_samples.txt")

# Identify constant/zero-variance columns
constant_cols_1 <- apply(normal_g4, 2, function(col) var(col) == 0)
constant_cols_2 <- apply(normal_noG4, 2, function(col) var(col) == 0)
constant_cols_3 <- apply(tumor_g4, 2, function(col) var(col) == 0)
constant_cols_4 <- apply(tumor_noG4, 2, function(col) var(col) == 0)

# Remove constant/zero-variance columns
mat_filtered_1 <- normal_g4[, !constant_cols_1]
mat_filtered_2 <- normal_noG4[, !constant_cols_2]
mat_filtered_3 <- tumor_g4[, !constant_cols_3]
mat_filtered_4 <- tumor_noG4[, !constant_cols_4]
dim(mat_filtered_1)  # 644 17310
dim(mat_filtered_2)  # 644 2921
dim(mat_filtered_3)  # 644 17333
dim(mat_filtered_4)  # 644 2935

# Calculate mean per tissue
mean_res <- data.frame(matrix(ncol = 5, nrow = nrow(mat_filtered_1)))
rownames(mean_res) <- metadata$subject_id
rownames(metadata) <- metadata$subject_id
colnames(mean_res) <- c("mean_normal_g4","mean_normal_noG4","mean_tumor_g4","mean_tumor_noG4","cancer_type")

for (i in seq(1, nrow(mean_res))){
  id_normal <- metadata[rownames(mean_res)[i],"sample_id_normal_overlap"]
  id_tumor <- metadata[rownames(mean_res)[i],"sample_id_tumor_overlap"]
  mean_res[i,"mean_normal_g4"] <- mean(as.numeric(mat_filtered_1[id_normal,]))
  mean_res[i,"mean_normal_noG4"] <- mean(as.numeric(mat_filtered_2[id_normal,]))
  mean_res[i,"mean_tumor_g4"] <- mean(as.numeric(mat_filtered_3[id_tumor,]))
  mean_res[i,"mean_tumor_noG4"] <- mean(as.numeric(mat_filtered_4[id_tumor,]))
  mean_res[i,"cancer_type"] <- metadata[rownames(mean_res)[i],"tumor_type.y"]
}

# Calculate coefficient of variation per sample
coef_res <- data.frame(matrix(ncol = 5, nrow = nrow(mat_filtered_1)))
rownames(coef_res) <- metadata$subject_id
rownames(metadata) <- metadata$subject_id
colnames(coef_res) <- c("coef_normal_g4","coef_normal_noG4","coef_tumor_g4","coef_tumor_noG4","cancer_type")

for (i in seq(1, nrow(coef_res))){
  id_normal <- metadata[rownames(coef_res)[i],"sample_id_normal_overlap"]
  id_tumor <- metadata[rownames(coef_res)[i],"sample_id_tumor_overlap"]
  coef_res[i,"coef_normal_g4"] <- (sd(as.numeric(mat_filtered_1[id_normal,])) / mean(as.numeric(mat_filtered_1[id_normal,]))) *100
  coef_res[i,"coef_normal_noG4"] <- (sd(as.numeric(mat_filtered_2[id_normal,])) / mean(as.numeric(mat_filtered_2[id_normal,]))) *100
  coef_res[i,"coef_tumor_g4"] <- (sd(as.numeric(mat_filtered_3[id_tumor,])) / mean(as.numeric(mat_filtered_3[id_tumor,]))) *100
  coef_res[i,"coef_tumor_noG4"] <- (sd(as.numeric(mat_filtered_4[id_tumor,])) / mean(as.numeric(mat_filtered_4[id_tumor,]))) *100
  coef_res[i,"cancer_type"] <- metadata[rownames(coef_res)[i],"tumor_type.y"]
}

# Sample sizes per cancer
cancer_counts <- table(mean_res$cancer_type)  # for UCEC, only 7 people

# Calculate Wilcoxon/T-test for each cancer
# wilcox_res <- data.frame(matrix(ncol = 1, nrow = nrow(cancer_counts)))
# rownames(wilcox_res) <- as.data.frame(cancer_counts)$Var1
# colnames(wilcox_res) <- "wilcox.p"
cancer_counts <- as.data.frame(cancer_counts)
rownames(cancer_counts) <- cancer_counts$Var1
colnames(cancer_counts) <- c("cancer_type","sample_size")
cancer_counts$t_test.p_NormG4_NormNoG4 <- NA
cancer_counts$t_test.p_TumG4_TumNoG4 <- NA
cancer_counts$t_test.p_NormG4_TumG4 <- NA
cancer_counts$t_test.p_NormNoG4_TumNoG4 <- NA

options(scipen = 999)  # decimals
# options(scipen = 0)  # scientific notation

## Test - paired t-test
coef_res <- read.table("C:/Users/ruthb/Documents/UTH-MDA/Master's files/2020-1/RAW DATA/g4_known-transc_tcga_expression/summary_coefficient_of_variation_results_TCGA_matched_subjects.txt")
i=1
coef_res_subset <- coef_res[coef_res$cancer_type %in% rownames(cancer_counts)[i],]
t.test(coef_res_subset[,"coef_normal_g4"],coef_res_subset[,"coef_normal_noG4"],paired=T)$p.value
##

for (i in seq(1, nrow(cancer_counts))){
  coef_res_subset <- coef_res[coef_res$cancer_type %in% rownames(cancer_counts)[i],]
  # cancer_counts[i,"t_test.p_NormG4_NormNoG4"] <- t.test(coef_res_subset[,"coef_normal_g4"],coef_res_subset[,"coef_normal_noG4"])$p.value
  # cancer_counts[i,"t_test.p_TumG4_TumNoG4"] <- t.test(coef_res_subset[,"coef_tumor_g4"],coef_res_subset[,"coef_tumor_noG4"])$p.value
  # cancer_counts[i,"t_test.p_NormG4_TumG4"] <- t.test(coef_res_subset[,"coef_normal_g4"],coef_res_subset[,"coef_tumor_g4"])$p.value
  # cancer_counts[i,"t_test.p_NormNoG4_TumNoG4"] <- t.test(coef_res_subset[,"coef_normal_noG4"],coef_res_subset[,"coef_tumor_noG4"])$p.value
  cancer_counts[i,"t_test.p_NormG4_NormNoG4"] <- t.test(coef_res_subset[,"coef_normal_g4"],coef_res_subset[,"coef_normal_noG4"],paired=T)$p.value
  cancer_counts[i,"t_test.p_TumG4_TumNoG4"] <- t.test(coef_res_subset[,"coef_tumor_g4"],coef_res_subset[,"coef_tumor_noG4"],paired=T)$p.value
  cancer_counts[i,"t_test.p_NormG4_TumG4"] <- t.test(coef_res_subset[,"coef_normal_g4"],coef_res_subset[,"coef_tumor_g4"],paired=T)$p.value
  cancer_counts[i,"t_test.p_NormNoG4_TumNoG4"] <- t.test(coef_res_subset[,"coef_normal_noG4"],coef_res_subset[,"coef_tumor_noG4"],paired=T)$p.value
}

write.table(cancer_counts, "summary_coef_PAIRED_t-tests_TCGA_matched_subjects.txt", row.names = F, col.names = T, sep="\t")
write.table(mean_res, "summary_mean_results_TCGA_matched_subjects.txt")
write.table(coef_res, "summary_coefficient_of_variation_results_TCGA_matched_subjects.txt")


# *Plotting boxplots ----
setwd("C:/Users/ruthb/Documents/UTH-MDA/Master's files/2020-1/RAW DATA/g4_known-transc_tcga_expression/")
library(ggplot2)

# selected_cancers <- c("LUAD","LUSC","HNSC","BLCA","ESCA","COAD","KIRP")
selected_cancers <- c("LUSC","HNSC","BLCA","KIRP","COAD")

# If using coefficient of variation
coef_res <- read.table("summary_coefficient_of_variation_results_TCGA_matched_subjects.txt")
sd_res <- coef_res
colnames(sd_res) <- c("var_normal_g4","var_normal_noG4","var_tumor_g4","var_tumor_noG4","cancer_type")

cancer_name <- selected_cancers[1]  # change accordingly
normal <- sd_res[sd_res$cancer_type %in% cancer_name,][,c("var_normal_g4","cancer_type")]
tumor <- sd_res[sd_res$cancer_type %in% cancer_name,][,c("var_tumor_g4","cancer_type")]
# normal <- sd_res[sd_res$cancer_type %in% cancer_name,][,c("var_normal_noG4","cancer_type")]
# tumor <- sd_res[sd_res$cancer_type %in% cancer_name,][,c("var_tumor_noG4","cancer_type")]
normal$sample_type <- "Normal"
tumor$sample_type <- "Tumor"
colnames(normal)[1] <- "sd"
colnames(tumor)[1] <- "sd"
normal$sample_id <- rownames(normal)
tumor$sample_id <- rownames(tumor)
myData <- rbind(normal, tumor)
# OPTIONAL - log transformation (important for ESCA)
# myData$sd <- log(myData$sd)
# # For ESCA only - Remove A43C and A4OR
# myData_0 <- myData
# myData <- myData[!(rownames(myData) %in% c("A43C","A4OR")),]

plot <- ggplot(myData, aes(x=sample_type, y=sd, fill = sample_type)) + 
  geom_boxplot(width = 0.45) +
  theme_classic() +
  #     geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.1) +
  #     ggtitle("Expression of genes that start at G4 DNA") +
  ylab("(SD/mean) *100") +
  # scale_x_discrete(labels = xlabs) +
  theme(plot.title = element_text(face = "bold.italic", color = "black", size = 8, hjust = 0.5)) +
        # axis.title.x = element_text(family = "Helvetica", face = "plain", color = "black", size = 8),
        # axis.title.y = element_text(family = "Helvetica", face = "plain", color = "black", size = 8),
        # axis.text.x = element_text(family = "Helvetica", face = "italic", color = "black", size = 7, angle = 45),
        # axis.text.y = element_text(family = "Helvetica", face = "plain", color = "black", size = 8),
        # legend.text = element_text(family = "Helvetica", face = "bold", color = "black", size = 7),
        # legend.title = element_text(family = "Helvetica", face = "italic", color = "black", size = 7)) +
  scale_y_continuous(breaks = seq(0, 2000, by = 500), limits = c(0, 2000)) #+
# stat_summary(fun.y = mean, geom = "point", shape = 23, size = 2, col = "black", fill = "black") #+ 
#stat_compare_means(label.x = 1.1, label.y = 15.0, size = 2.7, color = "red")

plot
# 3.03 x 3.84 portrait


# *Plotting paired line plots ----
# Sample data (replace this with your own paired data)
data <- merge(normal, tumor, by="sample_id")

# Create a paired line plot using ggplot2
data$Direction <- ifelse(data$sd.x < data$sd.y, "up", "down")
colnames(data)[2] <- "Normal"
colnames(data)[5] <- "Tumor"

# temp1 <- data[,1:4]
# colnames(temp1) <- c("sample_id","sd","cancer_type","sample_type")
# temp2 <- data[,c(1,5:7)]
# colnames(temp2) <- c("sample_id","sd","cancer_type","sample_type")
# myData2 <- rbind(temp1, temp2)

library(ggpubr)
color_palette <- c("up" = "red", "down" = "blue")
# color_palette2 <- c("Normal" = "blue", "Tumor" = "red")

# ggpaired: with color for up and down
ggpaired(data, cond1 = "Normal", cond2 = "Tumor",
  fill = "condition",
  #palette = "jcb",
  width = 0.5,
  color = "#333333", #point.size = 2,
  line.color = "Direction", geom = "line") +
  theme_classic() +
  # geom_boxplot(alpha = 0.2) +
  # theme(plot.title = element_text(face = "bold.italic", color = "black", size = 8, hjust = 0.5)) +
  scale_color_manual(values = color_palette) +
  scale_y_continuous(breaks = seq(0, 2000, by = 500), limits = c(0, 2000)) #+
  # theme(legend.position="none")
# 3.03 x 3.84 portrait

# # ggplot: without color for up and down
# ggplot(myData, aes(x = sample_type, y = sd)) +
#   # geom_boxplot(aes(fill = Condition), alpha = 0.2, col = "grey") +
#   geom_point(aes(col = sample_type)) +
#   geom_line(aes(group = sample_id))

# Analyzing n of increased and decreased
sd_res$inc_dec_G4 <- ifelse(sd_res$var_normal_g4 < sd_res$var_tumor_g4, "up", "down")
sd_res$inc_dec_noG4 <- ifelse(sd_res$var_normal_noG4 < sd_res$var_tumor_noG4, "up", "down")

direct <- "up"
canc <- "COAD"
# nrow(sd_res[sd_res$inc_dec_G4 %in% direct,][sd_res[sd_res$inc_dec_G4 %in% direct,]$cancer_type %in% canc,]) #&& sd_res$cancer_type %in% "LUSC"
nrow(sd_res[sd_res$inc_dec_noG4 %in% direct,][sd_res[sd_res$inc_dec_noG4 %in% direct,]$cancer_type %in% canc,]) #&& sd_res$cancer_type %in% "LUSC"

#----

# Inter-individual graphs for Uhlen ----
# Albino mentioned I should plot those for the healthy tissues of selected cancers (LUSC, HNSC, BLCA, KIRP, and COAD).
# lung, head and neck (salivarygland and tonsil), bladder (urinarybladder), kidney and colon (colon and colonrectum)
# Plot G4 vs. no G4. Also, when creating final figure panel, consider changing colors for normal vs. tumor comparisons

uhlen_g4 <- read.table("uhlen_g4.txt", header = T)
uhlen_no_g4 <- read.table("uhlen_no_g4.txt", header = T)
rownames(uhlen_g4) <- uhlen_g4$enstid
rownames(uhlen_no_g4) <- uhlen_no_g4$enstid

tissues <- as.data.frame(colnames(uhlen_g4)[-c(1,2)])
library(tidyr)
tissue_types <- separate(tissues, colnames(tissues), into = c("tissue_name", "id"), sep = "_")
tissue_types_2 <- as.data.frame(unique(tissue_types$tissue_name))  ## 33 different tissues

## Part 1 - 1 coef per gene ----
average_matrix_g4 <- as.data.frame(matrix(nrow=nrow(uhlen_g4), ncol=nrow(tissue_types_2)))
rownames(average_matrix_g4) <- uhlen_g4$enstid
colnames(average_matrix_g4) <- tissue_types_2[,1]

sd_matrix_g4 <- as.data.frame(matrix(nrow=nrow(uhlen_g4), ncol=nrow(tissue_types_2)))
rownames(sd_matrix_g4) <- uhlen_g4$enstid
colnames(sd_matrix_g4) <- tissue_types_2[,1]

coef_var_matrix_g4 <- as.data.frame(matrix(nrow=nrow(uhlen_g4), ncol=nrow(tissue_types_2)))
rownames(coef_var_matrix_g4) <- uhlen_g4$enstid
colnames(coef_var_matrix_g4) <- tissue_types_2[,1]

average_matrix_no_g4 <- as.data.frame(matrix(nrow=nrow(uhlen_no_g4), ncol=nrow(tissue_types_2)))
rownames(average_matrix_no_g4) <- uhlen_no_g4$enstid
colnames(average_matrix_no_g4) <- tissue_types_2[,1]

sd_matrix_no_g4 <- as.data.frame(matrix(nrow=nrow(uhlen_no_g4), ncol=nrow(tissue_types_2)))
rownames(sd_matrix_no_g4) <- uhlen_no_g4$enstid
colnames(sd_matrix_no_g4) <- tissue_types_2[,1]

coef_var_matrix_no_g4 <- as.data.frame(matrix(nrow=nrow(uhlen_no_g4), ncol=nrow(tissue_types_2)))
rownames(coef_var_matrix_no_g4) <- uhlen_no_g4$enstid
colnames(coef_var_matrix_no_g4) <- tissue_types_2[,1]

for (i in 1:length(tissue_types_2[,1])) {
  
  # Preparing subset for selected tissue
  tissue_temp <- tissue_types_2[i,1]
  df_subset_g4 <- colnames(uhlen_g4)[grep(tissue_temp, colnames(uhlen_g4))]
  df_subset_g4 <- uhlen_g4[, df_subset_g4]
  df_subset_no_g4 <- colnames(uhlen_no_g4)[grep(tissue_temp, colnames(uhlen_no_g4))]
  df_subset_no_g4 <- uhlen_no_g4[, df_subset_no_g4]
  
  # Calculating average per gene
  average_matrix_g4[,i] <- rowMeans(df_subset_g4)
  average_matrix_no_g4[,i] <- rowMeans(df_subset_no_g4)
  
  # Calculating sd per gene
  sd_matrix_g4[,i] <- apply(df_subset_g4, 1, sd)
  sd_matrix_no_g4[,i] <- apply(df_subset_no_g4, 1, sd)
  
  # Calculating coefficient of variation per gene
  coef_var_matrix_g4[,i] <- (sd_matrix_g4[,i] / average_matrix_g4[,i]) *100
  coef_var_matrix_no_g4[,i] <- (sd_matrix_no_g4[,i] / average_matrix_no_g4[,i]) *100
  
}

write.table(coef_var_matrix_g4, "summary_coefficient_of_variation_results_Uhlen_subjects_1_coef_per_gene_G4.txt")
write.table(coef_var_matrix_no_g4, "summary_coefficient_of_variation_results_Uhlen_subjects_1_coef_per_gene_no_G4.txt")

# Calculating Wilcoxon tests
wilcoxon_matrix <- as.data.frame(matrix(nrow=ncol(coef_var_matrix_g4), ncol=1))
rownames(wilcoxon_matrix) <- colnames(coef_var_matrix_g4)
colnames(wilcoxon_matrix) <- "wilcoxon.p"

for (i in 1:length(colnames(coef_var_matrix_g4))) {
  wilcoxon_matrix[i,1] <- wilcox.test(coef_var_matrix_g4[,i], coef_var_matrix_no_g4[,i])$p.value
}
# All results are significant

write.table(wilcoxon_matrix, "summary_wilcoxon_of_variation_results_Uhlen_subjects_1_coef_per_gene.txt")

number <- 33  # from 1 to 33
boxplot(coef_var_matrix_g4[,number], coef_var_matrix_no_g4[,number])
# In all cases variation is bigger for genes without G4 than for genes with G4


## Part 2 - 1 coef per individual ----
average_matrix <- as.data.frame(matrix(nrow=nrow(tissue_types), ncol=3))
rownames(average_matrix) <- tissues[,1]
colnames(average_matrix) <- c("pG4","no_pG4","tissue_name")
average_matrix$tissue_name <- tissue_types$tissue_name

sd_matrix <- as.data.frame(matrix(nrow=nrow(tissue_types), ncol=3))
rownames(sd_matrix) <- tissues[,1]
colnames(sd_matrix) <- c("pG4","no_pG4","tissue_name")
sd_matrix$tissue_name <- tissue_types$tissue_name

coef_var_matrix <- as.data.frame(matrix(nrow=nrow(tissue_types), ncol=3))
rownames(coef_var_matrix) <- tissues[,1]
colnames(coef_var_matrix) <- c("pG4","no_pG4","tissue_name")
coef_var_matrix$tissue_name <- tissue_types$tissue_name

# Calculating average per gene
average_matrix[,"pG4"] <- colMeans(uhlen_g4[,-c(1,2)])
average_matrix[,"no_pG4"] <- colMeans(uhlen_no_g4[,-c(1,2)])

# Calculating sd per gene
sd_matrix[,"pG4"] <- apply(uhlen_g4[,-c(1,2)], 2, sd)
sd_matrix[,"no_pG4"] <- apply(uhlen_no_g4[,-c(1,2)], 2, sd)

# Calculating coefficient of variation per gene
coef_var_matrix[,"pG4"] <- (sd_matrix[,"pG4"] / average_matrix[,"pG4"]) *100
coef_var_matrix[,"no_pG4"] <- (sd_matrix[,"no_pG4"] / average_matrix[,"no_pG4"]) *100


write.table(coef_var_matrix, "summary_coefficient_of_variation_results_Uhlen_subjects_1_coef_per_individual.txt")


# Plotting

library(ggplot2)

# selected_cancers <- c("LUAD","LUSC","HNSC","BLCA","ESCA","COAD","KIRP")
# selected_cancers <- c("LUSC","HNSC","BLCA","KIRP","COAD")
selected_cancers <- c("lung", "salivarygland", "tonsil", "urinarybladder", "kidney", "colon", "colonrectum")

# If using coefficient of variation
sd_res <- read.table("summary_coefficient_of_variation_results_Uhlen_subjects_1_coef_per_individual.txt")
# colnames(sd_res) <- c("var_normal_g4","var_normal_noG4","var_tumor_g4","var_tumor_noG4","cancer_type")

cancer_name <- selected_cancers[7]  # change accordingly

# Normal = g4; tumor = no_g4
normal <- sd_res[sd_res$tissue_name %in% cancer_name,][,c("pG4","tissue_name")]
tumor <- sd_res[sd_res$tissue_name %in% cancer_name,][,c("no_pG4","tissue_name")]
normal$sample_type <- "1 - With pG4"
tumor$sample_type <- "2 - No pG4"
colnames(normal)[1] <- "coef_var"
colnames(tumor)[1] <- "coef_var"
myData <- rbind(normal, tumor)
# OPTIONAL - log transformation (important for ESCA)
# myData$sd <- log(myData$sd)
# # For ESCA only - Remove A43C and A4OR
# myData_0 <- myData
# myData <- myData[!(rownames(myData) %in% c("A43C","A4OR")),]

plot <- ggplot(myData, aes(x=sample_type, y=coef_var, fill = sample_type)) + 
  geom_boxplot(width = 0.45) +
  theme_classic() +
  #     geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.1) +
  #     ggtitle("Expression of genes that start at G4 DNA") +
  ylab("(SD/mean) *100") +
  # scale_x_discrete(labels = xlabs) +
  theme(plot.title = element_text(face = "bold.italic", color = "black", size = 8, hjust = 0.5)) +
  # axis.title.x = element_text(family = "Helvetica", face = "plain", color = "black", size = 8),
  # axis.title.y = element_text(family = "Helvetica", face = "plain", color = "black", size = 8),
  # axis.text.x = element_text(family = "Helvetica", face = "italic", color = "black", size = 7, angle = 45),
  # axis.text.y = element_text(family = "Helvetica", face = "plain", color = "black", size = 8),
  # legend.text = element_text(family = "Helvetica", face = "bold", color = "black", size = 7),
  # legend.title = element_text(family = "Helvetica", face = "italic", color = "black", size = 7)) +
  scale_y_continuous(breaks = seq(0, 2000, by = 500), limits = c(0, 2000)) #+
# stat_summary(fun.y = mean, geom = "point", shape = 23, size = 2, col = "black", fill = "black") #+ 
#stat_compare_means(label.x = 1.1, label.y = 15.0, size = 2.7, color = "red")

plot
# 3.03 x 3.84

t.test(normal$coef_var, tumor$coef_var)$p.value
wilcox.test(normal$coef_var, tumor$coef_var)$p.value
# Good examples: "lung", "salivarygland", "kidney"

# Inter-individual graphs for TCGA ----
raw_files_path <- "C:/Users/ruthb/Documents/UTH-MDA/Master's files/2020-1/RAW DATA/g4_known-transc_tcga_expression"

# Loading matched samples
normal_g4 <- read.table(paste0(raw_files_path, "/normal_g4"))
normal_noG4 <- read.table(paste0(raw_files_path, "/normal_noG4"))
tumor_g4 <- read.table(paste0(raw_files_path, "/tumor_g4"))
tumor_noG4 <- read.table(paste0(raw_files_path, "/tumor_noG4"))

# Identify constant/zero-variance columns
constant_cols_1 <- apply(normal_g4, 2, function(col) var(col) == 0)
constant_cols_2 <- apply(normal_noG4, 2, function(col) var(col) == 0)
constant_cols_3 <- apply(tumor_g4, 2, function(col) var(col) == 0)
constant_cols_4 <- apply(tumor_noG4, 2, function(col) var(col) == 0)

# Remove constant/zero-variance columns
normal_g4 <- t(normal_g4[, !constant_cols_1])
normal_noG4 <- t(normal_noG4[, !constant_cols_2])
tumor_g4 <- t(tumor_g4[, !constant_cols_3])
tumor_noG4 <- t(tumor_noG4[, !constant_cols_4])

# Load metadata
metadata <- read.table(paste0(raw_files_path, "/organized_metadata_paired_samples.txt"))
tissues <- as.data.frame(unique(metadata$tumor_type.x))

## Part 1 - 1 coef per gene ----
matrix_normal_g4 <- as.data.frame(matrix(nrow=nrow(normal_g4), ncol=nrow(tissues)))
rownames(matrix_normal_g4) <- rownames(normal_g4)
colnames(matrix_normal_g4) <- tissues[,1]

matrix_normal_no_g4 <- as.data.frame(matrix(nrow=nrow(normal_noG4), ncol=nrow(tissues)))
rownames(matrix_normal_no_g4) <- rownames(normal_noG4)
colnames(matrix_normal_no_g4) <- tissues[,1]

matrix_tumor_g4 <- as.data.frame(matrix(nrow=nrow(tumor_g4), ncol=nrow(tissues)))
rownames(matrix_tumor_g4) <- rownames(tumor_g4)
colnames(matrix_tumor_g4) <- tissues[,1]

matrix_tumor_no_g4 <- as.data.frame(matrix(nrow=nrow(tumor_noG4), ncol=nrow(tissues)))
rownames(matrix_tumor_no_g4) <- rownames(tumor_noG4)
colnames(matrix_tumor_no_g4) <- tissues[,1]

for (i in 1:length(tissues[,1])) {
  
  # Preparing subset for selected tissue
  tissue_temp <- tissues[i,1]
  samples_temp_normal <- metadata[metadata$tumor_type.x %in% tissue_temp,]$sample_id_normal_overlap
  samples_temp_tumor <- metadata[metadata$tumor_type.x %in% tissue_temp,]$sample_id_tumor_overlap
  
  df_subset_normal_g4 <- normal_g4[,colnames(normal_g4) %in% samples_temp_normal]
  df_subset_normal_no_g4 <- normal_noG4[,colnames(normal_noG4) %in% samples_temp_normal]
  df_subset_tumor_g4 <- tumor_g4[,colnames(tumor_g4) %in% samples_temp_tumor]
  df_subset_tumor_no_g4 <- tumor_noG4[,colnames(tumor_noG4) %in% samples_temp_tumor]
  
  # Calculating average per gene
  average_1 <- as.data.frame(rowMeans(df_subset_normal_g4))
  average_2 <- as.data.frame(rowMeans(df_subset_normal_no_g4))
  average_3 <- as.data.frame(rowMeans(df_subset_tumor_g4))
  average_4 <- as.data.frame(rowMeans(df_subset_tumor_no_g4))
  
  # Calculating sd per gene
  sd_1 <- as.data.frame(apply(df_subset_normal_g4, 1, sd))
  sd_2 <- as.data.frame(apply(df_subset_normal_no_g4, 1, sd))
  sd_3 <- as.data.frame(apply(df_subset_tumor_g4, 1, sd))
  sd_4 <- as.data.frame(apply(df_subset_tumor_no_g4, 1, sd))
  
  # Calculating coefficient of variation per gene
  matrix_normal_g4[,i] <- (sd_1[,1] / average_1[,1]) *100
  matrix_normal_no_g4[,i] <- (sd_2[,1] / average_2[,1]) *100
  matrix_tumor_g4[,i] <- (sd_3[,1] / average_3[,1]) *100
  matrix_tumor_no_g4[,i] <- (sd_4[,1] / average_4[,1]) *100
  
}

write.table(matrix_normal_g4, paste0(path_tcga,"summary_coefficient_of_variation_results_TCGA_subjects_1_coef_per_gene_normal_G4.txt"))
write.table(matrix_normal_no_g4, paste0(path_tcga,"summary_coefficient_of_variation_results_TCGA_subjects_1_coef_per_gene_normal_no_G4.txt"))
write.table(matrix_tumor_g4, paste0(path_tcga,"summary_coefficient_of_variation_results_TCGA_subjects_1_coef_per_gene_tumor_G4.txt"))
write.table(matrix_tumor_no_g4, paste0(path_tcga,"summary_coefficient_of_variation_results_TCGA_subjects_1_coef_per_gene_tumor_no_G4.txt"))

# Calculating Wilcoxon tests
wilcoxon_matrix <- as.data.frame(matrix(nrow=nrow(tissues), ncol=4))
rownames(wilcoxon_matrix) <- tissues[,1]
colnames(wilcoxon_matrix) <- c("wilcox.p_NormG4_NormNoG4","wilcox.p_TumG4_TumNoG4","wilcox.p_NormG4_TumG4","wilcox.p_NormNoG4_TumNoG4")

for (i in 1:length(colnames(matrix_normal_g4))) {
  wilcoxon_matrix[i,1] <- wilcox.test(matrix_normal_g4[,i], matrix_normal_no_g4[,i])$p.value
  wilcoxon_matrix[i,2] <- wilcox.test(matrix_tumor_g4[,i], matrix_tumor_no_g4[,i])$p.value
  wilcoxon_matrix[i,3] <- wilcox.test(matrix_normal_g4[,i], matrix_tumor_g4[,i])$p.value
  wilcoxon_matrix[i,4] <- wilcox.test(matrix_normal_no_g4[,i], matrix_tumor_no_g4[,i])$p.value
}
# Most results are significant

write.table(wilcoxon_matrix, paste0(path_tcga,"summary_wilcoxon_of_variation_results_TCGA_subjects_1_coef_per_gene.txt"))

matrix_normal_g4 <- read.table(paste0(path_tcga,"summary_coefficient_of_variation_results_TCGA_subjects_1_coef_per_gene_normal_G4.txt"))
matrix_normal_no_g4 <- read.table(paste0(path_tcga,"summary_coefficient_of_variation_results_TCGA_subjects_1_coef_per_gene_normal_no_G4.txt"))
matrix_tumor_g4 <- read.table(paste0(path_tcga,"summary_coefficient_of_variation_results_TCGA_subjects_1_coef_per_gene_tumor_G4.txt"))
matrix_tumor_no_g4 <- read.table(paste0(path_tcga,"summary_coefficient_of_variation_results_TCGA_subjects_1_coef_per_gene_tumor_no_G4.txt"))

number <- 6  # from 1 to 15
boxplot(matrix_normal_g4[,number], matrix_tumor_g4[,number], names=c("Normal","Tumor"))
boxplot(matrix_tumor_g4[,number], matrix_tumor_no_g4[,number], names=c("Normal","Tumor"))
# Variation is bigger in tumor for both G4 and no G4


## Part 2 - 1 coef per individual ----
# Already calculated under "Boxplots comparing variability in G4 genes"


# Comparing pG4 vs. no pG4 for tumor and healthy ----
setwd("C:/Users/ruthb/Documents/UTH-MDA/Master's files/2020-2_and_2021-1/new tcga analysis/")

## A) Mean expression ----
# Looks like those values are already log2+1 transformed
tcga_normal_no_g4 <- read.table(paste0(path_tcga, "results_normal_nog4.txt"), header=T)
tcga_normal_g4 <- read.table(paste0(path_tcga, "results_normal_g4.txt"), header=T)
tcga_tumor_no_g4 <- read.table(paste0(path_tcga, "results_tumor_nog4.txt"), header=T)
tcga_tumor_g4 <- read.table(paste0(path_tcga, "results_tumor_g4.txt"), header=T)

# Remove inf values
# # nrow(tcga_normal_no_g4)
# # df <- tcga_normal_no_g4[is.finite(rowSums(tcga_normal_no_g4)), ]
# # nrow(df)
# tcga_normal_no_g4 <- tcga_normal_no_g4[is.finite(rowSums(tcga_normal_no_g4)), ]
# tcga_normal_g4 <- tcga_normal_g4[is.finite(rowSums(tcga_normal_g4)), ]
# tcga_tumor_no_g4 <- tcga_tumor_no_g4[is.finite(rowSums(tcga_tumor_no_g4)), ]
# tcga_tumor_g4 <- tcga_tumor_g4[is.finite(rowSums(tcga_tumor_g4)), ]

wilcox_matrix <- as.data.frame(matrix(nrow=ncol(tcga_normal_g4), ncol=6))
rownames(wilcox_matrix) <- colnames(tcga_normal_g4)
colnames(wilcox_matrix) <- c("Normal_p","Normal_mean_g4","Normal_mean_no_g4","Tumor_p","Tumor_mean_g4","Tumor_mean_no_g4")

for (i in 1:ncol(tcga_normal_g4)) {
  
  wilcox_matrix[i,"Normal_p"] <- wilcox.test(tcga_normal_g4[,i], tcga_normal_no_g4[,i])$p.value
  wilcox_matrix[i,"Normal_mean_g4"] <- mean(tcga_normal_g4[,i])
  wilcox_matrix[i,"Normal_mean_no_g4"] <- mean(tcga_normal_no_g4[,i])
  # wilcox_matrix[i,"Normal_p"] <- wilcox.test(tcga_normal_g4[tcga_normal_g4[,i] > 0,][,i], tcga_normal_no_g4[tcga_normal_no_g4[,i] > 0,][,i])$p.value
  # wilcox_matrix[i,"Normal_mean_g4"] <- mean(tcga_normal_g4[tcga_normal_g4[,i] > 0,][,i])
  # wilcox_matrix[i,"Normal_mean_no_g4"] <- mean(tcga_normal_no_g4[tcga_normal_no_g4[,i] > 0,][,i])
  
  wilcox_matrix[i,"Tumor_p"] <- wilcox.test(tcga_tumor_g4[,i], tcga_tumor_no_g4[,i])$p.value
  wilcox_matrix[i,"Tumor_mean_g4"] <- mean(tcga_tumor_g4[,i])
  wilcox_matrix[i,"Tumor_mean_no_g4"] <- mean(tcga_tumor_no_g4[,i])
  # wilcox_matrix[i,"Tumor_p"] <- wilcox.test(tcga_tumor_g4[tcga_tumor_g4[,i] > 0,][,i], tcga_tumor_no_g4[tcga_tumor_no_g4[,i] > 0,][,i])$p.value
  # wilcox_matrix[i,"Tumor_mean_g4"] <- mean(tcga_tumor_g4[tcga_tumor_g4[,i] > 0,][,i])
  # wilcox_matrix[i,"Tumor_mean_no_g4"] <- mean(tcga_tumor_no_g4[tcga_tumor_no_g4[,i] > 0,][,i])
  
}
write.table(wilcox_matrix, "wilcox_matrix_G4_vs_noG4.txt")

## *Plotting pairs of boxplots side by side ----
# Creating a long table with the data. One for normal (g4 vs non g4), one for tumor (g4 vs non g4)
library(reshape2)

########### OPTIONAL - Restricting to only a subset of genes ###################
setwd("C:/Users/ruthb/Documents/UTH-MDA/Master's files/2020-1/RAW DATA/g4_known-transc_tcga_expression")

cancers_of_interest <- c("LIHC","HNSC","ESCA","STAD","LUAD","LUSC","THCA","UCEC","BRCA","BLCA","PRAD","COAD","KICH","KIRC","KIRP")

new_tcga_avg_matrix_normal <- as.data.frame(matrix(nrow=20530,ncol=length(cancers_of_interest)+1))
colnames(new_tcga_avg_matrix_normal) <- c("entrez_id", cancers_of_interest)
new_tcga_avg_matrix_normal$entrez_id <- as.data.frame(read.table(paste0("script_normal_samples/LIHC__geneExpN.txt"), header=T))[,2]

new_tcga_avg_matrix_tumor <- as.data.frame(matrix(nrow=20530,ncol=length(cancers_of_interest)+1))
colnames(new_tcga_avg_matrix_tumor) <- c("entrez_id", cancers_of_interest)
new_tcga_avg_matrix_tumor$entrez_id <- as.data.frame(read.table(paste0("script_tumor_samples/LIHC__geneExpT.txt"), header=T))[,2]

# test1 <- as.data.frame(read.table(paste0("script_normal_samples/LIHC__geneExpN.txt"), header=T))
# test2 <- as.data.frame(read.table(paste0("script_normal_samples/HNSC__geneExpN.txt"), header=T))
# test3 <- as.data.frame(read.table(paste0("script_normal_samples/ESCA__geneExpN.txt"), header=T))

for (i in cancers_of_interest) {
  
  temp_matrix <- as.data.frame(read.table(paste0("script_normal_samples/",i,"__geneExpN.txt"), header=T))
  temp_matrix[,-c(1,2)] <- log2(temp_matrix[,-c(1,2)])
  new_tcga_avg_matrix_normal[,i] <- rowMeans(temp_matrix[,-c(1,2)])
  
  temp_matrix <- as.data.frame(read.table(paste0("script_tumor_samples/",i,"__geneExpT.txt"), header=T))
  temp_matrix[,-c(1,2)] <- log2(temp_matrix[,-c(1,2)])
  new_tcga_avg_matrix_tumor[,i] <- rowMeans(temp_matrix[,-c(1,2)])
  
}

setwd("C:/Users/ruthb/Documents/UTH-MDA/Master's files/2020-2_and_2021-1/new tcga analysis/")
write.table(new_tcga_avg_matrix_normal, "all_normal_avg_new_log2.txt", col.names=T, row.names=F)
write.table(new_tcga_avg_matrix_tumor, "all_tumor_avg_new_log2.txt", col.names=T, row.names=F)

# Load genes of interest (preferentially entrez ids)
selected_g4_genes <- read.table("C:/Users/ruthb/Documents/UTH-MDA/Cell cycle cancer genes expression/cell_cycle_g4_genes.txt")
selected_no_g4_genes <- read.table("C:/Users/ruthb/Documents/UTH-MDA/Cell cycle cancer genes expression/cell_cycle_no_g4_genes.txt")

rownames(new_tcga_avg_matrix_normal) <- new_tcga_avg_matrix_normal$entrez_id
rownames(new_tcga_avg_matrix_tumor) <- new_tcga_avg_matrix_tumor$entrez_id

tcga_normal_g4 <- new_tcga_avg_matrix_normal[new_tcga_avg_matrix_normal$entrez_id %in% selected_g4_genes$V1,][,-1]
tcga_normal_no_g4 <- new_tcga_avg_matrix_normal[new_tcga_avg_matrix_normal$entrez_id %in% selected_no_g4_genes$V1,][,-1]
tcga_tumor_g4 <- new_tcga_avg_matrix_tumor[new_tcga_avg_matrix_tumor$entrez_id %in% selected_g4_genes$V1,][,-1]
tcga_tumor_no_g4 <- new_tcga_avg_matrix_tumor[new_tcga_avg_matrix_tumor$entrez_id %in% selected_no_g4_genes$V1,][,-1]

# Replace -inf with zeroes
tcga_normal_no_g4[tcga_normal_no_g4 == -Inf] <- 0
tcga_tumor_no_g4[tcga_tumor_no_g4 == -Inf] <- 0
tcga_normal_g4[tcga_normal_g4 == -Inf] <- 0
tcga_tumor_g4[tcga_tumor_g4 == -Inf] <- 0

# Calculating the statistics
wilcox_matrix <- as.data.frame(matrix(nrow=ncol(tcga_normal_g4), ncol=6))
rownames(wilcox_matrix) <- colnames(tcga_normal_g4)
colnames(wilcox_matrix) <- c("Normal_p","Normal_mean_g4","Normal_mean_no_g4","Tumor_p","Tumor_mean_g4","Tumor_mean_no_g4")

for (i in 1:ncol(tcga_normal_g4)) {
  
  wilcox_matrix[i,"Normal_p"] <- wilcox.test(tcga_normal_g4[,i], tcga_normal_no_g4[,i])$p.value
  wilcox_matrix[i,"Normal_mean_g4"] <- mean(tcga_normal_g4[,i])
  wilcox_matrix[i,"Normal_mean_no_g4"] <- mean(tcga_normal_no_g4[,i])
  # wilcox_matrix[i,"Normal_p"] <- wilcox.test(tcga_normal_g4[tcga_normal_g4[,i] > 0,][,i], tcga_normal_no_g4[tcga_normal_no_g4[,i] > 0,][,i])$p.value
  # wilcox_matrix[i,"Normal_mean_g4"] <- mean(tcga_normal_g4[tcga_normal_g4[,i] > 0,][,i])
  # wilcox_matrix[i,"Normal_mean_no_g4"] <- mean(tcga_normal_no_g4[tcga_normal_no_g4[,i] > 0,][,i])
  
  wilcox_matrix[i,"Tumor_p"] <- wilcox.test(tcga_tumor_g4[,i], tcga_tumor_no_g4[,i])$p.value
  wilcox_matrix[i,"Tumor_mean_g4"] <- mean(tcga_tumor_g4[,i])
  wilcox_matrix[i,"Tumor_mean_no_g4"] <- mean(tcga_tumor_no_g4[,i])
  # wilcox_matrix[i,"Tumor_p"] <- wilcox.test(tcga_tumor_g4[tcga_tumor_g4[,i] > 0,][,i], tcga_tumor_no_g4[tcga_tumor_no_g4[,i] > 0,][,i])$p.value
  # wilcox_matrix[i,"Tumor_mean_g4"] <- mean(tcga_tumor_g4[tcga_tumor_g4[,i] > 0,][,i])
  # wilcox_matrix[i,"Tumor_mean_no_g4"] <- mean(tcga_tumor_no_g4[tcga_tumor_no_g4[,i] > 0,][,i])
  
}
write.table(wilcox_matrix, "wilcox_matrix_G4_vs_noG4_cell_cycle_genes_log2.txt")

################################################################################  


long_normal_1 <- melt(tcga_normal_g4)
long_normal_2 <- melt(tcga_normal_no_g4)
long_normal_1$gene_type <- "2.with_pG4"
long_normal_2$gene_type <- "1.no_pG4"
long_normal <- rbind(long_normal_1, long_normal_2)
# long_normal <- long_normal[long_normal$value > 0,]  # removing genes with expression=0
colnames(long_normal) <- c("cancer_type","gene_expr","gene_type")

long_tumor_1 <- melt(tcga_tumor_g4)
long_tumor_2 <- melt(tcga_tumor_no_g4)
long_tumor_1$gene_type <- "2.with_pG4"
long_tumor_2$gene_type <- "1.no_pG4"
long_tumor <- rbind(long_tumor_1, long_tumor_2)
# long_tumor <- long_tumor[long_tumor$value > 0,]  # removing genes with expression=0
colnames(long_tumor) <- c("cancer_type","gene_expr","gene_type")

# # Sample data (replace this with your own data)
# set.seed(42)
# library(ggplot2)
# 
# data <- data.frame(
#   tumor_type = rep(paste("Tumor", 1:15), each = 2),
#   data_type = rep(c("Normal", "Tumor"), times = 15),
#   value = c(rnorm(30, mean = 10, sd = 2),
#             rnorm(30, mean = 15, sd = 2))
# )

# Create boxplots side by side horizontally
my_palette_3 <- c("#33BFC4","#F2756D")#"#F2756D","#33BFC4")
cancer_order <- c("LIHC","LUAD","LUSC","UCEC","KIRC","ESCA","BLCA","BRCA",
                  "KIRP","THCA","STAD","COAD","HNSC","PRAD","KICH")
cancer_order <- rev(cancer_order)
long_normal$cancer_type <- factor(long_normal$cancer_type, levels = cancer_order)
long_tumor$cancer_type <- factor(long_tumor$cancer_type, levels = cancer_order)

ggplot(long_normal, aes(x = cancer_type, y = gene_expr, fill = gene_type)) +
  geom_boxplot(position = "dodge") +
  labs(x = "", y = "Mean expression log2(fpkm + 1)", title = "") +
  # scale_y_continuous(breaks = seq(0, 25, by = 5), limits = c(0, 25)) +
  theme_minimal() +
  coord_flip() +
  scale_fill_manual(values = my_palette_3) +
  theme(axis.text.x = element_text(size = 10.5, color = "black")) +
  theme(axis.text.y = element_text(size = 10.5, color = "black"))
# 8.35 x 4.3 (landscape) or 4.5 x 4.5 (portrait)


## B) Coef. of variation ----
tcga_normal_g4 <- read.table(paste0(path_tcga,"summary_coefficient_of_variation_results_TCGA_subjects_1_coef_per_gene_normal_G4.txt"))
tcga_normal_no_g4 <- read.table(paste0(path_tcga,"summary_coefficient_of_variation_results_TCGA_subjects_1_coef_per_gene_normal_no_G4.txt"))
tcga_tumor_g4 <- read.table(paste0(path_tcga,"summary_coefficient_of_variation_results_TCGA_subjects_1_coef_per_gene_tumor_G4.txt"))
tcga_tumor_no_g4 <- read.table(paste0(path_tcga,"summary_coefficient_of_variation_results_TCGA_subjects_1_coef_per_gene_tumor_no_G4.txt"))

## *Plotting pairs of boxplots side by side ----
# Creating a long table with the data. One for normal (g4 vs non g4), one for tumor (g4 vs non g4)
library(reshape2)

long_normal_1 <- melt(tcga_normal_g4)
long_normal_2 <- melt(tcga_normal_no_g4)
long_normal_1$gene_type <- "2.with_pG4"
long_normal_2$gene_type <- "1.no_pG4"
long_normal <- rbind(long_normal_1, long_normal_2)
# long_normal <- long_normal[long_normal$value > 0,]  # removing genes with expression=0
colnames(long_normal) <- c("cancer_type","coef_var","gene_type")

long_tumor_1 <- melt(tcga_tumor_g4)
long_tumor_2 <- melt(tcga_tumor_no_g4)
long_tumor_1$gene_type <- "2.with_pG4"
long_tumor_2$gene_type <- "1.no_pG4"
long_tumor <- rbind(long_tumor_1, long_tumor_2)
# long_tumor <- long_tumor[long_tumor$value > 0,]  # removing genes with expression=0
colnames(long_tumor) <- c("cancer_type","coef_var","gene_type")

long_g4_1 <- melt(tcga_normal_g4)
long_g4_2 <- melt(tcga_tumor_g4)
long_g4_1$gene_type <- "2.normal"
long_g4_2$gene_type <- "1.tumor"
long_g4 <- rbind(long_g4_1, long_g4_2)
colnames(long_g4) <- c("cancer_type","coef_var","gene_type")

long_noG4_1 <- melt(tcga_normal_no_g4)
long_noG4_2 <- melt(tcga_tumor_no_g4)
long_noG4_1$gene_type <- "2.normal"
long_noG4_2$gene_type <- "1.tumor"
long_noG4 <- rbind(long_noG4_1, long_noG4_2)
colnames(long_noG4) <- c("cancer_type","coef_var","gene_type")

# Create boxplots side by side horizontally
my_palette_3 <- c("#33BFC4","#F2756D")#"#F2756D","#33BFC4")
cancer_order <- c("LUAD","COAD","LUSC","KIRC","KIRP","STAD","PRAD","HNSC",
                  "KICH","BRCA","BLCA","LIHC","THCA","UCEC","ESCA")
cancer_order <- rev(cancer_order)
long_normal$cancer_type <- factor(long_normal$cancer_type, levels = cancer_order)
long_tumor$cancer_type <- factor(long_tumor$cancer_type, levels = cancer_order)
long_g4$cancer_type <- factor(long_g4$cancer_type, levels = cancer_order)
long_noG4$cancer_type <- factor(long_noG4$cancer_type, levels = cancer_order)

df <- long_noG4
# Sanity check
# boxplot(df[df$cancer_type %in% "STAD" & df$gene_type %in% "2.normal",]$coef_var,
#         df[df$cancer_type %in% "STAD" & df$gene_type %in% "1.tumor",]$coef_var)

ggplot(df[!is.na(df$coef_var),], aes(x = cancer_type, y = coef_var, fill = gene_type)) +
  #coord_cartesian(ylim=c(0, 300)) +  # zoom without removing datapoints from calculation
  geom_boxplot(position = "dodge") +
  labs(x = "", y = "Coef. of variation of expression across individuals", title = "") +
  #scale_y_continuous(breaks = seq(0, 300, by = 100), limits = c(0, 300)) +  # removes data points from mean calculation
  theme_minimal() +
  coord_flip() +
  scale_fill_manual(values = my_palette_3) +
  theme(axis.text.x = element_text(size = 10.5, color = "black")) +
  theme(axis.text.y = element_text(size = 10.5, color = "black"))
# 8.35 x 4.3 (landscape) or 4.5 x 4.5 (portrait)
# 4.5 x 7 (landscape) when there are too many outliers



source('csg_heatmap.R')
library(heatmaply)

# Load the CSV data file located at the specified directory

df_data_org <- read.csv("/Users/saurabhbhardwaj/Desktop/CONTI/sample_data.csv", header = TRUE)

# The 'select_top_genes' function provides the top quality features as per the chosen reference value
# Inputs:
# df_data_org: A dataframe where the rows are the features and columns are the samples. The first column is the name of the features.
# k: A list of samples in each group.
# ref: The reference pattern. The length of this vector should be equal to the number of groups in the data.
# num_genes: The number of top quality genes to be selected based on cosine similarity.

# The function returns a list containing two elements:
# 1. A dataframe with the top 'num_genes' genes sorted by cosine similarity to the reference pattern and their respective group-wise expression values.
# 2. A list containing the indices of the samples in each group.

cos_matrix_sort_label_u <- select_top_genes(df_data_org, k = list(4,3,3), ref = c(1, 0, 0), num_genes = 12)

# # path for reordering of samples
metadata_path_arranged <- read.csv("/Users/saurabhbhardwaj/Desktop/CONTI/sample_order.csv", header = TRUE)
# Assign cosine matrix sorted by label to tot_sg variable
tot_sg <- cos_matrix_sort_label_u
data_tmm_path_arranged <- df_data_org[,colnames(metadata_path_arranged)]
idx <- match(cos_matrix_sort_label_u[[1]]$Genes, data_tmm_path_arranged$X)
tot_sg[[1]][,2:ncol(df_data_org)] <- data_tmm_path_arranged[idx,2:ncol(df_data_org)]
colName <- colnames(data_tmm_path_arranged)[2:ncol(df_data_org)]
geneName <- tot_sg[[1]][,1]
tot_sg1 <- tot_sg[[1]][,2:ncol(tot_sg[[1]])]

tot_sg_scale_log_center <- heatmap_data_modification(tot_sg1)


# Plot Heatmap
data_plot <- tot_sg_scale_log_center
thld <- 2.5
data_plot[data_plot > thld] <- thld
data_plot[data_plot < -thld] <- -thld

c <- redblue(300)
#install.packages("heatmaply")


heatmaply((data_plot),colors = rgb(c), cluster = FALSE, dendrogram = c("none") )


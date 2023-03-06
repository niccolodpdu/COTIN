
#' Function to calculate group means
#' ---------------------------------------------
#' This function "calculate_mean" takes two arguments, a matrix X and a vector k. 
#' The function calculates the mean of each group defined by the values in 
#' vector k for each row of matrix X.
#'
#' @param X a matrix
#' @param k a vector
#' @return a matrix of group means
#' @examples
#' calculate_mean(matrix(1:6, nrow = 2), c(2, 1))
#' ---------------------------------------------

calculate_mean <- function(X, k){
            nrow <- nrow(X)
            ncol <- length(k)
            gr_mean <- matrix(nrow = nrow, ncol = ncol)
            cumsum_k <- cumsum(k)
            for(i in 1:nrow) {
                        for(j in 1:ncol){
                                    if(j==1){
                                                gr_mean[i,j] <- mean(X[i, 1:k[[j]]])
                                    }else{
                                                gr_mean[i,j] <- mean(X[i, (cumsum_k[j-1]+1):cumsum_k[j]])
                                    }
                        }
            }
            return(gr_mean)
}

#' Function to calculate the cosine similarity
#'
#' This code defines a function called "getCosineSimilarity" 
#' that takes in two vectors, x and y, and returns their cosine similarity. 
#'
#' @param x a numeric vector
#' @param y a numeric vector
#' @return a numeric value between -1 and 1
#' @examples
#' getCosineSimilarity(c(1, 2, 3), c(3, 2, 1))

getCosineSimilarity <- function(x,y){
            if (!is.vector(x) || !is.vector(y)) {
                        stop("x and y have to be vectors!")
            }
            if (length(x) != length(y)) {
                        stop("x and y have to be the same length!")
            }
            xy <- sum(x * y)
            nx <- sqrt(sum(x^2))
            ny <- sqrt(sum(y^2))
            nxny <- nx*ny
            Cs <- xy/nxny
            return(Cs)
}

#' To create indices that can be used to sample data points from a dataset. 
#'
#' The function takes a vector k as an input which specifies the number of 
#' data points to be sampled from each group.
#'
#' @param k a vector
#' @return a list of indices for each group
#' @examples
#' create_indices(c(2, 1))

create_indices <- function(k) {
            idx_Sample <- list()
            start <- 1
            for (i in 1:length(k)) {
                        end <- start + k[[i]] - 1
                        idx_Sample[[i]] <- c(start:end)
                        start <- end + 1
            }
            return(idx_Sample)
}

#' To generate a colormap with shades of red and blue colors. 
#'
#' The colormap is used to create a gradient from the blue color to the red color. 
#' The function takes an optional argument m which specifies the number of colors in the colormap.
#'
#' @param m an integer, number of colors in the colormap
#' @return a matrix with three columns (r, g, b) representing the color map
#' @examples
#' redblue()
#' 
redblue <- function(m = nrow(colormap())) {
            if (m %% 2 == 0) {
                        # From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
                        m1 <- m * 0.5
                        r <- (0:(m1 - 1))/max(m1 - 1, 1)
                        g <- r
                        r <- c(r, rep(1, m1))
                        g <- c(g, rev(g))
                        b <- rev(r)
            } else {
                        # From [0 0 1] to [1 1 1] to [1 0 0];
                        m1 <- floor(m * 0.5)
                        r <- (0:(m1 - 1))/max(m1, 1)
                        g <- r
                        r <- c(r, rep(1, m1 + 1))
                        g <- c(g, 1, rev(g))
                        b <- rev(r)
            }
            cbind(r, g, b)
}

#' Select top genes based on cosine similarity
#'
#'
#' @param df_data_org A data.frame with the gene expression data. 
#' @param k The number of groups to be defined for group-wise mean expression calculation.
#' @param ref A numeric vector representing the reference pattern against which the cosine similarity is calculated.
#' @param num_genes The number of genes to be selected based on cosine similarity score.
#' @return A list containing a dataframe with the selected genes and their group-wise expression values sorted by cosine similarity and a list of indices of the samples in each group.
#' @examples
#' select_top_genes(df_data_org = df, k = 2, ref = c(1, 0), num_genes = 5)


# The function takes the following inputs:
# df_data_org: A dataframe where the rows are the features and columns are the samples. The first column is the name of the features.
# k: A list of samples in each group.
# ref: The reference pattern. The length of this vector should be equal to the number of groups in the data.
# num_genes: The number of top quality genes to be selected based on cosine similarity.

# The function returns a list containing two elements:
# 1. A dataframe with the top num_genes genes sorted by cosine similarity to the reference pattern and their respective group-wise expression values.
# 2. A list containing the indices of the samples in each group.

select_top_genes <- function(df_data_org, k, ref, num_genes) {
            
            # Extract the expression data from the input dataframe.
            X <- df_data_org[,2:ncol(df_data_org)]
            X <- as.matrix(X)
            
            # Get the number of rows and columns in the expression data.
            nrow <- nrow(X)
            ncol <- ncol(X)
            
            # Get the total number of samples in the data.
            num_samples <- ncol(df_data_org)-1   
            
            # Calculate the group-wise mean expression for each feature.
            gr_mean <- calculate_mean(X, k)
            gr_mean_df <- as.data.frame(gr_mean)
            
            # Add the feature names as a column to the group-wise mean expression dataframe.
            gr_mean_df <- cbind(df_data_org[,1],gr_mean_df)
            colnames(gr_mean_df) <- c("Genes","col1","col2","col3")
            
            # Transpose the group-wise mean expression data for cosine similarity calculation.
            X <- gr_mean_df[,2:ncol(gr_mean_df)]
            X <- t(X)
            
            # Calculate the cosine similarity between each sample and the reference pattern.
            cos_sim_u <- rep(0,ncol(X))
            for (i in 1:ncol(X)) {
                        cos_sim_u[i] <- getCosineSimilarity(X[,i],ref)
            }
            cos_sim_u <- as.matrix(cos_sim_u)
            cos_sim_u <- as.data.frame(cos_sim_u)
            
            # Combine the cosine similarity values with the gene names in a dataframe.
            cos_matrix_u <- data.frame(Genes=gr_mean_df[,1], cos_sim_u)
            
            # Sort the dataframe by cosine similarity and select the top num_genes genes.
            cos_matrix_sort_u <- cos_matrix_u[order(-cos_matrix_u$V1),]
            cos_matrix_sort_u <- cos_matrix_sort_u[1:num_genes, ]
            
            # Get the indices of the selected genes in the original dataframe.
            idx <- match(cos_matrix_sort_u$Genes, gr_mean_df$Genes)
            
            # Combine the sorted cosine similarity dataframe with the original dataframe to get the group-wise expression values for the selected genes.
            cos_matrix_sort_label_u <- cbind(cos_matrix_sort_u, gr_mean_df[idx,2:ncol(gr_mean_df)])
            
            # Get the indices of the samples in each group.
            idx_Sample <- create_indices(k)
            
            # Return a list containing the sorted cosine similarity dataframe and the list of sample indices in each group.
            return(list(cos_matrix_sort_label_u, idx_Sample))
            
}


#' Modify Data for Heatmap Visualization
#'
#' This function modifies the input data for visualization in a heatmap by performing scaling, log transformation, and centering.
#'
#' @param tot_sg1 A numeric matrix containing the input data to be modified.
#'
#' @return A modified numeric matrix containing the scaled, log-transformed, and centered data.
#'
#' @examples
#' tot_sg1 <- matrix(runif(100), nrow = 10)
#' modified_data <- heatmap_data_modification(tot_sg1)
#'
# Function to modify data as per new heatmap design
# Inputs:
#   tot_sg1 - data frame containing numeric data
# Outputs:
#   tot_sg_scale_log_center - modified data frame

heatmap_data_modification <- function(tot_sg1) {
            # Calculate row sums and factor for scaling
            tot_sg_sum <- rowSums(tot_sg1)
            factor <- 1000000 / tot_sg_sum
            
            # Scale data and set values less than 1 to 1
            tot_sg_scale <- tot_sg1 * factor
            tot_sg_scale[tot_sg_scale < 1] <- 1
            
            # Log transform scaled data and calculate mean and standard deviation
            tot_sg_scale_log <- log(tot_sg_scale)
            mean_gb <- mean(as.matrix(tot_sg_scale_log))
            std_gb <- sd(as.matrix(tot_sg_scale_log, na.rm = TRUE))
            
            # Center data around mean and standard deviation
            tot_sg_scale_log_center <- (tot_sg_scale_log - mean_gb) / std_gb
            
            # Return modified data
            return(tot_sg_scale_log_center)
}



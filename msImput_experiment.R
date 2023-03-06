library(Ternary)
library(hrbrthemes)
library(ggplot2)
library(DirichletReg)
library(DescTools)
library(psych)
library(scales)
library(devEMF)
library(rgl)
library(magick)
library(ggtern)
library(plyr)
source('Cosbin_functions.R')
source('msImput.R')
#source('Cosbin_functions_visualization.R')
#source('Cosbin_functions_evaluation.R')


### Description ###

# 1. Simulation data, then virtually mask some features out if their missing rates are over than 
#    a given threshold in all groups
# 2. Pre-imputation, and compare the SGs that are imputed here but previously masked-out in step 1 
#    (quantitatively, by cosine value: deviation-in-degree of SGs)
# 3. Scatter simplex
# 4. Cosbin v.s. Totalcount again: by accessing the deviation-in-degree over iCEGs in normalized data



###################################### Generate ideal simulation data ###################################################

# 30000 values (30 samples x 1000 genes) in total

###### Simulation setting

MAR=60000
Total_missing=70000

SG_threshold<-0.98
iCEG_threshold<-0.95
missing_rate_threshold<-0.6
MAR_option<-'0' # MAR happens in SG only ('SG only') or not (any other value)

dirichlet_param<-2

nGroup <- 5
nGene <- 1000
psDEG <-  c(0.06,0.09,0.12,0.15,0.18) # absolute percentage over all genes
sum(psDEG) # check the proportion of SG here

low_expression_SG<-0.1
high_expression_SG<-1

sd_scaler=3 # for COTIN imputation

###### Starting generating simulation data

ref <- vector("list", length = nGroup)
for (i in 1:nGroup) {
  ref[[i]] <- rep(low_expression_SG, nGroup)
  ref[[i]][i] <- high_expression_SG
}

non_DEG <- rdirichlet(nGene-sum(psDEG)*nGene, rep(dirichlet_param, nGroup))

### SG
data_SG <- NULL
for (i in 1:nGroup) {
  data_SG <-
    rbind(data_SG, matrix(rep(ref[[i]], nGene * psDEG[i]), ncol = nGroup, byrow = TRUE))
}

### sDEG (iDEG + noise)
noise <- rnorm(nGroup * nGene * sum(psDEG), mean = 0, sd = 0.04)
data_SG <- abs(data_SG + noise)

data <- rbind(data_SG, non_DEG)


### Add replicates

nRep <- rep(30, nGroup) 

data_grnd_full <- NULL
for (i in 1:nGroup) {
  data_grnd_full <-
    cbind(data_grnd_full, matrix(rep(data[, i], nRep[i]), ncol = nRep[i], byrow = FALSE))
}

dim(data_grnd_full)
noise <- rnorm(nRep * nGroup * nGene, mean = 0, sd = 0.01)
data_grnd_full <- abs(data_grnd_full + noise)

for(i in 1:nGene){
  data_grnd_full[i,] <- data_grnd_full[i,] * 1000
}

### Ground truth data (super sample, to be compared with the normalized result)
data_grnd_super_sample <- NULL
start <- 1
for(rep in nRep){
  data_grnd_super_sample <- cbind(data_grnd_super_sample, rowMeans(data_grnd_full[, start : (start + rep - 1)]))
  start <- start + rep
}

### Add label
ind_CEG <- which(cos_iCEG(data_grnd_super_sample) >= iCEG_threshold)
ind_sDEG <- seq(1,sum(psDEG)*nGene,1)
ind_neither <- setdiff(seq(1, nGene, 1), c(ind_CEG, ind_sDEG))

ref <- vector("list", length = nGroup)
for (i in 1:nGroup) {
  ref[[i]] <- rep(0, nGroup)
  ref[[i]][i] <- 1
}


ind_sDEG2 <- list()
for (i in 1:nGroup) {
  ind_sDEG2[[i]] <- which(cos_ref(data_grnd_super_sample, ref[[i]]) >= SG_threshold)
}

lbl <- rep(0, nGene)
lbl[ind_CEG] <- 'iCEG (labeled)'
lbl[ind_neither] <- 'DEG'
lbl[ind_sDEG] <- 'sDEG'

number_iCEG_grnd<-length(ind_CEG)

lbl2 <- lbl
lbl3 <- lbl
for (i in 1:nGroup) {
  lbl2[ind_sDEG2[[i]]] <- paste0('sDEG G', i)
  lbl3[ind_sDEG2[[i]]] <- i
}

lbl2[lbl2=="sDEG"]<-'DEG'

### Ground truth label
CEG_grnd <- rep(0, nGene)
CEG_grnd[ind_CEG] <- 1
message("Check the number of iCEG here")
sum(CEG_grnd)

sDEG_grnd <- rep(0, nGene)
sDEG_grnd[ind_sDEG] <- 1

DEG_grnd <- rep(1, nGene)
DEG_grnd <- DEG_grnd-(sDEG_grnd + CEG_grnd)
sum(DEG_grnd)

DEG_and_sDEG_grnd <- rep(1, nGene)
DEG_and_sDEG_grnd <- DEG_and_sDEG_grnd-CEG_grnd
sum(DEG_and_sDEG_grnd)


#################### Plot 1. Scatter simplex plot of the original data matrix (without missing values) ############

# Over 3 groups:
X <- data_grnd_super_sample
Xproj <- X
Xproj <- apply(Xproj, 1, function(x) x / sum(x))

A <- diag(dim(Xproj)[1])
K = dim(A)[2]
PS <- t(matrix(c(cos((seq_len(K) - 1) * 2 * pi / K),
                 sin((seq_len(K) - 1) * 2 * pi / K)), K))
library(corpcor)
mg <-list()
for (i in 1:nGroup){
  mg[[i]]<-which(lbl3==i)
}

mg_grnd<-mg

tmp <- PS %*% pseudoinverse(A)
Xp <- tmp %*% Xproj

# color-coding for SGs
plot(t(PS), xlab = 'Ground Truth', ylab = NA, asp = 1, pch = '.')
points(Xp[1,], Xp[2,], col = rgb(80, 80, 80, max = 255, alpha = 125),
       xlab = NA, ylab = NA, asp = 1)

mg.col = c('dark green','orange','yellow','dodgerblue','blue')
for(i in seq_along(mg)){
  points(Xp[1,mg[[i]]], Xp[2,mg[[i]]],
         col=mg.col[i])
}
points(t(PS), pch = 17)


######################### Introducing missing values ############################

### MAR
removed<-data_grnd_full

if (MAR_option == 'SG only'){
  index<-sample(which(data_grnd_full[1:(sum(psDEG)*nGene),]>0),MAR)
} else {
  index<-sample(which(data_grnd_full>0),MAR)
}

for (i in 1:MAR){
  removed[index]<-0
}
sum(removed == 0)  

### LLOD
removed[removed < quantile(removed,(Total_missing/(dim(data_grnd_full)[1]*dim(data_grnd_full)[2])))]<-0
sum(removed == 0)  


input_data <- removed
before_norm<-input_data

before_remove<-input_data
after_remove<-NULL

## Mask-out
missing_rate<-apply(after_norm,1,function(x)
  sum(x==0)/length(x))

after_remove<-after_norm[missing_rate<missing_rate_threshold,]


#################### Plot 2. Red circle on those got masked-out #####################
X <- data_grnd_super_sample
Xproj <- X
Xproj <- apply(Xproj, 1, function(x) x / sum(x))

A <- diag(dim(Xproj)[1])
K = dim(A)[2]
PS <- t(matrix(c(cos((seq_len(K) - 1) * 2 * pi / K),
                 sin((seq_len(K) - 1) * 2 * pi / K)), K))
library(corpcor)
mg <-list()
for (i in 1:nGroup){
  mg[[i]]<-which(lbl3==i)
}

mg_grnd<-mg

tmp <- PS %*% pseudoinverse(A)
Xp <- tmp %*% Xproj

# color-coding for SGs
plot(t(PS), xlab = 'Ground Truth Masked-out', ylab = NA, asp = 1, pch = '.')
points(Xp[1,], Xp[2,], col = rgb(80, 80, 80, max = 255, alpha = 125),
       xlab = NA, ylab = NA, asp = 1)


index_masked_out<-which(missing_rate>=missing_rate_threshold)

for(i in 1:length(index_masked_out)){
  points(Xp[1,index_masked_out[i]], Xp[2,index_masked_out[i]],
         col='hotpink')
}
points(t(PS), pch = 17)



mg_grnd_masked_out<-intersect(index_masked_out,unlist(mg_grnd))
non_SG_masked_out<-setdiff(index_masked_out,mg_grnd_masked_out)

for(i in 1:length(mg_grnd_masked_out)){
  points(Xp[1,mg_grnd_masked_out[i]], Xp[2,mg_grnd_masked_out[i]],
         col='red')
}
points(t(PS), pch = 17)


legend('bottomright',"Legend",c("Masked-out non-SG","Masked-out SG"),
       fill=c('hotpink','red'),cex=0.42,pt.cex = 1)



################# The msImput, The Mean and The Min/2: Their performance in imputing masked-out genes ###################

# COTIN MsImput


COTIN_imputed<-msImput(after_norm,nRep,'global',missing_rate_threshold=1, sd_scaler=sd_scaler)


# See which mgs got masked-out

COTIN_imputed_masked_out<-COTIN_imputed[mg_grnd_masked_out,]
after_missing_masked_out<-after_norm[mg_grnd_masked_out,]
before_missing_masked_out<-data_grnd_full[mg_grnd_masked_out,]

before_missing_super_sample <- NULL
start <- 1
for(rep in nRep){
  before_missing_super_sample <- cbind(before_missing_super_sample, 
                                       rowMeans(before_missing_masked_out[, start : (start + rep - 1)])) 
  start <- start + rep
}



## Functions for evaluation


rmse <- function(ground_truth,imputed) sqrt(mean((ground_truth - imputed)^2))
nrmse <- function(ground_truth,imputed) rmse(ground_truth,imputed)/sd(ground_truth)
euclidean <- function(a, b) sqrt(sum((a - b)^2))



####  Min/2 imputed  ####
Half_min<-min(after_norm[after_norm!=0])/2
Half_min_imputed<-after_norm
Half_min_imputed[Half_min_imputed==0]<-Half_min

# SG
SG_Half_min_error<-euclidean(Half_min_imputed[mg_grnd_masked_out,],before_missing_masked_out)
Averaged_SG_Half_min_error<-SG_Half_min_error/length(mg_grnd_masked_out)

Half_min_imputed_super_sample <- NULL
start <- 1
for(rep in nRep){
  Half_min_imputed_super_sample <- cbind(Half_min_imputed_super_sample, 
                                                 rowMeans(Half_min_imputed[, start : (start + rep - 1)])) 
  start <- start + rep
}

Half_min_imputed_super_sample_masked_out<-Half_min_imputed_super_sample[mg_grnd_masked_out,]

Deviation_in_degree_SG<-NULL
for (i in 1:dim(Half_min_imputed_super_sample_masked_out)[1]){
  temp_cosine<-mycosine(Half_min_imputed_super_sample_masked_out[i,],before_missing_super_sample[i,])
  Deviation_in_degree_SG<-c(Deviation_in_degree_SG,(180 * acos(temp_cosine)) / pi)
}
Averaged_Half_min_Deviation_in_degree_SG<-sum(Deviation_in_degree_SG)/length(mg_grnd_masked_out)


# non SG
Averaged_non_SG_Half_min_error<-euclidean(Half_min_imputed[non_SG_masked_out,],data_grnd_full[non_SG_masked_out,])/
  length(non_SG_masked_out)

Half_min_imputed_super_sample_non_SG_masked_out<-Half_min_imputed_super_sample[non_SG_masked_out,]

Deviation_in_degree_non_SG<-NULL
for (i in 1:length(non_SG_masked_out)){
  temp_cosine<-mycosine(Half_min_imputed_super_sample_non_SG_masked_out[i,],data_grnd_super_sample[non_SG_masked_out,][i,])
  Deviation_in_degree_non_SG<-c(Deviation_in_degree_non_SG,(180 * acos(temp_cosine)) / pi)
}

Averaged_Half_min_Deviation_in_degree_non_SG<-sum(Deviation_in_degree_non_SG)/length(non_SG_masked_out)



# Together
Half_min_together_error <- euclidean(data_grnd_full[c(mg_grnd_masked_out,non_SG_masked_out),],
                                         Half_min_imputed[c(mg_grnd_masked_out,non_SG_masked_out),])
Averaged_Half_min_together_error <- Half_min_together_error/length(c(mg_grnd_masked_out,non_SG_masked_out))

Averaged_Half_min_Deviation_in_degree<-(sum(Deviation_in_degree_non_SG)+
                                      sum(Deviation_in_degree_SG))/length(c(mg_grnd_masked_out,non_SG_masked_out))


# rmse & nrmse
rmse_Half_min_sg <- rmse(before_missing_masked_out,Half_min_imputed[mg_grnd_masked_out,])
rmse_Half_min_non_sg <- rmse(data_grnd_full[non_SG_masked_out,],Half_min_imputed[non_SG_masked_out,])
rmse_Half_min_all <- rmse(data_grnd_full[c(mg_grnd_masked_out,non_SG_masked_out),],
                      Half_min_imputed[c(mg_grnd_masked_out,non_SG_masked_out),])

nrmse_Half_min_sg <- nrmse(before_missing_masked_out,Half_min_imputed[mg_grnd_masked_out,])
nrmse_Half_min_non_sg <- nrmse(data_grnd_full[non_SG_masked_out,],Half_min_imputed[non_SG_masked_out,])
nrmse_Half_min_all <- nrmse(data_grnd_full[c(mg_grnd_masked_out,non_SG_masked_out),],
                        Half_min_imputed[c(mg_grnd_masked_out,non_SG_masked_out),])


####  Mean imputed  ####
Mean<-mean(after_norm[after_norm!=0])
Mean_imputed<-after_norm
Mean_imputed[Mean_imputed==0]<-Mean

# SG
SG_Mean_error<-euclidean(Mean_imputed[mg_grnd_masked_out,],before_missing_masked_out)
Averaged_SG_Mean_error<-SG_Mean_error/length(mg_grnd_masked_out)

Mean_imputed_super_sample <- NULL
start <- 1
for(rep in nRep){
  Mean_imputed_super_sample <- cbind(Mean_imputed_super_sample, 
                                     rowMeans(Mean_imputed[, start : (start + rep - 1)])) 
  start <- start + rep
}

Mean_imputed_super_sample_masked_out<-Mean_imputed_super_sample[mg_grnd_masked_out,]

Deviation_in_degree_SG<-NULL
for (i in 1:dim(Mean_imputed_super_sample_masked_out)[1]){
  temp_cosine<-mycosine(Mean_imputed_super_sample_masked_out[i,],before_missing_super_sample[i,])
  Deviation_in_degree_SG<-c(Deviation_in_degree_SG,(180 * acos(temp_cosine)) / pi)
}

Mean_Deviation_in_degree_SG<-sum(Deviation_in_degree_SG)
Averaged_Mean_Deviation_in_degree_SG<-Mean_Deviation_in_degree_SG/length(mg_grnd_masked_out)


# non-SG
Averaged_non_SG_Mean_error<-euclidean(Mean_imputed[non_SG_masked_out,],data_grnd_full[non_SG_masked_out,])/
  length(non_SG_masked_out)

Mean_imputed_super_sample_non_SG_masked_out<-Mean_imputed_super_sample[non_SG_masked_out,]

Deviation_in_degree_non_SG<-NULL
for (i in 1:length(non_SG_masked_out)){
  temp_cosine<-mycosine(Mean_imputed_super_sample_non_SG_masked_out[i,],data_grnd_super_sample[non_SG_masked_out,][i,])
  Deviation_in_degree_non_SG<-c(Deviation_in_degree_non_SG,(180 * acos(temp_cosine)) / pi)
}

Averaged_Mean_Deviation_in_degree_non_SG<-sum(Deviation_in_degree_non_SG)/length(non_SG_masked_out)

# Together
Mean_together_error <- euclidean(data_grnd_full[c(mg_grnd_masked_out,non_SG_masked_out),],
                                     Mean_imputed[c(mg_grnd_masked_out,non_SG_masked_out),])
Averaged_Mean_together_error <- Mean_together_error/length(c(mg_grnd_masked_out,non_SG_masked_out))

Averaged_Mean_Deviation_in_degree<-(sum(Deviation_in_degree_non_SG)+
                                      sum(Deviation_in_degree_SG))/length(c(mg_grnd_masked_out,non_SG_masked_out))

# rmse & nrmse
rmse_Mean_sg <- rmse(before_missing_masked_out,Mean_imputed[mg_grnd_masked_out,])
rmse_Mean_non_sg <- rmse(data_grnd_full[non_SG_masked_out,],Mean_imputed[non_SG_masked_out,])
rmse_Mean_all <- rmse(data_grnd_full[c(mg_grnd_masked_out,non_SG_masked_out),],
                      Mean_imputed[c(mg_grnd_masked_out,non_SG_masked_out),])

nrmse_Mean_sg <- nrmse(before_missing_masked_out,Mean_imputed[mg_grnd_masked_out,])
nrmse_Mean_non_sg <- nrmse(data_grnd_full[non_SG_masked_out,],Mean_imputed[non_SG_masked_out,])
nrmse_Mean_all <- nrmse(data_grnd_full[c(mg_grnd_masked_out,non_SG_masked_out),],
                        Mean_imputed[c(mg_grnd_masked_out,non_SG_masked_out),])



####  COTIN imputed  ####
SG_CONTI_error<-euclidean(COTIN_imputed_masked_out,before_missing_masked_out)
Averaged_SG_CONTI_error<-SG_CONTI_error/length(mg_grnd_masked_out)

COTIN_super_sample <- NULL
start <- 1
for(rep in nRep){
  COTIN_super_sample <- cbind(COTIN_super_sample, rowMeans(COTIN_imputed[, start : (start + rep - 1)])) 
  start <- start + rep
}

COTIN_super_sample_masked_out<-COTIN_super_sample[mg_grnd_masked_out,]
  
Deviation_in_degree_SG<-NULL
for (i in 1:dim(COTIN_super_sample_masked_out)[1]){
  temp_cosine<-mycosine(COTIN_super_sample_masked_out[i,],before_missing_super_sample[i,])
  Deviation_in_degree_SG<-c(Deviation_in_degree_SG,(180 * acos(temp_cosine)) / pi)
}

COTIN_Deviation_in_degree_SG<-sum(Deviation_in_degree_SG)
Averaged_COTIN_Deviation_in_degree_SG<-COTIN_Deviation_in_degree_SG/length(mg_grnd_masked_out)


# non-SG
Averaged_non_SG_COTIN_error<-euclidean(COTIN_imputed[non_SG_masked_out,],data_grnd_full[non_SG_masked_out,])/
  length(non_SG_masked_out)

COTIN_super_sample_non_SG_masked_out<-COTIN_super_sample[non_SG_masked_out,]

Deviation_in_degree_non_SG<-NULL
for (i in 1:length(non_SG_masked_out)){
  temp_cosine<-mycosine(COTIN_super_sample_non_SG_masked_out[i,],data_grnd_super_sample[non_SG_masked_out,][i,])
  Deviation_in_degree_non_SG<-c(Deviation_in_degree_non_SG,(180 * acos(temp_cosine)) / pi)
}

Averaged_COTIN_Deviation_in_degree_non_SG<-sum(Deviation_in_degree_non_SG)/length(non_SG_masked_out)

# Together
COTIN_together_error <- euclidean(data_grnd_full[c(mg_grnd_masked_out,non_SG_masked_out),],
                                      COTIN_imputed[c(mg_grnd_masked_out,non_SG_masked_out),])
Averaged_COTIN_together_error <- COTIN_together_error/length(c(mg_grnd_masked_out,non_SG_masked_out))

Averaged_COTIN_Deviation_in_degree<-(sum(Deviation_in_degree_non_SG)+
                                       sum(Deviation_in_degree_SG))/length(c(mg_grnd_masked_out,non_SG_masked_out))

# rmse & nrmse
rmse_CONTI_sg <- rmse(before_missing_masked_out,COTIN_imputed[mg_grnd_masked_out,])
rmse_CONTI_non_sg <- rmse(data_grnd_full[non_SG_masked_out,],COTIN_imputed[non_SG_masked_out,])
rmse_CONTI_all <- rmse(data_grnd_full[c(mg_grnd_masked_out,non_SG_masked_out),],
                       COTIN_imputed[c(mg_grnd_masked_out,non_SG_masked_out),])

nrmse_CONTI_sg <- nrmse(before_missing_masked_out,COTIN_imputed[mg_grnd_masked_out,])
nrmse_CONTI_non_sg <- nrmse(data_grnd_full[non_SG_masked_out,],COTIN_imputed[non_SG_masked_out,])
nrmse_CONTI_all <- nrmse(data_grnd_full[c(mg_grnd_masked_out,non_SG_masked_out),],
                       COTIN_imputed[c(mg_grnd_masked_out,non_SG_masked_out),])


################################ Cosbin v.s. Totalcount: Deviation-in-degree over iCEGs ###############################

# Cosbin + COTIN imputed
cosbin_super_sample <- NULL
start <- 1
for(rep in nRep){
  cosbin_super_sample <- cbind(cosbin_super_sample, rowMeans(COTIN_imputed[, start : (start + rep - 1)])) 
  # data 3 will have a column number the same as your group number
  start <- start + rep
}

# Cosbin functions
cosbin_out <- cosbin(cosbin_super_sample)
cosbin_out_full <- cosbin_convert(cosbin_out, COTIN_imputed) 

# Totalcount + COTIN imputed
totalcount_out<-totalcount(COTIN_imputed)$data


# Compare the Deviation-in-degree over iCEGs
grnd_iCEG<-data_grnd_full[ind_CEG,]

cosbin_iCEG<-cosbin_out_full[ind_CEG,]
totalcount_iCEG<-totalcount_out[ind_CEG,]


# Cosbin's
Deviation_in_degree_iCEG<-NULL
for (i in 1:dim(grnd_iCEG)[1]){
  temp_cosine<-mycosine(cosbin_iCEG[i,],grnd_iCEG[i,])
  Deviation_in_degree_iCEG<-c(Deviation_in_degree_iCEG,(180 * acos(temp_cosine)) / pi)
}

Cosbin_deviation_iCEG<-sum(Deviation_in_degree_iCEG)
Averaged_Cosbin_deviation_iCEG<-Cosbin_deviation_iCEG/dim(grnd_iCEG)[1]

Averaged_Cosbin_iCEG_error<-euclidean(cosbin_iCEG,data_grnd_full[ind_CEG,])/length(ind_CEG)

# Totalcount's
Deviation_in_degree_iCEG<-NULL
for (i in 1:dim(grnd_iCEG)[1]){
  temp_cosine<-mycosine(totalcount_iCEG[i,],grnd_iCEG[i,])
  Deviation_in_degree_iCEG<-c(Deviation_in_degree_iCEG,(180 * acos(temp_cosine)) / pi)
}

Totalcount_deviation_iCEG<-sum(Deviation_in_degree_iCEG)
Averaged_Totalcount_deviation_iCEG<-Totalcount_deviation_iCEG/dim(grnd_iCEG)[1]

Averaged_Totalcount_iCEG_error<-euclidean(totalcount_iCEG,data_grnd_full[ind_CEG,])/length(ind_CEG)




################################   Scatter Simplex Plot   #####################################


### COTIN first
COTIN_super_sample_simplex_plot <- NULL
start <- 1
for(rep in nRep){
  COTIN_super_sample_simplex_plot <- cbind(COTIN_super_sample_simplex_plot, 
                                           rowMeans(COTIN_imputed[, start : (start + rep - 1)])) 
  start <- start + rep
}

ind_CEG <- which(cos_iCEG(COTIN_super_sample_simplex_plot) >= iCEG_threshold)

ref <- vector("list", length = nGroup)
for (i in 1:nGroup) {
  ref[[i]] <- rep(0, nGroup)
  ref[[i]][i] <- 1
}


ind_sDEG2 <- list()
for (i in 1:nGroup) {
  ind_sDEG2[[i]] <- which(cos_ref(COTIN_super_sample_simplex_plot, ref[[i]]) >= SG_threshold)
}

sum(lengths(ind_sDEG2)) # This gives the number of identified sDEGs in the data

lbl <- rep(0, nGene)
lbl[ind_CEG] <- 'iCEG (labeled)'

iCEG_COTIN<-length(lbl[ind_CEG]) # This gives the number of identified iCEGs in the data

lbl2 <- lbl
lbl3 <- lbl
for (i in 1:nGroup) {
  lbl2[ind_sDEG2[[i]]] <- paste0('sDEG G', i)
  lbl3[ind_sDEG2[[i]]] <- i
}

lbl2[lbl2=="sDEG"]<-'DEG'


X <- COTIN_super_sample_simplex_plot
Xproj <- X
Xproj <- apply(Xproj, 1, function(x) x / sum(x))

A <- diag(dim(Xproj)[1])
K = dim(A)[2]
PS <- t(matrix(c(cos((seq_len(K) - 1) * 2 * pi / K),
                 sin((seq_len(K) - 1) * 2 * pi / K)), K))
library(corpcor)
mg <-list()
for (i in 1:nGroup){
  mg[[i]]<-which(lbl3==i)
}

mg_CONTI<-mg

tmp <- PS %*% pseudoinverse(A)
Xp <- tmp %*% Xproj

# color-coding for SGs
plot(t(PS), xlab = 'COTIN Imputed', ylab = NA, asp = 1, pch = '.')
points(Xp[1,], Xp[2,], col = rgb(80, 80, 80, max = 255, alpha = 125),
       xlab = NA, ylab = NA, asp = 1)

mg.col = c('dark green','orange','yellow','dodgerblue','blue')
for(i in seq_along(mg)){
  points(Xp[1,mg[[i]]], Xp[2,mg[[i]]],
         col=mg.col[i])
}
points(t(PS), pch = 17)

# color-coding for iCEGs
# plot(t(PS), xlab = NA, ylab = NA, asp = 1, pch = '.')
# points(Xp[1,], Xp[2,], col = rgb(80, 80, 80, max = 255, alpha = 125),
#        xlab = NA, ylab = NA, asp = 1)

# for(i in 1:length(ind_CEG)){
#   points(Xp[1,ind_CEG[i]], Xp[2,ind_CEG[[i]]],
#          col='green')
# }
# points(t(PS), pch = 17)

mg_temp<-intersect(index_masked_out,unlist(mg_grnd))
mg_intersect_COTIN<-intersect(unlist(mg_CONTI),mg_temp)

# color-coding for SGs
plot(t(PS), xlab = 'Masked-out genes recovered by COTIN', ylab = NA, asp = 1, pch = '.')
points(Xp[1,], Xp[2,], col = rgb(80, 80, 80, max = 255, alpha = 125),
       xlab = NA, ylab = NA, asp = 1)

mg.col = c('dark green','orange','yellow','dodgerblue','blue')
# for(i in seq_along(mg_grnd)){
#   points(Xp[1,mg_grnd[[i]]], Xp[2,mg_grnd[[i]]],
#          col=mg.col[i])
# }
# points(t(PS), pch = 17)

for(i in 1:length(mg_grnd_masked_out)){
  points(Xp[1,mg_grnd_masked_out[i]], Xp[2,mg_grnd_masked_out[i]],
         col='red')
}
points(t(PS), pch = 17)

for(i in 1:length(non_SG_masked_out)){
  points(Xp[1,non_SG_masked_out[i]], Xp[2,non_SG_masked_out[i]],
         col='pink')
}
points(t(PS), pch = 17)

# legend('bottomright',"Legend",c("Group 1 SG","Group 2 SG","Group 3 SG","Group 4 SG","Group 5 SG","Retrived masked-out SG",
#                                 "Retrived masked-out non-SG"),
#        fill=c('dark green','orange','yellow','dodgerblue','blue','red','pink'),cex=0.42,pt.cex = 1)

legend('bottomright',"Legend",c("Retrived masked-out SG",
                                "Retrived masked-out non-SG"),
       fill=c('red','pink'),cex=0.42,pt.cex = 1)



### Min/2 imputed

ind_sDEG2 <- list()
for (i in 1:nGroup) {
  ind_sDEG2[[i]] <- which(cos_ref(Half_min_imputed_super_sample, ref[[i]]) >= SG_threshold)
}

sum(lengths(ind_sDEG2)) # This gives the number of identified sDEGs in the data

lbl <- rep(0, nGene)
lbl[ind_CEG] <- 'iCEG (labeled)'

iCEG_Half_min<-length(lbl[ind_CEG]) # This gives the number of identified iCEGs in the data

lbl2 <- lbl
lbl3 <- lbl
for (i in 1:nGroup) {
  lbl2[ind_sDEG2[[i]]] <- paste0('sDEG G', i)
  lbl3[ind_sDEG2[[i]]] <- i
}

lbl2[lbl2=="sDEG"]<-'DEG'




X <- Half_min_imputed_super_sample
Xproj <- X
Xproj <- apply(Xproj, 1, function(x) x / sum(x))

A <- diag(dim(Xproj)[1])
K = dim(A)[2]
PS <- t(matrix(c(cos((seq_len(K) - 1) * 2 * pi / K),
                 sin((seq_len(K) - 1) * 2 * pi / K)), K))
library(corpcor)
mg <-list()
for (i in 1:nGroup){
  mg[[i]]<-which(lbl3==i)
}

mg_Half_min<-mg

tmp <- PS %*% pseudoinverse(A)
Xp <- tmp %*% Xproj

# color-coding for SGs
plot(t(PS), xlab = 'Half Minimum Imputed', ylab = NA, asp = 1, pch = '.')
points(Xp[1,], Xp[2,], col = rgb(80, 80, 80, max = 255, alpha = 125),
       xlab = NA, ylab = NA, asp = 1)

mg.col = c('dark green','orange','yellow','dodgerblue','blue')
for(i in seq_along(mg)){
  points(Xp[1,mg[[i]]], Xp[2,mg[[i]]],
         col=mg.col[i])
}
points(t(PS), pch = 17)

# color-coding for iCEGs
# plot(t(PS), xlab = NA, ylab = NA, asp = 1, pch = '.')
# points(Xp[1,], Xp[2,], col = rgb(80, 80, 80, max = 255, alpha = 125),
#        xlab = NA, ylab = NA, asp = 1)
# 
# for(i in 1:length(ind_CEG)){
#   points(Xp[1,ind_CEG[i]], Xp[2,ind_CEG[[i]]],
#          col='green')
# }
# points(t(PS), pch = 17)


mg_temp<-intersect(index_masked_out,unlist(mg_grnd))
mg_intersect_Half_min<-intersect(unlist(mg_Half_min),mg_temp)

# color-coding for SGs
plot(t(PS), xlab = 'Masked-out genes recovered by half minimum', ylab = NA, asp = 1, pch = '.')
points(Xp[1,], Xp[2,], col = rgb(80, 80, 80, max = 255, alpha = 125),
       xlab = NA, ylab = NA, asp = 1)

# mg.col = c('dark green','orange','yellow','dodgerblue','blue')
# for(i in seq_along(mg_grnd)){
#   points(Xp[1,mg_grnd[[i]]], Xp[2,mg_grnd[[i]]],
#          col=mg.col[i])
# }
# points(t(PS), pch = 17)

for(i in 1:length(mg_grnd_masked_out)){
  points(Xp[1,mg_grnd_masked_out[i]], Xp[2,mg_grnd_masked_out[i]],
         col='red')
}
points(t(PS), pch = 17)

for(i in 1:length(non_SG_masked_out)){
  points(Xp[1,non_SG_masked_out[i]], Xp[2,non_SG_masked_out[i]],
         col='pink')
}
points(t(PS), pch = 17)

# legend('bottomright',"Legend",c("Group 1 SG","Group 2 SG","Group 3 SG","Group 4 SG","Group 5 SG","Retrived masked-out SG",
#                                 "Retrived masked-out non-SG"),
#        fill=c('dark green','orange','yellow','dodgerblue','blue','red','pink'),cex=0.42,pt.cex = 1)

legend('bottomright',"Legend",c("Retrived masked-out SG",
                                "Retrived masked-out non-SG"),
       fill=c('red','pink'),cex=0.42,pt.cex = 1)



### Mean imputed

ind_sDEG2 <- list()
for (i in 1:nGroup) {
  ind_sDEG2[[i]] <- which(cos_ref(Mean_imputed_super_sample, ref[[i]]) >= SG_threshold)
}

sum(lengths(ind_sDEG2)) # This gives the number of identified sDEGs in the data

lbl <- rep(0, nGene)
lbl[ind_CEG] <- 'iCEG (labeled)'

iCEG_Mean<-length(lbl[ind_CEG]) # This gives the number of identified iCEGs in the data

lbl2 <- lbl
lbl3 <- lbl
for (i in 1:nGroup) {
  lbl2[ind_sDEG2[[i]]] <- paste0('sDEG G', i)
  lbl3[ind_sDEG2[[i]]] <- i
}

lbl2[lbl2=="sDEG"]<-'DEG'




X <- Mean_imputed_super_sample
Xproj <- X
Xproj <- apply(Xproj, 1, function(x) x / sum(x))

A <- diag(dim(Xproj)[1])
K = dim(A)[2]
PS <- t(matrix(c(cos((seq_len(K) - 1) * 2 * pi / K),
                 sin((seq_len(K) - 1) * 2 * pi / K)), K))
library(corpcor)
mg <-list()
for (i in 1:nGroup){
  mg[[i]]<-which(lbl3==i)
}

mg_Mean<-mg

tmp <- PS %*% pseudoinverse(A)
Xp <- tmp %*% Xproj

# color-coding for SGs
plot(t(PS), xlab = 'Mean Imputed', ylab = NA, asp = 1, pch = '.')
points(Xp[1,], Xp[2,], col = rgb(80, 80, 80, max = 255, alpha = 125),
       xlab = NA, ylab = NA, asp = 1)

mg.col = c('dark green','orange','yellow','dodgerblue','blue') # using the ground truth labels
for(i in seq_along(mg_grnd)){
  points(Xp[1,mg_grnd[[i]]], Xp[2,mg_grnd[[i]]],
         col=mg.col[i])
}
points(t(PS), pch = 17)

# color-coding for iCEGs
# plot(t(PS), xlab = NA, ylab = NA, asp = 1, pch = '.')
# points(Xp[1,], Xp[2,], col = rgb(80, 80, 80, max = 255, alpha = 125),
#        xlab = NA, ylab = NA, asp = 1)
# 
# for(i in 1:length(ind_CEG)){
#   points(Xp[1,ind_CEG[i]], Xp[2,ind_CEG[[i]]],
#          col='green')
# }
# points(t(PS), pch = 17)


mg_temp<-intersect(index_masked_out,unlist(mg_grnd))
mg_intersect_Mean<-intersect(unlist(mg_Mean),mg_temp)

# color-coding for SGs
plot(t(PS), xlab = 'Masked-out genes recovered by mean', ylab = NA, asp = 1, pch = '.')
points(Xp[1,], Xp[2,], col = rgb(80, 80, 80, max = 255, alpha = 125),
       xlab = NA, ylab = NA, asp = 1)


for(i in 1:length(mg_grnd_masked_out)){
  points(Xp[1,mg_grnd_masked_out[i]], Xp[2,mg_grnd_masked_out[i]],
         col='red')
}
points(t(PS), pch = 17)

for(i in 1:length(non_SG_masked_out)){
  points(Xp[1,non_SG_masked_out[i]], Xp[2,non_SG_masked_out[i]],
         col='pink')
}
points(t(PS), pch = 17)

legend('bottomright',"Legend",c("Retrived masked-out SG",
                                "Retrived masked-out non-SG"),
       fill=c('red','pink'),cex=0.42,pt.cex = 1)



######################################### Report the quantitative Results ######################################

# Mechanism-specific missing rates
MAR/(dim(before_norm)[1]*dim(before_norm)[2]) # MAR rate
(Total_missing-MAR)/(dim(before_norm)[1]*dim(before_norm)[2]) # LLOD rate

### Comparison to Min/2 and Mean
missing_rate_threshold # Threshold for masking-out
1000-dim(after_remove)[1] # Number of masked-out features (may include non-SGs)
length(non_SG_masked_out) # Number of masked-out non-SGs
length(mg_grnd_masked_out) # Number of masked-out SGs
sum(after_missing_masked_out == 0) # Number of missing values in the masked-out SGs

# Euclidean distance for imputing the missing values in the masked-out SGs
SG_Mean_error
SG_Half_min_error
SG_CONTI_error

Averaged_SG_Mean_error
Averaged_SG_Half_min_error
Averaged_SG_CONTI_error

# Deviation of masked-out SGs
Averaged_Mean_Deviation_in_degree_SG
Averaged_Half_min_Deviation_in_degree_SG
Averaged_COTIN_Deviation_in_degree_SG

# What about masked-out non-SG
Averaged_non_SG_Half_min_error
Averaged_non_SG_Mean_error
Averaged_non_SG_COTIN_error

Averaged_Half_min_Deviation_in_degree_non_SG
Averaged_Mean_Deviation_in_degree_non_SG
Averaged_COTIN_Deviation_in_degree_non_SG

# Together: masked-out SGs + non-SGs
COTIN_together_error
Mean_together_error
Half_min_together_error

Averaged_COTIN_together_error
Averaged_Mean_together_error
Averaged_Half_min_together_error

Averaged_COTIN_Deviation_in_degree
Averaged_Mean_Deviation_in_degree
Averaged_Half_min_Deviation_in_degree

# Deviation of iCEGs
number_iCEG_grnd
Cosbin_deviation_iCEG
Averaged_Cosbin_deviation_iCEG

Totalcount_deviation_iCEG
Averaged_Totalcount_deviation_iCEG

Averaged_Cosbin_iCEG_error
Averaged_Totalcount_iCEG_error

# Checking iCEG number
number_iCEG_grnd
iCEG_COTIN
iCEG_Half_min
iCEG_Mean

# RMSE & NRMSE
rmse_CONTI_sg
rmse_CONTI_non_sg
rmse_CONTI_all
nrmse_CONTI_sg
nrmse_CONTI_non_sg
nrmse_CONTI_all

rmse_Half_min_sg
rmse_Half_min_non_sg
rmse_Half_min_all
nrmse_Half_min_sg
nrmse_Half_min_non_sg
nrmse_Half_min_all

rmse_Mean_sg
rmse_Mean_non_sg
rmse_Mean_all
nrmse_Mean_sg
nrmse_Mean_non_sg
nrmse_Mean_all

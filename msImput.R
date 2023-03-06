#' @param readin Read-in data, with features on the row and samples (grouped) on the columns. All missing values (NA) should 
#' be imputed beforehand.
#' @param nRep The number of samples in each group, e.g., c(5,4,3,2), corresponding to the readin data.
#' @param min_option The options for the minimum value, with 'global' (default) meaning to use the global minimum value across the 
#' whole data matrix, or 'local' meaning to use the feature-wise minimum value.
#' @param missing_rate_threshold A value between 0 and 1. The upper threshold of sample-wise missing-rate allowed for features 
#' in any of each group, with 1 as default (all features are allowed).
#' @param sd_scaler For scaling up the standard deviation of the distribution obtained from the non-missing value. Default
#' is 1 (no scaling).
#' 
#' @return the normalized data (with selected features. If you want all features included, set missing_rate_threshold to 1).
#' 
#' @example 
#' msImputed <- msImput(readin=data, nRep=c(5,4,3,2), min_option='global', missing_rate_threshold=0.8, sd_scaler=5)
#' 
#' @End

msImput<-function(readin,nRep,min_option='global',missing_rate_threshold=1, sd_scaler=1){
  
  before_remove<-readin
  rownames<-row.names(readin)
  
  message("msImput: Pre-processing")
  
  # turn characters into numbers
  
  for (i in 1:dim(before_remove)[2]){   
    temp<-as.numeric(before_remove[,i])
    before_remove[,i]<-temp
  }
 
  
  
  # temporarily replace NA with 0
  
  before_remove[is.na(before_remove)]<-0   
  after_remove<-NULL 
  index_name<-NULL
  
  
  # select features with missing rate less than ?% in any of each group
  
  message("msImput: Selecting features with missing rate less than the threshold you choose in any of each group")
  
  nRep_sum<-vector() ## suppose nRep is 21 39 34 2, nRep_sum is 21 60 94 96
  for (i in 1:length(nRep)){
    nRep_sum[i]<-sum(nRep[1:i])
  }
  nRep_plus<-c(0,nRep_sum)  ## nRep_plus will be 0 21 60 94 96
  non_0<-list() ## this stores the number of non-0 values in each group
  
  for (j in 1:dim(before_remove)[1]){
    for (i in 1:length(nRep)){
      non_0[i]<-sum(before_remove[j,(nRep_plus[i]+1):(nRep_plus[i+1])] !=0)
    }
    
    if (any(as.matrix(non_0)>as.matrix(nRep)*(1-missing_rate_threshold))){ 
      after_remove<-rbind(after_remove,before_remove[j,])
      index_name<-c(index_name,j)
    }
  }
  
  before_impute<-after_remove
  
  row.names(before_impute)<-rownames[index_name]
  
  
  # Imputation
  
  after_impute<-NULL
  protein_segment<-list()
  log_segment<-list()
  overall_LLOD_cells<-0
  
  if (min_option == 'local'){
    message("msImput: Using weighted feature-wise minimum value & group-mean for feature-wise pre-imputation")
    
    for (j in 1:dim(before_impute)[1]){
      
      this_protein<-before_impute[j,]
      
      if (length(this_protein)==length(this_protein[this_protein==0])){
        message("msImput: Features with all value missing detected! Skip to the next feature. Consider removing.")
        after_impute<-rbind(after_impute,this_protein)
        overall_LLOD_cells<-overall_LLOD_cells+length(this_protein)
        
      } else {
      
        feature_min<-min(this_protein[this_protein!=0])
        feature_min_log<-log(feature_min)
        
        for (i in 1:length(nRep)){
          non_0[i]<-sum(before_impute[j,(nRep_plus[i]+1):(nRep_plus[i+1])] !=0)
          protein_segment[[i]]<-before_impute[j,(nRep_plus[i]+1):(nRep_plus[i+1])]

          if (non_0[i]==0){
            protein_segment[[i]][protein_segment[[i]]==0]<-feature_min/2
            missing_LLOD_rate<-1
            
          } else if (non_0[i]==1){
            protein_segment[[i]][protein_segment[[i]]==0]<-feature_min/2
            missing_LLOD_rate<-(length(protein_segment[[i]])-1)/length(protein_segment[[i]])
            
          } else {
            log_segment[[i]]<-log(protein_segment[[i]][protein_segment[[i]]!=0])
            missing_LLOD_rate<-pnorm(feature_min_log,mean = mean(log_segment[[i]]), sd = sd_scaler*sd(log_segment[[i]]))
            protein_segment[[i]][protein_segment[[i]]==0]<-
              missing_LLOD_rate * feature_min/2 + 
              (1 - missing_LLOD_rate) * mean(protein_segment[[i]][protein_segment[[i]]!=0])
          }
          overall_LLOD_cells<-overall_LLOD_cells + missing_LLOD_rate*length(protein_segment[[i]])
        }
      
        after_impute<-rbind(after_impute,unlist(protein_segment))
      }
    }
  }
  
  
  
  if (min_option == 'global'){
    message("msImput: Using weighted global minimum value & feature-wise group-mean for feature-wise pre-imputation")
    for (j in 1:dim(before_impute)[1]){
      
      this_protein<-before_impute[j,]
      
      if (length(this_protein)==length(this_protein[this_protein==0])){
        message("msImput: Features with all value missing detected! Skip to the next feature. Consider removing.")
        after_impute<-rbind(after_impute,this_protein)
        overall_LLOD_cells<-overall_LLOD_cells+length(this_protein)
        
      } else {
        
        global_min<-min(before_impute[before_impute!=0])
        global_min_log<-log(global_min)
        
        for (i in 1:length(nRep)){
          non_0[i]<-sum(before_impute[j,(nRep_plus[i]+1):(nRep_plus[i+1])] !=0)
          protein_segment[[i]]<-before_impute[j,(nRep_plus[i]+1):(nRep_plus[i+1])]
        }
        
        for (i in 1:length(nRep)){
          if (non_0[i]==0){
            protein_segment[[i]][protein_segment[[i]]==0]<-global_min/2
            missing_LLOD_rate<-1
            
          } else if (non_0[i]==1){
            protein_segment[[i]][protein_segment[[i]]==0]<-global_min/2
            missing_LLOD_rate<-(length(protein_segment[[i]])-1)/length(protein_segment[[i]])
            
          } else {
            log_segment[[i]]<-log(protein_segment[[i]][protein_segment[[i]]!=0])
            missing_LLOD_rate<-pnorm(global_min_log,mean = mean(log_segment[[i]]), sd = sd_scaler*sd(log_segment[[i]]))
            protein_segment[[i]][protein_segment[[i]]==0]<-
              missing_LLOD_rate * global_min/2 + 
              (1 - missing_LLOD_rate) * mean(protein_segment[[i]][protein_segment[[i]]!=0])
          }
          overall_LLOD_cells<-overall_LLOD_cells + missing_LLOD_rate*length(protein_segment[[i]])
        }
        
        after_impute<-rbind(after_impute,unlist(protein_segment))
      }
    }
  }
  

  row.names(after_impute)<-rownames[index_name]

  message("msImput: Done.")
  
  return(after_impute)
  
}
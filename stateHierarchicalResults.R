# setwd("~/Desktop/rep-5-n-247-357/result_dir/")
args <- commandArgs()
otu_table_file <- args[6]
meta_table_file <- args[7]
kldm_res_dir <- args[8]
attr_file <- args[9]
target_dir <- paste(kldm_res_dir, "/hierarchical_results/", sep = "")
# otu_table_file <- "../otu_table"
# meta_table_file <- "../meta_table"
# target_dir <- "./hierarchical_results/"
# attr_file <- "../attr_table"
# meta_name_file <- "../meta_name_file"
# otu_name_file <- "../otu_name_file"

otu_table <- as.matrix(read.table(file = otu_table_file, sep = " "))
otu_table_r <- log(otu_table + 1) / rowSums(log(otu_table + 1))
meta_table <- as.matrix(read.table(file = meta_table_file, sep =  " "))
attr_table <- as.matrix(read.table(file = attr_file, sep = " "))
# meta_name <- read.table(file = meta_name_file)
# otu_name <- read.table(file = otu_name_file)

attr_sums <- colSums(attr_table)

library(stringr)

saveInfo <- function(merge_dir, stat_info, asso_info, basic_info) {
  meta_mean <- stat_info[[1]]
  meta_sd <- stat_info[[2]]
  otu_mean <- stat_info[[3]]
  otu_sd <- stat_info[[4]]
  attr_mean1 <- stat_info[[5]]
  attr_mean2 <- stat_info[[6]]
  
  # mkdir stat
  stat_dir <- paste(merge_dir, "/stat/", sep = "")
  dir.create(stat_dir)
  
  meta_file <- paste(stat_dir, "/meta_info", sep = "")
  otu_file <- paste(stat_dir, "/otu_info", sep = "")
  attr_file <- paste(stat_dir, "/attr_info", sep = "")
  
  meta_info <- rbind(meta_mean, meta_sd)
  otu_info <- rbind(otu_mean, otu_sd)
  attr_info <- rbind(attr_mean1, attr_mean2)
  
  write.table(x = meta_info, file = meta_file, sep = " ", row.names = FALSE, col.names = FALSE)
  write.table(x = otu_info, file = otu_file, sep = " ", row.names = FALSE, col.names = FALSE)
  write.table(x = attr_info, file = attr_file, sep = " ", row.names = FALSE, col.names = FALSE)
  
  otu_otu_m <- asso_info[[1]]
  ef_otu_m <- asso_info[[2]]
  otu_otu_file <- paste(stat_dir, "/otu_otu_association", sep = "")
  ef_otu_file <- paste(stat_dir, "/ef_otu_association", sep = "")
  write.table(x = otu_otu_m, file = otu_otu_file, sep = " ", row.names = FALSE, col.names = FALSE)
  write.table(x = ef_otu_m, file = ef_otu_file, sep = " ", row.names = FALSE, col.names = FALSE)
  
  per_sample_size <- basic_info[1]
  per_otu_size <- basic_info[2]
  
  basic_file <- paste(stat_dir, "/basic_info", sep = "")
  write.table(x = basic_info, file = basic_file, sep = "\n", row.names = FALSE, col.names = FALSE)
  
}

convertAssociation <- function(otu_index, otu_otu, ef_otu) {
  p <- ncol(otu_table)
  q <- ncol(meta_table)
  
  otu_otu_full <- matrix(0, p, p)
  ef_otu_full <- matrix(0, q, p)
  
  per_p <- nrow(otu_otu)
  for (i in 1:per_p) {
    index1 <- otu_index[i]
    for (j in 1:per_p) {
      index2 <- otu_index[j]
      otu_otu_full[index1, index2] <- otu_otu[i, j]
    }
    
    for (k in 1:q) {
      ef_otu_full[k, index1] <- ef_otu[k, i]
    }
  }
  
  return(list(otu_otu_full, ef_otu_full))
}

parseDirRec <- function(per_dir) {
  # merged cluster
  per_merge_dir <- paste(per_dir, "/merged/", sep = "")
  if (dir.exists(per_merge_dir)) {
    print(per_merge_dir)
    per_result_list <- ImportCluster(per_merge_dir, 1)
    per_sample_index <- unlist(per_result_list[[9]]) + 1
    per_otu_index <- unlist(per_result_list[[1]]) + 1
    # per_meta_mean <- per_result_list[[3]]
    per_stat <- statisticFeature(per_sample_index)
    
    per_otu_otu <- as.matrix(per_result_list[[7]])
    per_ef_otu <- as.matrix(per_result_list[[8]])
    per_asso_stat <- convertAssociation(per_otu_index, per_otu_otu, per_ef_otu)
    
    per_basic <- c(length(per_sample_index), length(per_otu_index))
    
    saveInfo(per_merge_dir, per_stat, per_asso_stat, per_basic)
  }
  
  # add cluster 
  per_add_dir <- paste(per_dir, "/added/", sep = "")
  if (dir.exists(per_add_dir)) {
    parseDirRec(per_add_dir)
  }
  
  # origin cluster
  per_origin_dir <- paste(per_dir, "/original/", sep = "")
  if (dir.exists(per_origin_dir)) {
    parseDirRec(per_origin_dir)
  }
}

statisticFeature <- function(sample_index) {
  sub_otu_r <- otu_table_r[sample_index,]
  sub_meta <- meta_table[sample_index,]
  sub_attr <- attr_table[sample_index,]
  
  meta_mean <- colMeans(sub_meta)
  meta_sd <- apply(X = sub_meta, MARGIN = 2, FUN = sd)
  otu_r_mean <- colMeans(sub_otu_r)
  otu_sd <- apply(X = sub_otu_r, MARGIN = 2, FUN = sd)
  
  sub_attr_sums <- colSums(sub_attr)
  sub_size <- length(sample_index)
  
  attr_ratio1 <- sub_attr_sums / sub_size
  attr_ratio2 <- sub_attr_sums / attr_sums
  
  per_stat <- list(meta_mean, meta_sd, otu_r_mean, otu_sd, attr_ratio1, attr_ratio2)
  return(per_stat)
}

ImportCluster <- function(result_dir, i) {
  # otu list
  otu_index_file <- paste(result_dir, "/OTU_Index_", as.character(i), sep = "")
  otu_index <- read.table(file = otu_index_file, header = FALSE, sep = " ")[1,]
  
  # B
  b_file <- paste(result_dir, "/B_", as.character(i), sep = "")
  B <- read.table(file = b_file, header = FALSE, sep = " ")
  
  # theta
  theta_file <- paste(result_dir, "/Theta_", as.character(i), sep = "")
  Theta <- read.table(file = theta_file, header = FALSE, sep = " ")
  
  # otu otu association
  association1_file <- paste(result_dir, "/OTU_OTU_Association_", as.character(i), sep = "")
  otu_otu_association <- read.table(file = association1_file, header = FALSE, sep = " ")
  
  # meta otu association
  association2_file <- paste(result_dir, "/Meta_OTU_Association_", as.character(i), sep = "")
  meta_otu_association <- read.table(file = association2_file, header = FALSE, sep = " ")
  
  # meta mean
  meta_mean_file <- paste(result_dir, "/Meta_Mean_", as.character(i), sep = "")
  meta_mean <- read.table(file = meta_mean_file, header = FALSE, sep = " ")[1,]
  
  # meta cov
  meta_cov_file <- paste(result_dir, "/Meta_Cov_", as.character(i), sep = "")
  meta_cov <- read.table(file = meta_cov_file, header = FALSE, sep = " ")
  
  # sample num
  num_file <- paste(result_dir, "/Sample_Num_", as.character(i), sep = "")
  sample_num <- read.table(file = num_file, header = FALSE, sep = " ")[1,1]
  
  # sample index
  sample_index_file <- paste(result_dir, "/Sample_Index_", as.character(i), sep = "")
  sample_index <- read.table(file = sample_index_file, header = FALSE, sep = " ")[1,]
  
  per_list <- list(otu_index, sample_num, meta_mean, meta_cov, B, Theta, otu_otu_association, meta_otu_association, sample_index)

  return(per_list)
}

for (per_dir in list.dirs(target_dir, recursive = FALSE)) {
  per_number <- as.numeric(str_match(per_dir, "-(\\d+)")[[2]])
  
  parseDirRec(per_dir)
}


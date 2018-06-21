# setwd("~/Desktop/")
args <- commandArgs()
# result directory of kmldm
# result_dir <- "./result_dir/"
result_dir <- args[6]
# obtain the cluster number
number_file <- paste(result_dir, "/Cluster_Number", sep = "")
cluster_number <- read.table(file = number_file, header = FALSE, sep = " ")[1,1]
# print(cluster_number)
parameters <- list()
# import kmldm result
for (i in 1:cluster_number) {
  # print(i)
  per_result <- list()
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
  parameters[[i]] <- per_list
  # print("end")
}

data_file <- paste(result_dir, "/parameters.RData", sep = "")
save(parameters, file = data_file)

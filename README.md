# kLDM [![language](https://img.shields.io/badge/language-c++-brightgreen.svg)]() [![version](https://img.shields.io/badge/version-v1.0-blue.svg)]() [![topic](https://img.shields.io/badge/metagenomics-association_inference-00CED1.svg)]()
kLDM is designed to infer multiple association networks based on variation of environmental factors. 

![](https://github.com/tinglab/kLDM/blob/master/pictures/model_explain.jpg)

kLDM assumes that interactions among microbes are regulated by environmental factors and they tend to be stable when environmental conditions are identical but may change with variation in values of environmental factors. Under the same environmental condition, values of environmental factors are fluctuant in a small range, while core microbes are the same, and associations within microbial community are stable. Between different environmental conditions, values of environmental factors, species of microbes and their associations can be different.

## Pre-requisites
- gcc 4.8 (support c++11 and openmp)
- make
- lapack, blas library
```
# centos
sudo yum install lapack lapack-devel blas blas-devel
# or ubuntu
# sudo apt-get install libblas-dev
# sudo apt-get install liblapack-dev
```

## Installation
Using the kLDM in Linux OS servers is relatively easy.
### Linux OS
Just compile source and the executable file 'run' will be generated.
- The Eigen library for linear algebra was put into the directory 'kLDM/eigen', which will be used during the compiling process ~ If the Eigen library has been installed into the computer, You can change the Makefile to use own Eigen.
```
# compile kLDM
cd ./src
make
```

## Usage
### Input
To run kLDM, you need to specify four parameters:
* Three input files: otu_table, meta_table and matrix_shape

File Name | Content
----------|--------
otu_table | N * P matrix, every line is a sample and P microbes' counts are seperated by a spacing (' ').
meta_table | N * Q matrix like otu_table
matrix_shape | three integers: N, P, Q, which are seperated by line feed ('\n').
    
* Output directory
### Example
* Commands
```
cd ./datasets/crc-data/
../../src/run ./crc_otu_table ./crc_fit_meta_table ./crc_matrix_shape ./result
```
* Estimated parameters will be saved into 'result' directory as below:
![](https://github.com/tinglab/kLDM/blob/master/pictures/result_example.png)
### Output
kLDM estimates EF conditions and association networks under every EF condition via a split-merge algorithm:
![](https://github.com/tinglab/kLDM/blob/master/pictures/sm-process.jpg)
* Output Explaination in the 'result' directory

File Name | Content
----------|--------
Cluster_Number | the number of estimated EF conditions are saved into the file.
Sample_Index_k (e.g. k=1, k=2) | indexes of samples belongs to the cluster k in original 'otu_table' (the minimum index is 0).
Sample_Num_k | the number of samples of the cluster k.
OTU_Index_k | indexes of OTUs belongs to the cluster k in original 'otu_table'. **:pill: Not all OTUs will be included in the cluster k, and OTUs' order can be changed in different clusters.**
Meta_OTU_Association_k | estimated EF-OTU associations of the cluster k.
OTU_OTU_Association_k | estimated OTU-OTU associatons of the cluster k.
B_k | records the matrix which parameterizes EF-OTU associations of the cluster k.
Theta_k | records the matrix which parameterizes OTU-OTU associations of the cluster k.
Meta_Mean_k | records the mean values of EFs of the cluster k.
Meta_Cov_k | records the covariance of EFs of the cluster k.

**:smiling_imp: WARNING: the minimum index is 0 in all index files ('Sample_Index_k', 'OTU_Index_k'), if you process them in R language, remember to add one (+1) to get the right index.**

* You can use the script 'ImportKmldmResult.R' to combine files within the 'result' to produce a 'parameters.RData' file, which can be processed in RStudio conveniently.
```
../../ImportKmldmResult.R ./result
```

## Datasets
Three datasets used in the paper.

Name | Directory | Description | Link
-----|-----------|------------ | ----
Colorectal cancer gut dataset | ./datasets/crc-data/ | 16s rRNA gene sequencing, 490 samples with 117 OTUs and 5 EFs | [Baxter, et al. (2016)](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-016-0290-3) 
TARA Oceans eukaryotic dataset | ./datasets/tara-data/ | 18s rRNA gene sequencing, 221 samples with 67 OTUs and 17 EFs | [Limamendez, et al., 2015](http://science.sciencemag.org/content/348/6237/1262073)
American Gut project dataset | ./datasets/american-gut-data | 16s rRNA gene sequencing, 11946 samples with 216 OTUs and 22 EFs related to human lifestyle | ftp://ftp.microbio.me/AumericanGut, [Official Website](http://americangut.org/)

## .Biom file processing
The script './processBiom.py' can be used to filter OTUs and Samples in the .biom file.

## :dart: New Feature -- Saving Intermediate Clustering Process and Visualization
### Intermediate Merge Process
Except the final estimated EF-condtions, or clusters, kLDM will also record sub-clusters and the merge process in the directory '{output_directory}/hierarchical_results/'. Every time kLDM will merge clusters from 'added/merged' and 'original/merged' into one cluster, which is saved into new 'merged' directory. The structure of directories is shown as below, 
```
|-- results
      |-- hierarchical_results
            |-- cluster-1
            |        |-- added
            |        |       |-- added
            |        |       |-- merged
            |        |       |-- original
            |        |
            |        |-- merged
            |        |       |-- B_1
            |        |       |-- Meta_Cov_1
            |        |       |-- Meta_Mean_1
            |        |       |-- Meta_OTU_Association_1
            |        |       |-- OTU_Index_1
            |        |       |-- OTU_OTU_Association_1
            |        |       |-- Sample_Index_1
            |        |       |-- Sample_Num_1
            |        |       |-- Theta_1
            |        |       
            |        |-- original
            |               |-- added
            |               |-- merged
            |               |-- original
            |        
            |-- cluster-2
            |-- cluster-3
            ...
```

### Hierarchical Server
In order to visualize all clusters and be convenient for comparing abundance and OTU-OTU and EF-OTU associations among clusters, a web application is developed based python Flask framework and its code is at 'hierarchical_server' directory.
* Pre-requisites (python packages): flask, argparse, numpy, scipy
* How to run
```
cd ./datasets/
Rscript ../stateHierarchicalResults.R ./american-gut-data/otu_table ./american-gut-data/meta_table ./american-gut-data/result/ ./american-gut-data/attr_table
python ../hierarchical_server/hierarchical_server.py --on ./american-gut-data/otu_annotation.txt -mn ./american-gut-data/meta_annotation.txt -an ./american-gut-data/attr_name -ms ./american-gut-data/matrix_shape -res ./american-gut-data/result/ -ot ./american-gut-data/otu_table -mt ./american-gut-data/meta_table -at ./american-gut-data/attr_table
```
* Open the browser and enter http://locahost:8087 
* The result of kLDM on American Gut project has been shown at http://diaglecture.com:8087

![](https://github.com/tinglab/kLDM/blob/master/pictures/hierarchical_server.png)

## Contact
If you have any question, please send email to **Yuqing Yang** (yyq410@163.com).

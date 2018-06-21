# kLDM [![language](https://img.shields.io/badge/language-c++-brightgreen.svg)]() [![version](https://img.shields.io/badge/version-v1.0-blue.svg)]()
kLDM is designed to infer multiple association networks based on variation of environmental factors. 

![](https://github.com/tinglab/kLDM/blob/master/pictures/model_explain.jpg)

## Pre-requisites
- C++ 11
- Openmp
- make
- lapack, blas library

## Installation
Just compile source and the executable file 'run' will be generated.
```
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
OTU_Index_k | indexes of OTUs belongs to the cluster k in original 'otu_table'. **:pill: Not all OTUs will be included in the cluster k**
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

# kLDM
kLDM is designed to infer multiple association networks based on variation of environmental factors. [![language](https://img.shields.io/badge/language-c++-brightgreen.svg)]() [![version](https://img.shields.io/badge/version-v1.0-blue.svg)]()

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



## Datasets

Name | Directory | Description | Link
-----|-----------|------------ | ----
Colorectal cancer gut dataset | ./datasets/crc-data/ | 16s rRNA gene sequencing, 490 samples with 117 OTUs and 5 EFs | [Baxter, et al. (2016)](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-016-0290-3) 
TARA Oceans eukaryotic dataset | ./datasets/tara-data/ | 18s rRNA gene sequencing, 221 samples with 67 OTUs and 17 EFs | [Limamendez, et al., 2015](http://science.sciencemag.org/content/348/6237/1262073)
American Gut project dataset | ./datasets/american-gut-data | 16s rRNA gene sequencing, 11946 samples with 216 OTUs and 22 EFs related to human lifestyle | ftp://ftp.microbio.me/AumericanGut, [Official Website](http://americangut.org/)

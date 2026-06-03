## scEntropyR
# Description
This is a tool for estimating cellular entropy through single-cell data analysis.It integrates multiple methods of cell entropy, making the calculation easier and more convenient.
# Installation
Before installing the scgmt package, you need to install the following dependent R packages.
~~~
  Seurat, 
  reshape2, 
  stringr, 
  ggplot2, 
  gridExtra, 
  Biobase, 
  graph,
  BioNet, 
  entropy, 
  cluster, 
  grid, 
  princurve, 
  splines, 
  mgcv,
  lmtest, 
  igraph, 
  tidyverse, 
  CytoTRACE, 
  SCENT, 
  AnnotationDbi,
  org.Hs.eg.db, 
  clusterProfiler, 
  qlcMatrix, 
  Matrix
~~~
Install scEntropyR package:
~~~
devtools::install_github("ZhaoLabs-SJTU/scEntropyR")
~~~
# Usage
Here, we use the pbmc3k dataset that comes with SeuratData.

# Prepare Data
~~~
library(SeuratData)
InstallData("pbmc3k")
data("pbmc3k")
rds <- pbmc3k.final
~~~

# Main calculation process
~~~
library(Seurat)
library(scEntropyR)
rds <- run_entropy(rds,method="slice",species="Human")
rds <- run_entropy(rds,method="cytotrace")
rds <- run_entropy(rds,method="stemid")
rds <- run_entropy(rds,method="scent_ccat",species="Human")
rds <- run_entropy(rds,method="scent_sr",species="Human")


# The final cell entropy score will be stored in the meta.data of the Seurat object
colnames(rds@meta.data)
~~~

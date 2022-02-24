
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scROSHI

<!-- badges: start -->

<!-- badges: end -->

scROSHI identifies cell types based on expression profiles of single
cell analysis by utilizing previously obtained cell type specific gene
sets. It takes into account the hierarchical nature of cell type
relationship and does not require training or annotated data.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ETH-NEXUS/scROSHI")
```

## Example

This is a basic example for the scROSHI function

scROSHI requires three input objects:

*sce\_data*

A SingleCellExperiment object containing the expression profiles of the
single cell analysis

*celltype\_lists*

Marker gene list for all cell types. It can be provided as a list of
genes with cell types as names or as a path to a file containing the
marker genes. Supported file formats are .gmt or .gmx files.

*type\_config*

Config file to define major cell types and hierarchical subtypes. It
should be provided as a two-column data.frame where the first column are
the major cell types and the second column are the subtypes. If several
subtypes exists they should be separated by comma.

``` r
library(scROSHI)
#results <- scROSHI(sce_data = sce,
#                   celltype_lists = path_gmx,
#                   type_config = celltype_config)
```

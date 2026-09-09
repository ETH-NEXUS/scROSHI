<!-- README.md is generated from README.Rmd. Please edit that file -->



# scROSHI

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/scROSHI)](https://cran.r-project.org/package=scROSHI)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

scROSHI identifies cell types based on expression profiles of single cell analysis by utilizing previously obtained cell type specific gene sets. It takes into account the hierarchical nature of cell type relationship and does not require training or annotated data.
A detailed description of the method can be found at:
Michael Prummer, Anne Bertolini, Lars Bosshard, Florian Barkmann, Josephine Yates, Valentina Boeva, The Tumor Profiler Consortium , Daniel Stekhoven, Franziska Singer, scROSHI: robust supervised hierarchical identification of single cells, NAR Genomics and Bioinformatics, Volume 5, Issue 2, June 2023, lqad058, https://doi.org/10.1093/nargab/lqad058

## Installation

From CRAN: [scROSHI](https://cran.r-project.org/package=scROSHI)

```r
install.packages("scROSHI")
```

You can install the development version from [GitHub](https://github.com/ETH-NEXUS/scROSHI) (required R version >= 3.6) with:

``` r
# install.packages("devtools")
devtools::install_github("ETH-NEXUS/scROSHI")
```

## Example

This is a basic example for the scROSHI function

scROSHI requires three input objects:

*sce_data*

A `SingleCellExperiment` object containing the expression profiles of the single cell analysis.

- `dimnames` need to be specified and rownames need to match the gene names in `celltype_list`.
- A column named `barcodes` in `colData` of the SCE object is required.

*celltype_lists*

Marker gene list for all cell types. It can be provided as a list of genes with cell types as names or as a path to a file containing the marker genes. Supported file formats are .gmt or .gmx files.

*type_config*

Config file to define major cell types and hierarchical subtypes. It should be provided as a two-column data.frame where the first column are the major cell types and the second column are the subtypes. If several subtypes exists they should be separated by comma.


``` r
library(scROSHI)
data("test_sce_data")
data("config")
data("marker_list")

results <- scROSHI(sce_data = test_sce_data,
                  celltype_lists = marker_list,
                  type_config = config)
table(results$celltype_final)
#> 
#>                      B.cells                B.cells.naive 
#>                            2                          183 
#>            B.cells.precursor              Dendritic.cells 
#>                           43                           37 
#>                    Monocytes                     NK.cells 
#>                          237                          231 
#>                 Plasma.cells Plasmacytoid.dendritic.cells 
#>                           11                           10 
#>                      T.cells                  T.cells.CD4 
#>                           65                          398 
#>                  T.cells.CD8                    uncertain 
#>                           85                           14
```

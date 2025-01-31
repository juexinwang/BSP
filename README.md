# Big-Small Patch (BSP)

[![DOI](https://zenodo.org/badge/603900314.svg)](https://zenodo.org/doi/10.5281/zenodo.10019315)

Big-small patch is a granularity-guided, data-driven, and parameter-free model for identifying spatial variable genes in 2D and 3D high-throughput spatial transcriptomics data.

![BSP](flowchart.png)

## Update: please use scBSP as a better implementation of BSP algorithm

This package utilizes a granularity-based dimension-agnostic tool, single-cell big-small patch (scBSP), implementing sparse matrix operation and KD-tree/balltree method for distance calculation, for the identification of spatially variable genes on
large-scale data. A corresponding Python library is available at [https://pypi.org/project/scbsp](https://pypi.org/project/scbsp/).

# Installation
This package can be installed on R CRAN
```
# Install sparseMatrixStats if not already installed

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

if (!requireNamespace("sparseMatrixStats", quietly = TRUE)) {
    BiocManager::install("sparseMatrixStats")
}

# Install scBSP from CRAN
install.packages("scBSP")

```

# Tutorial
A detailed tutorial is available at [here](https://castleli.github.io/scBSP/scBSP.html)

# Reference
Li, Jinpu, Yiqing Wang, Mauminah Azam Raina, Chunhui Xu, Li Su, Qi Guo, Qin Ma, Juexin Wang, and Dong Xu. "scBSP: A fast and accurate tool for identifying spatially variable genes from spatial transcriptomic data." bioRxiv (2024).

Wang, Juexin, Jinpu Li, Skyler T. Kramer, Li Su, Yuzhou Chang, Chunhui Xu, Michael T. Eadon, Krzysztof Kiryluk, Qin Ma, and Dong Xu. "Dimension-agnostic and granularity-based spatially variable gene identification using BSP." Nature Communications 14, no. 1 (2023): 7367.


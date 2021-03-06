
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Modelling eye-level visibility of urban green space: Optimising city-wide point-based viewshed computations through prototyping

<!-- badges: start -->

[![DOI:10.5194/agile-giss-3-27-2022](https://zenodo.org/badge/DOI/10.5194/agile-giss-3-27-2022.svg)](https://doi.org/10.5194/agile-giss-3-27-2022)
[![Sample
Data](https://badgen.net/badge/Sample%20Data/10.5281%252Fzenodo.6421423/blue?)](https://doi.org/10.5281/zenodo.6421423)
<!-- badges: end -->

This is the github repository for our publication “Modelling eye-level
visibility of urban green space: Optimising city-wide point-based
viewshed computations through prototyping”. Here we provide reproducible
workflows for our analyses. Data supporting this publication is
accessible via the following DOI
([10.5281/zenodo.6421423](https://doi.org/10.5281/zenodo.6421423)).

## Installation

You can install protoVS from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("STBrinkmann/protoVS")
```

## Workflows

#### 1. Data

The elevation data required downloading and pre-processing. In the
<a href="docs/workflows/00_Data/README.md">data workflow</a> we have
provided R scripts that describe the complete process. Processed input
data has also been uploaded on Zenodo for easy access.

#### 2. Evaluation

Reproducible workflows have been provided for the sections
<a href="docs/workflows/04_1_Experiment_Viewshed/README.md">4.1
Viewshed</a> and
<a href="docs/workflows/04_2_Experiment_VGVI/README.md">4.2 Greenness
visibility</a>.

## Citation

``` r
citation("protoVS")
#> 
#> To cite protoVS in publications use:
#> 
#>   Brinkmann, S. T., Kremer, D., and Walker, B. B.: Modelling eye-level
#>   visibility of urban green space: Optimising city-wide point-based
#>   viewshed computations through prototyping, AGILE GIScience Ser., 3,
#>   27, https://doi.org/10.5194/agile-giss-3-27-2022, 2022.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {Modelling eye-level visibility of urban green space: Optimising city-wide point-based viewshed computations through prototyping},
#>     author = {Sebastian T. Brinkmann and Dominik Kremer and Blake Byron Walker},
#>     journal = {AGILE: GIScience Series},
#>     year = {2022},
#>     volume = {3},
#>     number = {27},
#>     url = {https://doi.org/10.5194/agile-giss-3-27-2022},
#>   }
```

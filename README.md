# fastLDSC

## Fast Matrix-Based LD Score Regression (fastLDSC)

**fastLDSC** relies on GWAS summary statistics to efficiently estimate heritability, genetic covariance, and environmental covariance. Unlike the original command-line implementation, fastLDSC utilizes optimized matrix operations and a C++ backend to significantly accelerate the block-jackknife estimation process. It supports both symmetric analysis (genetic correlation matrices) and rectangular modes (cross-trait analysis), enabling scalable multi-trait analysis for large-scale biobank data.

## Installation

To install the latest version of the fastLDSC package from GitHub, run the following code in R:

```R
library(devtools)
install_github("borangao/fastLDSC")

## Quick Start

See [Tutorial](https://borangao.github.io/meSuSie_Analysis/) for detailed documentation and examples.

## Issues
All feedback, bug reports and suggestions are warmly welcomed! Please make sure to raise issues with a detailed and reproducible exmple and also please provide the output of your sessionInfo() in R! 

How to cite `fastLDSC`
-------------------
Kirtikanth Kalapatapu, Lulu Shang#, Boran Gao#. fastLDSC: Scalable estimation of high-dimensional genetic and environmental covariance using LD score regression.

## Contact
Please contact lshang@mdanderson.org or gao824@purdue.edu

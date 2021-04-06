# Biased Pre-Processing Strategies for Imbalanced Spatio-Temporal Forecasting (Extension)

This repository contains the research compendium of the journal article:

Oliveira, M., Moniz, N., Torgo, L., & Costa, V. S. Biased Resampling Strategies for Imbalanced Spatio-Temporal Forecasting. International Journal of Data Science and Analytics. Accepted.

The journal article is an extension to the conference paper:

Oliveira, M., Moniz, N., Torgo, L., & Costa, V. S. (2019, October). Biased Resampling Strategies for Imbalanced Spatio-Temporal Forecasting. In 2019 IEEE International Conference on Data Science and Advanced Analytics (DSAA) (pp. 100-109). IEEE. doi: [10.1109/DSAA.2019.00024](https://doi.org/10.1109/dsaa.2019.00024) 

You are free to use and/or adapt the code we freely provide. However, we do require that if you do that you cite the paper where these results and code were published.

If you adapt the code to your own needs, you are also required to maintain information on your code concerning the original source of the code (e.g. the URL of this repository) and a reference to the original paper.

## Prerequisites

To install this package, run:

```
library(devtools)  # You need to install this package!
install_github("mrfoliveira/STResampling-JDSA2020",ref="master")
```

If this previous install_github calls somehow fail (there are reports of problems with different libcurl library version on Linux hosts) you may try in alternative the following in R:

```
library(devtools)
install_git("mrfoliveira/STResampling-JDSA2020",ref="master")
```

## Reproducing experiments

To obtain all results shown in the article, run the following lines from the package directory:

```
library(STResamplingJDSA)
PATH <- system.file("inst/", package="STResamplingJDSA")
source(paste0(PATH, "/generate_inds.R"))
source(paste0(PATH, "/exps_internalTuning.R"))
source(paste0(PATH, "/exps_externalPrequential.R"))
```

To generate an HTML report containing all figures and tables in the article, run:

```
library(STResamplingJDSA)
knitr::knit(system.file("inst/report.Rmd", package="STResampling-JDSA2020"))
```

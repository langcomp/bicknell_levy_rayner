Data and analysis for Bicknell, Levy, and Rayner paper
======================================================

## Re-running the analysis
This repository contains all the data (in da1 format) and R analysis code used in the Bicknell, Levy, and Rayner paper. The .Rout and .Rerr files show the R output of each file, so there may be no need to re-run anything (depending on your goals).

If you wish to re-run the analysis, first install all R packages that are used in the various .R files. All required packages are currently available on CRAN except:

* ggintervals, which you can install via devtools::install_github("langcomp/ggintervals")
* lmer.slimr, which you can install via devtools::install_github("kbicknell/lmer.slimr")

After installing all required packages, you should be able to run `make clean && make` to replicate the analysis in a unix/linux/macOS environment. On a Mac, this will require installing the developer tools. If you don't wish to install those, you can simply run each of the R scripts, in the order they're listed in the Makefile. For example:

```
Rscript exclusions.R > exclusions.Rout 2> exclusions.Rerr
Rscript get_saccade_data.R > get_saccade_data.Rout 2> get_saccade_data.Rerr
...
```

Finally, note that the analysis fits quite a lot of big mixed-effects regression models. It is more practical to run it on a big compute server, rather than a laptop.

## Analysis structure

The three experiments were "code-named" Yosemite (expt 1), YosDos (expt 2), and YosDosReplication (expt 3). This was because the experimental shift manipulation caused readers to sometimes overshoot their intended targets, and Yosemite Sam is a cartoon character who shoots too much, i.e., overshoots. Credit goes to Liz Schotter for the codename.

The data for each experiment is in its own subdirectory, named after the experiment's codename.

One final note is that there are two overall analyses reported in the paper: the primary one in the main text and a variant of that analysis making far fewer data exclusions reported in supplementary material. The results in this git repository are for the primary analysis. If you wish to re-run the analysis on the "full" dataset (making fewer exclusions), you need to change `full_dataset <- FALSE` to `full_dataset <- TRUE` in get_data.R.

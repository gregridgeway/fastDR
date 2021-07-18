fastDR
======

R package for fast and easy doubly robust estimation of treatments effects

To install from GitHub, make sure you have installed `remotes`

    install.packages("remotes")

Then install `fastDR` from GitHub:

    remotes::install_github("gregridgeway/fastDR")

`fastDR` requires a modified version of the survey package. The R `survey` package normalizes the weights in a particular way (normalizes across the whole dataset) that can cause numerical instability. `fastDR()` normalizes within the treatment and control groups so that the largest weight is 1 within each. To obtain a patched version of the survey package that avoids the weight rescaling:

    remotes::install_github("gregridgeway/survey","patch-1")

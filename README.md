fastDR
======

R package for fast and easy doubly robust estimation of treatments effects

To install from GitHub, make sure you have installed `remotes`

    install.packages("remotes")

Then install `fastDR` from GitHub:

    remotes::install_github("gregridgeway/fastDR")

`fastDR` requires the `survey` package version 4.0 or later. Earlier `survey` package versions normalized the weights in a particular way (across the whole dataset) that could cause numerical instability. `fastDR()` needs the weights to be normalized within the treatment and control groups so that the largest weight is 1 within each. As of version 4.0, there is a `rescale=FALSE` option for the survey package that resolves this issue.

However, when `rescale=FALSE` survey 4.0 neglects to copy the weights into the dataset sent to `glm()`. You need a patched version of the survey package available with

    remotes::install_github("gregridgeway/survey","patch-2")
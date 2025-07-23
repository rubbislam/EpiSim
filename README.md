# Simulation-based Inference for Epidemiological Dynamics

### Summary

These are notes, scripts, and quarto files from a short course given by [Spencer J. Fox](https://spncrfx.wordpress.com/) and [Qianying Lin](qianylin.com) at [Summer Institute in Statistics and Modeling in Infectious Diseases (SISMID)](https://sph.emory.edu/SISMID/index.html), starting from 2024. All codes and notes are adapted from a [previous version](https://kingaa.github.io/sbied/) developed by [Aaron King](https://kinglab.eeb.lsa.umich.edu/king/) and [Ed Ionides](https://ionides.github.io/). From the [course page](https://rubbislam.quarto.pub/episim/) published based on the current repo, panels lead to topics corresponding to sessions in the course. Within each panel, **R** scripts and the corresponding dataset are included to reproduce results discussed in course. These scripts help students understand, digest, explore, and modify the models, results, and analyses out of their interest and curiosity.

------------------------------------------------------------------------

### Reqauired software

The codes require, at a minimum, [**R**](https://cran.r-project.org/) version 4.2 and [**pomp**](https://kingaa.github.io/pomp/) version 6.3. Windows users must also have the appropriate version of [**Rtools**](https://cran.r-project.org/bin/windows/Rtools/) installed.

Several other packages are required to perform data manipulation, visualization, and analyses:

-   tidyverse
-   ggplot2
-   foreach
-   doParallel
-   doRNG

------------------------------------------------------------------------

### License

![CC-BY-NC](https://i.creativecommons.org/l/by-nc/4.0/88x31.png)

This work is licensed under the [Creative Commons Attribution-NonCommercial 4.0 International License (CC BY-NC 4.0)](https://creativecommons.org/licenses/by-nc/4.0/). Under its terms, you are free to:

-   Share — copy and redistribute the material in any medium or format
-   Adapt — remix, transform, and build upon the material

under the following terms:

-   Attribution — You must give appropriate credit, provide a link to the license, and indicate if changes were made. You may do so in any reasonable manner, but not in any way that suggests the licensor endorses you or your use.
-   NonCommercial — You may not use the material for commercial purposes.
-   No additional restrictions — You may not apply legal terms or technological measures that legally restrict others from doing anything the license permits.

The licensor cannot revoke these freedoms as long as you follow the license terms.

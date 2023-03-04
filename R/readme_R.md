# `R` folder

This folder contains all the scripts to reproduce the results reported in the paper.

The folder is structured as follows:

    .
    ├── dmhs-scripts   # contains the files for the dm-int method
    ├── ma-scripts     # contains the files for the method by Ma et al.
    ├── sensitivity    # contains the files for the simulation study (see `output/sensitivity/numbering.txt` for details on the numbering)
    └── clus.R
    └── genscen.R
    └── misspec.R
    └── parallel_tt_ppmx.R
    └── prior_nclu.R
    └── readme_R.md

-   file `clus.R` reproduce the simulation studies for scenarios S1, S2, and S3 reported in Table 5 in the Supplementary Material;
-   file `genscen.R` generates all the scenarios for the simulation studies reported in Table 1 of the main paper and Table 4 of the Supplementary Material;
-   file `misspec.R` reproduce the simulation studies for scenarios S4, and S5 reported in Table 6 in the Supplementary Material;
-   file `parallel_tt_ppmx.R` reproduce the results of t-ppmx on all Scenarios 1a-3b reported in the main paper (Table 1) and S1a-S3b reported in the Supplementary Material (Table 4), according to the loaded input data;
-   file `prior_nclu.R` reproduce Figure 1 (saved in `figs/plot_pnc.pdf`) and Table 1 of the Supplementary Material.

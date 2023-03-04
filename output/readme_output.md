# `output` folder

This folder contains all the results reported in Table 1 and 2 in the paper, as well as those reported in Table 1, 2, and, 4 in the Supplementary Material.

The folder is structured as follows:

    .
    ├── lgg-da              # contains the output the LGG data analysis
    ├── sensitivity         # contains the output the sensitivity analysis
    ├── simulation-study    # contains the output the simulation studies
    |   ├── main            # simulation studies reported in the main paper
    |   └── supp            # simulation studies reported in the Supplementary Material
    └── readme_output.md

## `lgg-da` folder

It contains outputs regarding the case study.

-   `clu_lgg.RData` contains the metrics evaluating the clustering produced by t-ppmx. It can be obtained running the script `LGG-data-analysis/lgg10fDA_dd.R`
-   `lgg_rep4clunks.RData` contains the output of the t-ppmx function on the whole dataset
-   `lgg10fDA_feb9.RData` contains the output of the t-ppmx function
-   `lggfixedclunks.RData` contains the output of the t-ppmx function given considering the clustering as fixed
-   `pre_meas_lgg.RData` contains the metrics evaluating t-ppmx's prediction perfromance reported in Table 2
-   `resPPMX_lgg.RData` contains the metrics evaluating t-ppmx's goodness of fit

## `sensitivity` folder

It contains the outputs of the sensitivity analysis. Files `res1.RData`-`res9.RData` stores the metrics evaluating the prediction. These results are displayed in Table 1 of the Supplementary Material. Files `clu1.RData`-`clu9.RData` stores the metrics evaluating the clustering. These results are displayed in Table 2 of the Supplementary Material. The file `numbering.txt` explains the relationship between file numbering and parameter elicitation, as discussed in Section A of the Supplementary Material.

## `simulation-study` folder

It contains the output of the simulations study reported in the paper (`main` folder) and in the Supplementary Material B2 (`supp` folder). Each file name in this folder consists of three parts separated by an underscore. The first part of the file name corresponds to the scenario being evaluated, while the second part indicates the method being employed. The third part of the file name indicates whether the file contains metrics related to clustering or prediction. For example, the file named `scen1a_ppmx_res.RData` in the `simulation-study/main` folder contains metrics that evaluate the predictive performance of t-ppmx on Scenario 1a, as reported in the main paper.

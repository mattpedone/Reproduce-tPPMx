# `LGG-data-analysis` folder

This repository contains `R` scripts for analyzing LGG data. The main scripts to reproduce results in Table 2 are:

-   `lgg_samplerJF.R` that reproduce results for dm-int;
-   `lgg10fDA_dd.R` that reproduce results for t-ppmx;
-   `ma_lgg_10f.R` that reproduce results for pam-bp, km-bp, hc-bp, according to the argument `clusterAlg` of the `ConsensusClusterPlus` function at lines 211 and 266.

Input data for these scripts are stored in the `data` folder and the outputs are saved in the `output/lgg-da` folder. Additionally, the script `lgg10fDA_dd.R` generates Figures 2 and 3 in the Supplementary Material, as described in the `figs/readme_figs.md` file.

The file `fixed_clustering.R` takes as input the output of `lgg10fDA_dd.R` and generates Figures 1, 2, and 3 in the paper and Figure 4 in the Supplementary Material. 

The file `PSRF.R` compute the potential scale reduction factor and generates Table 7 in the Supplementary Material. 

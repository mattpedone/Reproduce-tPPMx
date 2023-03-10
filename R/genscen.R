rm(list=ls())
library(treatppmx)
set.seed(121)
### ----
# This script is used to generate all the scenarios for the main simulation study
###

# Simulation Study - paper
# Scenario 1a
scen1a <- treatppmx::genmech_het(npred = 25, nset = 50, overlap = 1.0)
save(scen1a, file = "data/scen1a.RData")

# Scenario 1b
scen1b <- treatppmx::genmech_het(npred = 50, nset = 50, overlap = 1.0)
save(scen1b, file = "data/scen1b.RData")

# Scenario 2a
scen2a <- treatppmx::genmech_het(npred = 25, nset = 50, overlap = .9)
save(scen2a, file = "data/scen2a.RData")

# Scenario 2b
scen2b <- treatppmx::genmech_het(npred = 50, nset = 50, overlap = .9)
save(scen2b, file = "data/scen2b.RData")

# Scenario 3a
scen3a <- treatppmx::genmech_het(npred = 25, nset = 50, overlap = .8)
save(scen3a, file = "data/scen3a.RData")

# Scenario 3b
scen3b <- treatppmx::genmech_het(npred = 50, nset = 50, overlap = .8)
save(scen3b, file = "data/scen3b.RData")

# Simulation Study - supplementary (linear)
# Scenario 1aL
scen1aL <- treatppmx::genmech_hetL(npred = 25, nset = 50, overlap = 1.0)
save(scen1aL, file = "data/scen1aL.RData")

# Scenario 1bL
scen1bL <- treatppmx::genmech_hetL(npred = 50, nset = 50, overlap = 1.0)
save(scen1bL, file = "data/scen1bL.RData")

# Scenario 2aL
scen2aL <- treatppmx::genmech_hetL(npred = 25, nset = 50, overlap = .9)
save(scen2aL, file = "data/scen2aL.RData")

# Scenario 2bL
scen2bL <- treatppmx::genmech_hetL(npred = 50, nset = 50, overlap = .9)
save(scen2bL, file = "data/scen2bL.RData")

# Scenario 3aL
scen3aL <- treatppmx::genmech_hetL(npred = 25, nset = 50, overlap = .8)
save(scen3aL, file = "data/scen3aL.RData")

# Scenario 3bL
scen3bL <- treatppmx::genmech_hetL(npred = 50, nset = 50, overlap = .8)
save(scen3bL, file = "data/scen3bL.RData")

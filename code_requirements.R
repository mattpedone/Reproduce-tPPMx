if (!require("renv"))
  install.packages("renv")
dependencies <- renv::dependencies(path = "../Reproduce-tPPMx")
unique(dependencies$Package)

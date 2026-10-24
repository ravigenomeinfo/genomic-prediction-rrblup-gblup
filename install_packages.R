# Installs the R packages needed by this project (only those that are missing).
pkgs <- c("rrBLUP", "ggplot2", "scales", "testthat")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
message("All required packages are installed.")

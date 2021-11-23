# https://statsandr.com/blog/an-efficient-way-to-install-and-load-r-packages/

# Package names
packages <- c("doRNG", "foreach", "doParallel", "doParallel", "properties", "tidyr", "dplyr", "plyr", 
              "reshape2", "magrittr", "ggplot2", "pomp")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages], repos = "https://cloud.r-project.org/")
}
################################################
## Install packages for the analysis 
## Sylvia Ranjeva, 2019
################################################
#!/usr/bin/Rscript

## Create the personal library if it doesn't exist. Ignore a warning if the directory already exists.
dir.create(Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)


## Install multiple packages
install.packages(c("plyr",
                   "dplyr",
                   "tidyr",
                   "reshape2",
                   "stringr",
                   "ggplot2",
                   "cowplot",
                   "corrplot",
                   "viridis",
                   "polycor",
                   "gridExtra",
                   "devtools",
                   "magrittr"),
                 dependencies = T,
                 Sys.getenv("R_LIBS_USER"), 
                 repos = 'http://cran.us.r-project.org'
)

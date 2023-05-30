library(devtools)
library(roxygen2)

#setwd()

# create the folder for our packages
devtools::create("geneObjects")

# add the .R files containing the functions and objects for our package
## "gene.R", "constructors.R", "accessors.R", "class_specific.R"
# as well as the testthat .R file and folder, and finally the vignette

# produce the documentation via roxygen2 
current.node = as.package("geneObjects")
load_all(current.node$path)
document(current.node)
build_vignettes("geneObjects")

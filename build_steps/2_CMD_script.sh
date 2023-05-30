## to build the .tar.gz file

R CMD build geneObjects

## check the version name

R CMD check geneObjects_1.0.0.tar.gz

## check the errors/warnings/notes, eventually fix them

R CMD INSTALL geneObjects_1.0.0.tar.gz
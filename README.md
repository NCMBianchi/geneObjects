# geneObjects
A test GitHub repository for an R package (as part of the Scientific Programming exam).

## R package content
### Classes
This package provides a set of **S4 classes** for different gene types (i.e. protein-coding genes, long non-coding RNA genes and microRNA genes were used as example, but more could be implemented in the future) which are inherited from a virtual '*Gene*' class that represent genes in general.
The three subclasses inherit several information slots from the '*Gene*' class:
| slot            | type         |
|-----------------|--------------|
| ID              | character    |
| HUGO symbol     | character    |
| gene name       | character    |
| description     | character    |
| gene structure  | GRanges      |
Each subclass (i.e. *ProteinCodingGene*, *LncRNAGene* and *MicroRNAGene* classes) instead has two unique information slots specific for it:
| subclass            | slot              | type         |
|---------------------|-------------------|--------------|
| ProteinCodingGene   | protein ID        | character    |
| ProteinCodingGene   | protein sequence  | character    |
| LncRNAGene          | lncRNA ID         | character    |
| LncRNAGene          | RNA sequence      | character    |
| MicroRNAGene        | microRNA ID       | character    |
| MicroRNAGene        | microRNA sequence | character    |
Examples of further subclasses would be *siRNA*, *snRNA*, *sncRNA* or *matureTranscript*.
### Methods: constructors and accessors
This package also provides methods that act as **constructors** (i.e. *ProteinCodingGene()*, *LncRNAGene()* and *MicroRNAGene()* methods) for each of the three subclasses, based on each combination of shared slots and unique ones – as well as **accessors** (i.e. *getID()* and *setID()*) that allow for easy retrieval of the gene 'ID' for an existing object of either one of the subclasses, and for easy replacement of said gene 'ID', without the user having to access the variable directly (e.g. *gene@symbol*).
Methods for class *Gene* where avoided as that virtual class is only created to act as a root for subclass development.
Moreover, for the sake of simplicity of this package, the only accessors that were implemented were for the 'ID' slot - yet many more could be implemented: the most meaningful ones would be 'HUGO symbol' and 'gene name' from the *Gene* class, while each of the subclass specific IDs (i.e. *protein ID*, *lncRNA ID* and *microRNA ID*) would also make sense as targets for accessor methods to develop in the future.
### Functions
So far only a simple *lengthProduct()* function was implemented in this package: given a certain subclass object, it returns the length of the character string in the product sequence slot (i.e. 'protein sequence', 'RNA sequence' and 'microRNA sequence').
Examples of further functions would be:
- an *alignSequence()* function that takes two subclass objects (checking that they share the same type of subclass) and returns an alignment similar to what many other packages do: https://stackoverflow.com/questions/4497747/how-to-perform-basic-multiple-sequence-alignments-in-r
- a *termSearch()* function that takes a list of subclass objects (any subclasses) and a character string, in order to perform a search for a given term in the 'description' slot of each of the objects in the list, and return the name of the objects that contain that given term: this could even be implement using RegEX syntax, to allow for more customised searches.
### Vignettes and manuals
Documentation was created automatically by 'Roxygen2' during the package implementation, while a simple - yet thorough - vignette was manually written to show how to use the methods and functions of this package.
### Test scripts
Scripts redacted according to syntax required by the 'testthat' package and function were prepared and launched, returning to error or warnings. They were also provided, according to repositories' good practices.

## 'build_steps' folder
This extra folder contains three files:
- '1_build.R': the script used to create the package directory and automatically create documentation and manuals (requires 'Roxygen2' and 'devtools')
- '2_CMD_script.sh': the script used to build the .tar.gz compressed file and to test the package
- '3_git.sh': the script used to push the package on the GitHub repository 
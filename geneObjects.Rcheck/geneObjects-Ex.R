pkgname <- "geneObjects"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('geneObjects')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("Gene-class")
### * Gene-class

flush(stderr()); flush(stdout())

### Name: Gene-class
### Title: Gene class
### Aliases: Gene-class

### ** Examples

# Not usually instantiated directly, used as a superclass for other gene classes.




cleanEx()
nameEx("LncRNAGene-class")
### * LncRNAGene-class

flush(stderr()); flush(stdout())

### Name: LncRNAGene-class
### Title: Long non-coding RNA gene class
### Aliases: LncRNAGene-class

### ** Examples

gene_structure <- GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges::IRanges(start = 3000, end = 4000), strand = "-")
LncRNAGene(ID = "ENSG000002", HUGO_symbol = "LINC00001",
           gene_name = "Long intergenic non-protein coding RNA 1",
           description = "Long non-coding RNA",
           gene_structure = gene_structure,
           lncRNA_ID = "ENST000004", RNA_sequence = "ACUGCUAGCUAGUCA...")




cleanEx()
nameEx("LncRNAGene")
### * LncRNAGene

flush(stderr()); flush(stdout())

### Name: LncRNAGene
### Title: Constructor function for LncRNAGene class
### Aliases: LncRNAGene

### ** Examples

gene_structure <- GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges::IRanges(3000, 4000), strand = "-")
LncRNAGene(ID = "ENSG000002", HUGO_symbol = "LINC00001",
           gene_name = "Long intergenic non-protein coding RNA 1",
           description = "Long non-coding RNA",
           gene_structure = gene_structure,
           lncRNA_ID = "ENST000004", RNA_sequence = "ACUGCUAGCUAGUCA...")



cleanEx()
nameEx("MicroRNAGene-class")
### * MicroRNAGene-class

flush(stderr()); flush(stdout())

### Name: MicroRNAGene-class
### Title: MicroRNA gene class
### Aliases: MicroRNAGene-class

### ** Examples

gene_structure <- GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges::IRanges(start = 5000, end = 6000), strand = "+")
MicroRNAGene(ID = "ENSG000003", HUGO_symbol = "MIR00001",
             gene_name = "MicroRNA 1",
             description = "MicroRNA",
             gene_structure = gene_structure,
             microRNA_ID = "ENST000005", microRNA_seed_sequence = "UGAGGUAGUAGGUUGUAUGGUAG")




cleanEx()
nameEx("MicroRNAGene")
### * MicroRNAGene

flush(stderr()); flush(stdout())

### Name: MicroRNAGene
### Title: Constructor function for MicroRNAGene class
### Aliases: MicroRNAGene

### ** Examples

gene_structure <- GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges::IRanges(5000, 6000), strand = "+")
MicroRNAGene(ID = "ENSG000003", HUGO_symbol = "MIR00001",
             gene_name = "MicroRNA 1",
             description = "MicroRNA",
             gene_structure = gene_structure,
             microRNA_ID = "ENST000005", microRNA_seed_sequence = "UGAGGUAGUAGGUUGUAUGGUAG")



cleanEx()
nameEx("ProteinCodingGene-class")
### * ProteinCodingGene-class

flush(stderr()); flush(stdout())

### Name: ProteinCodingGene-class
### Title: Protein-coding gene class
### Aliases: ProteinCodingGene-class

### ** Examples

gene_structure <- GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges::IRanges(start = 1000, end = 2000), strand = "+")
ProteinCodingGene(ID = "ENSG000001", HUGO_symbol = "BRCA1",
                  gene_name = "Breast cancer 1",
                  description = "Breast cancer type 1 susceptibility protein",
                  gene_structure = gene_structure,
                  protein_ID = "ENSP000003", protein_sequence = "MENSDRNSIKVAV...")




cleanEx()
nameEx("ProteinCodingGene")
### * ProteinCodingGene

flush(stderr()); flush(stdout())

### Name: ProteinCodingGene
### Title: Constructor function for ProteinCodingGene class
### Aliases: ProteinCodingGene

### ** Examples

gene_structure <- GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges::IRanges(1000, 2000), strand = "+")
ProteinCodingGene(ID = "ENSG000001", HUGO_symbol = "BRCA1",
                  gene_name = "Breast cancer 1",
                  description = "Breast cancer type 1 susceptibility protein",
                  gene_structure = gene_structure,
                  protein_ID = "ENSP000003", protein_sequence = "MENSDRNSIKVAV...")



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')

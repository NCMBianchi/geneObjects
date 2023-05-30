## ----style, echo = FALSE, results = 'asis'------------------------------------
library(BiocStyle)

## ---- echo = FALSE------------------------------------------------------------
library(knitr)
library(GenomicRanges)

## ---- echo = FALSE------------------------------------------------------------
library(geneObjects)

## -----------------------------------------------------------------------------
gene_structure <- GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges::IRanges(start = 5000, end = 6000), strand = "+")

## -----------------------------------------------------------------------------
# create a ProteinCodingGene object
pc_gene <- ProteinCodingGene(ID = "ENSG000001", HUGO_symbol = "BRCA1", gene_name = "Breast cancer 1", description = "Breast cancer type 1 susceptibility protein", gene_structure = gene_structure, protein_ID = "ENSP000003", protein_sequence = "MENSDRNSIKVAV")

## -----------------------------------------------------------------------------
# create a LncRNAGene object
lnc_gene <- LncRNAGene(ID = "ENSG000002", HUGO_symbol = "LINC00001", gene_name = "Long intergenic non-protein coding RNA 1", description = "Long non-coding RNA", gene_structure = gene_structure, lncRNA_ID = "ENST000004", RNA_sequence = "ACUGCUAGCUAGUCA")

## -----------------------------------------------------------------------------
# create a MicroRNAGene object
mir_gene <-  MicroRNAGene(ID = "ENSG000003", HUGO_symbol = "MIR00001", gene_name = "MicroRNA 1", description = "MicroRNA", gene_structure = gene_structure, microRNA_ID = "ENST000005", microRNA_seed_sequence = "UGAGGUAGUAGGUUGUAUGGUAG")

## -----------------------------------------------------------------------------
# set ID
setID(pc_gene) <- "Gene1"

# get ID
getID(pc_gene)

## -----------------------------------------------------------------------------
# compute length of gene product
lengthProduct(pc_gene)


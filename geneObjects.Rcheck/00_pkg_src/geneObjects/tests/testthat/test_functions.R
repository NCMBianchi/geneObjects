context("Testing lengthProduct methods")

gene_structure <- GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges::IRanges(start = 1000, end = 2000), strand = "+")

test_that("lengthProduct works correctly for ProteinCodingGene objects", {
  # Prepare an example ProteinCodingGene object for the tests
  test_gene <- ProteinCodingGene(ID = "ENSG000001", HUGO_symbol = "BRCA1", gene_name = "Breast cancer 1", description = "Breast cancer type 1 susceptibility protein", gene_structure = gene_structure, protein_ID = "ENSP000003", protein_sequence = "MENSDRNSI")
  
  # Test lengthProduct
  expect_equal(lengthProduct(test_gene), 9)
})

test_that("lengthProduct works correctly for LncRNAGene objects", {
  # Prepare an example LncRNAGene object for the tests
  test_gene <- LncRNAGene(ID = "ENSG000002", HUGO_symbol = "LINC00001", gene_name = "Long intergenic non-protein coding RNA 1", description = "Long non-coding RNA", gene_structure = gene_structure, lncRNA_ID = "ENST000004", RNA_sequence = "ACUGCUAGCUAGUCA")
  
  # Test lengthProduct
  expect_equal(lengthProduct(test_gene), 15)
})

test_that("lengthProduct works correctly for MicroRNAGene objects", {
  # Prepare an example MicroRNAGene object for the tests
  test_gene <- MicroRNAGene(ID = "ENSG000003", HUGO_symbol = "MIR00001", gene_name = "MicroRNA 1", description = "MicroRNA", gene_structure = gene_structure, microRNA_ID = "ENST000005", microRNA_seed_sequence = "UGAGGUAGUAGGUUGUAUGGUAG")
  
  # Test lengthProduct
  expect_equal(lengthProduct(test_gene), 23)
})
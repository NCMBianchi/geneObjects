context("Testing ID methods")

test_that("ID works correctly for Gene object's subclasses", {
  # Prepare example Gene objects for the tests
  gene_structure <- GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges::IRanges(start = 1000, end = 2000), strand = "+")
  test_pcg <- ProteinCodingGene(ID = "pcg1", HUGO_symbol = "BRCA1", gene_name = "Breast cancer 1", description = "Breast cancer type 1 susceptibility protein", gene_structure = gene_structure, protein_ID = "ENSP000003", protein_sequence = "MENSDRNSI")
  test_lng <- LncRNAGene(ID = "lng1", HUGO_symbol = "LINC00001", gene_name = "Long intergenic non-protein coding RNA 1", description = "Long non-coding RNA", gene_structure = gene_structure, lncRNA_ID = "ENST000004", RNA_sequence = "ACUGCUAGCUAGUCA")
  test_mrg <- MicroRNAGene(ID = "mrg1", HUGO_symbol = "MIR00001", gene_name = "MicroRNA 1", description = "MicroRNA", gene_structure = gene_structure, microRNA_ID = "ENST000005", microRNA_seed_sequence = "UGAGGUAGUAGGUUGUAUGGUAG")
  
  # Test getter and setter for ProteinCodingGene object
  expect_equal(getID(test_pcg), "pcg1")
  setID(test_pcg) <- "pcg2"
  expect_equal(getID(test_pcg), "pcg2")
  
  # Test getter and setter for LncRNAGene object
  expect_equal(getID(test_lng), "lng1")
  setID(test_lng) <- "lng2"
  expect_equal(getID(test_lng), "lng2")
  
  # Test getter and setter for MicroRNAGene object
  expect_equal(getID(test_mrg), "mrg1")
  setID(test_mrg) <- "mrg2"
  expect_equal(getID(test_mrg), "mrg2")
})
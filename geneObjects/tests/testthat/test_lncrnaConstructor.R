context("LncRNAGene constructor")

test_that("LncRNAGene constructor works correctly", {
  gene_structure <- GRanges("chr1", IRanges(3000, 4000), strand = "-")
  gene <- LncRNAGene(ID = "ENSG000002", HUGO_symbol = "LINC00001",
                     gene_name = "Long intergenic non-protein coding RNA 1",
                     description = "Long non-coding RNA",
                     gene_structure = gene_structure,
                     lncRNA_ID = "ENST000004", RNA_sequence = "ACUGCUAGCUAGUCA...")
  expect_s4_class(gene, "LncRNAGene")
  expect_equal(gene@ID, "ENSG000002")
  expect_equal(gene@HUGO_symbol, "LINC00001")
  expect_equal(gene@gene_name, "Long intergenic non-protein coding RNA 1")
  expect_equal(gene@description, "Long non-coding RNA")
  expect_equal(gene@gene_structure, gene_structure)
  expect_equal(gene@lncRNA_ID, "ENST000004")
  expect_equal(gene@RNA_sequence, "ACUGCUAGCUAGUCA...")
})

test_that("LncRNAGene constructor throws error with incorrect inputs", {
  gene_structure <- GRanges("chr1", IRanges(3000, 4000), strand = "-")
  
  # testing with incorrect gene_structure input
  expect_error(LncRNAGene(ID = "ENSG000002", HUGO_symbol = "LINC00001",
                          gene_name = "Long intergenic non-protein coding RNA 1",
                          description = "Long non-coding RNA",
                          gene_structure = "incorrect_input",
                          lncRNA_ID = "ENST000004", RNA_sequence = "ACUGCUAGCUAGUCA..."),
               "'gene_structure' must be of class 'GRanges'")
  
  # testing with missing lncRNA_ID input
  expect_error(LncRNAGene(ID = "ENSG000002", HUGO_symbol = "LINC00001",
                          gene_name = "Long intergenic non-protein coding RNA 1",
                          description = "Long non-coding RNA",
                          gene_structure = gene_structure,
                          RNA_sequence = "ACUGCUAGCUAGUCA..."),
               "'lncRNA_ID' must be of class 'character'")
  
  # testing with incorrect ID input
  expect_error(LncRNAGene(ID = 123, HUGO_symbol = "LINC00001",
                          gene_name = "Long intergenic non-protein coding RNA 1",
                          description = "Long non-coding RNA",
                          gene_structure = gene_structure,
                          lncRNA_ID = "ENST000004", RNA_sequence = "ACUGCUAGCUAGUCA..."),
               "'ID' must be of class 'character'")
  
  # testing with missing HUGO_symbol input
  expect_error(LncRNAGene(ID = "ENSG000002",
                          gene_name = "Long intergenic non-protein coding RNA 1",
                          description = "Long non-coding RNA",
                          gene_structure = gene_structure,
                          lncRNA_ID = "ENST000004", RNA_sequence = "ACUGCUAGCUAGUCA..."),
               "'HUGO_symbol' must be of class 'character'")
  
  # testing with missing gene_name input
  expect_error(LncRNAGene(ID = "ENSG000002", HUGO_symbol = "LINC00001",
                          description = "Long non-coding RNA",
                          gene_structure = gene_structure,
                          lncRNA_ID = "ENST000004", RNA_sequence = "ACUGCUAGCUAGUCA..."),
               "'gene_name' must be of class 'character'")
  
  # testing with missing description input
  expect_error(LncRNAGene(ID = "ENSG000002", HUGO_symbol = "LINC00001",
                          gene_name = "Long intergenic non-protein coding RNA 1",
                          gene_structure = gene_structure,
                          lncRNA_ID = "ENST000004", RNA_sequence = "ACUGCUAGCUAGUCA..."),
               "'description' must be of class 'character'")
  
  # testing with missing RNA_sequence input
  expect_error(LncRNAGene(ID = "ENSG000002", HUGO_symbol = "LINC00001",
                          gene_name = "Long intergenic non-protein coding RNA 1",
                          description = "Long non-coding RNA",
                          gene_structure = gene_structure,
                          lncRNA_ID = "ENST000004"),
               "'RNA_sequence' must be of class 'character'")
})
context("MicroRNAGene constructor")

test_that("MicroRNAGene constructor works correctly", {
  gene_structure <- GRanges("chr1", IRanges(5000, 6000), strand = "+")
  gene <- MicroRNAGene(ID = "ENSG000003", HUGO_symbol = "MIR00001",
                       gene_name = "MicroRNA 1",
                       description = "MicroRNA",
                       gene_structure = gene_structure,
                       microRNA_ID = "ENST000005", microRNA_seed_sequence = "UGAGGUAGUAGGUUGUAUGGUAG")
  expect_s4_class(gene, "MicroRNAGene")
  expect_equal(gene@ID, "ENSG000003")
  expect_equal(gene@HUGO_symbol, "MIR00001")
  expect_equal(gene@gene_name, "MicroRNA 1")
  expect_equal(gene@description, "MicroRNA")
  expect_equal(gene@gene_structure, gene_structure)
  expect_equal(gene@microRNA_ID, "ENST000005")
  expect_equal(gene@microRNA_seed_sequence, "UGAGGUAGUAGGUUGUAUGGUAG")
})

test_that("MicroRNAGene constructor throws error with incorrect inputs", {
  gene_structure <- GRanges("chr1", IRanges(5000, 6000), strand = "+")
  
  # testing with incorrect gene_structure input
  expect_error(MicroRNAGene(ID = "ENSG000003", HUGO_symbol = "MIR00001",
                            gene_name = "MicroRNA 1",
                            description = "MicroRNA",
                            gene_structure = "incorrect_input",
                            microRNA_ID = "ENST000005", microRNA_seed_sequence = "UGAGGUAGUAGGUUGUAUGGUAG"),
               "'gene_structure' must be of class 'GRanges'")
  
  # testing with missing microRNA_ID input
  expect_error(MicroRNAGene(ID = "ENSG000003", HUGO_symbol = "MIR00001",
                            gene_name = "MicroRNA 1",
                            description = "MicroRNA",
                            gene_structure = gene_structure,
                            microRNA_seed_sequence = "UGAGGUAGUAGGUUGUAUGGUAG"),
               "'microRNA_ID' must be of class 'character'")
  
  # testing with incorrect ID input
  expect_error(MicroRNAGene(ID = 123, HUGO_symbol = "MIR00001",
                            gene_name = "MicroRNA 1",
                            description = "MicroRNA",
                            gene_structure = gene_structure,
                            microRNA_ID = "ENST000005", microRNA_seed_sequence = "UGAGGUAGUAGGUUGUAUGGUAG"),
               "'ID' must be of class 'character'")
  
  # testing with missing HUGO_symbol input
  expect_error(MicroRNAGene(ID = "ENSG000003",
                            gene_name = "MicroRNA 1",
                            description = "MicroRNA",
                            gene_structure = gene_structure,
                            microRNA_ID = "ENST000005", microRNA_seed_sequence = "UGAGGUAGUAGGUUGUAUGGUAG"),
               "'HUGO_symbol' must be of class 'character'")
  
  # testing with missing gene_name input
  expect_error(MicroRNAGene(ID = "ENSG000003", HUGO_symbol = "MIR00001",
                            description = "MicroRNA",
                            gene_structure = gene_structure,
                            microRNA_ID = "ENST000005", microRNA_seed_sequence = "UGAGGUAGUAGGUUGUAUGGUAG"),
               "'gene_name' must be of class 'character'")
  
  # testing with missing description input
  expect_error(MicroRNAGene(ID = "ENSG000003", HUGO_symbol = "MIR00001",
                            gene_name = "MicroRNA 1",
                            gene_structure = gene_structure,
                            microRNA_ID = "ENST000005", microRNA_seed_sequence = "UGAGGUAGUAGGUUGUAUGGUAG"),
               "'description' must be of class 'character'")
  
  # testing with missing microRNA_seed_sequence input
  expect_error(MicroRNAGene(ID = "ENSG000003", HUGO_symbol = "MIR00001",
                            gene_name = "MicroRNA 1",
                            description = "MicroRNA",
                            gene_structure = gene_structure,
                            microRNA_ID = "ENST000005"),
               "'microRNA_seed_sequence' must be of class 'character'")
})
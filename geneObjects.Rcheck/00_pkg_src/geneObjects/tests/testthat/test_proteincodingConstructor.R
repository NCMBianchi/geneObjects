context("ProteinCodingGene constructor")

test_that("ProteinCodingGene constructor works correctly", {
  gene_structure <- GRanges("chr1", IRanges(1000, 2000), strand = "+")
  gene <- ProteinCodingGene(ID = "ENSG000001", HUGO_symbol = "BRCA1",
                            gene_name = "Breast cancer 1",
                            description = "Breast cancer type 1 susceptibility protein",
                            gene_structure = gene_structure,
                            protein_ID = "ENSP000003", protein_sequence = "MENSDRNSIKVAV...")
  expect_s4_class(gene, "ProteinCodingGene")
  expect_equal(gene@ID, "ENSG000001")
  expect_equal(gene@HUGO_symbol, "BRCA1")
  expect_equal(gene@gene_name, "Breast cancer 1")
  expect_equal(gene@description, "Breast cancer type 1 susceptibility protein")
  expect_equal(gene@gene_structure, gene_structure)
  expect_equal(gene@protein_ID, "ENSP000003")
  expect_equal(gene@protein_sequence, "MENSDRNSIKVAV...")
})

test_that("ProteinCodingGene constructor throws error with incorrect inputs", {
  gene_structure <- GRanges("chr1", IRanges(1000, 2000), strand = "+")
  
  # testing with incorrect gene_structure input
  expect_error(ProteinCodingGene(ID = "ENSG000001", HUGO_symbol = "BRCA1",
                                 gene_name = "Breast cancer 1",
                                 description = "Breast cancer type 1 susceptibility protein",
                                 gene_structure = "incorrect_input",
                                 protein_ID = "ENSP000003", protein_sequence = "MENSDRNSIKVAV..."),
               "'gene_structure' must be of class 'GRanges'")
  
  # testing with missing protein_ID input
  expect_error(ProteinCodingGene(ID = "ENSG000001", HUGO_symbol = "BRCA1",
                                 gene_name = "Breast cancer 1",
                                 description = "Breast cancer type 1 susceptibility protein",
                                 gene_structure = gene_structure,
                                 protein_sequence = "MENSDRNSIKVAV..."),
               "'protein_ID' must be of class 'character'")
  
  # testing with incorrect ID input
  expect_error(ProteinCodingGene(ID = 123, HUGO_symbol = "BRCA1",
                                 gene_name = "Breast cancer 1",
                                 description = "Breast cancer type 1 susceptibility protein",
                                 gene_structure = gene_structure,
                                 protein_ID = "ENSP000003", protein_sequence = "MENSDRNSIKVAV..."),
               "'ID' must be of class 'character'")
  
  # testing with missing HUGO_symbol input
  expect_error(ProteinCodingGene(ID = "ENSG000001",
                                 gene_name = "Breast cancer 1",
                                 description = "Breast cancer type 1 susceptibility protein",
                                 gene_structure = gene_structure,
                                 protein_ID = "ENSP000003", protein_sequence = "MENSDRNSIKVAV..."),
               "'HUGO_symbol' must be of class 'character'")
  
  # testing with missing gene_name input
  expect_error(ProteinCodingGene(ID = "ENSG000001", HUGO_symbol = "BRCA1",
                                 description = "Breast cancer type 1 susceptibility protein",
                                 gene_structure = gene_structure,
                                 protein_ID = "ENSP000003", protein_sequence = "MENSDRNSIKVAV..."),
               "'gene_name' must be of class 'character'")
  
  # testing with missing description input
  expect_error(ProteinCodingGene(ID = "ENSG000001", HUGO_symbol = "BRCA1",
                                 gene_name = "Breast cancer 1",
                                 gene_structure = gene_structure,
                                 protein_ID = "ENSP000003", protein_sequence = "MENSDRNSIKVAV..."),
               "'description' must be of class 'character'")
  
  # testing with missing protein_sequence input
  expect_error(ProteinCodingGene(ID = "ENSG000001", HUGO_symbol = "BRCA1",
                                 gene_name = "Breast cancer 1",
                                 description = "Breast cancer type 1 susceptibility protein",
                                 gene_structure = gene_structure,
                                 protein_ID = "ENSP000003"),
               "'protein_sequence' must be of class 'character'")
})
#' Constructor function for ProteinCodingGene class
#'
#' @description This function creates a new instance of the ProteinCodingGene class.
#' @param ID The gene ID (Ensembl ID or NCBI gene ID).
#' @param HUGO_symbol The HUGO gene symbol.
#' @param gene_name The full gene name.
#' @param description A brief description of the gene.
#' @param gene_structure The gene structure represented as a GRanges object.
#' @param protein_ID The protein ID.
#' @param protein_sequence The protein sequence.
#' @usage ProteinCodingGene(ID, HUGO_symbol, gene_name, description, gene_structure, protein_ID, protein_sequence)
#' @return An instance of the ProteinCodingGene class.
#' @author Niccolò Bianchi\cr UniMi + PoliMI\cr Maintainer: Niccolò
#' Bianchi\cr E-Mail: <niccolo.bianchi2@@studenti.unimi.it>
#' @references \url{https://www.ensembl.org/index.html?redirect=no}\cr
#' @examples
#' gene_structure <- GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges::IRanges(1000, 2000), strand = "+")
#' ProteinCodingGene(ID = "ENSG000001", HUGO_symbol = "BRCA1",
#'                   gene_name = "Breast cancer 1",
#'                   description = "Breast cancer type 1 susceptibility protein",
#'                   gene_structure = gene_structure,
#'                   protein_ID = "ENSP000003", protein_sequence = "MENSDRNSIKVAV...")
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom methods is
#' @importFrom methods new
#' @export
ProteinCodingGene <- function(ID, HUGO_symbol, gene_name, description, gene_structure, protein_ID, protein_sequence){
  if(missing(ID) || !is.character(ID)) stop("'ID' must be of class 'character'")
  if(missing(HUGO_symbol) || !is.character(HUGO_symbol)) stop("'HUGO_symbol' must be of class 'character'")
  if(missing(gene_name) || !is.character(gene_name)) stop("'gene_name' must be of class 'character'")
  if(missing(description) || !is.character(description)) stop("'description' must be of class 'character'")
  if(missing(gene_structure) || !is(gene_structure, "GRanges")) stop("'gene_structure' must be of class 'GRanges'")
  if(missing(protein_ID) || !is.character(protein_ID)) stop("'protein_ID' must be of class 'character'")
  if(missing(protein_sequence) || !is.character(protein_sequence)) stop("'protein_sequence' must be of class 'character'")
  
  new("ProteinCodingGene", 
      ID = ID, 
      HUGO_symbol = HUGO_symbol, 
      gene_name = gene_name, 
      description = description, 
      gene_structure = gene_structure, 
      protein_ID = protein_ID, 
      protein_sequence = protein_sequence
  )
}

#' Constructor function for LncRNAGene class
#'
#' @description This function creates a new instance of the LncRNAGene class.
#' @param ID The gene ID (Ensembl ID or NCBI gene ID).
#' @param HUGO_symbol The HUGO gene symbol.
#' @param gene_name The full gene name.
#' @param description A brief description of the gene.
#' @param gene_structure The gene structure represented as a GRanges object.
#' @param lncRNA_ID The lncRNA ID.
#' @param RNA_sequence The RNA sequence.
#' @usage LncRNAGene(ID, HUGO_symbol, gene_name, description, gene_structure, lncRNA_ID, RNA_sequence)
#' @return An instance of the LncRNAGene class.
#' @author Niccolò Bianchi\cr UniMi + PoliMI\cr Maintainer: Niccolò
#' Bianchi\cr E-Mail: <niccolo.bianchi2@@studenti.unimi.it>
#' @references \url{https://www.ensembl.org/index.html?redirect=no}\cr
#' @examples
#' gene_structure <- GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges::IRanges(3000, 4000), strand = "-")
#' LncRNAGene(ID = "ENSG000002", HUGO_symbol = "LINC00001",
#'            gene_name = "Long intergenic non-protein coding RNA 1",
#'            description = "Long non-coding RNA",
#'            gene_structure = gene_structure,
#'            lncRNA_ID = "ENST000004", RNA_sequence = "ACUGCUAGCUAGUCA...")
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom methods is
#' @importFrom methods new
#' @export
LncRNAGene <- function(ID, HUGO_symbol, gene_name, description,
                       gene_structure, lncRNA_ID, RNA_sequence) {
  if (missing(ID) || !is(ID, "character")) stop("'ID' must be of class 'character'")
  if (missing(HUGO_symbol) || !is(HUGO_symbol, "character")) stop("'HUGO_symbol' must be of class 'character'")
  if (missing(gene_name) || !is(gene_name, "character")) stop("'gene_name' must be of class 'character'")
  if (missing(description) || !is(description, "character")) stop("'description' must be of class 'character'")
  if (missing(gene_structure) || !is(gene_structure, "GRanges")) stop("'gene_structure' must be of class 'GRanges'")
  if (missing(lncRNA_ID) || !is(lncRNA_ID, "character")) stop("'lncRNA_ID' must be of class 'character'")
  if (missing(RNA_sequence) || !is(RNA_sequence, "character")) stop("'RNA_sequence' must be of class 'character'")
  
  new("LncRNAGene", ID = ID, HUGO_symbol = HUGO_symbol, gene_name = gene_name,
      description = description, gene_structure = gene_structure,
      lncRNA_ID = lncRNA_ID, RNA_sequence = RNA_sequence)
}

#' Constructor function for MicroRNAGene class
#'
#' @description This function creates a new instance of the MicroRNAGene class.
#' @param ID The gene ID (Ensembl ID or NCBI gene ID).
#' @param HUGO_symbol The HUGO gene symbol.
#' @param gene_name The full gene name.
#' @param description A brief description of the gene.
#' @param gene_structure The gene structure represented as a GRanges object.
#' @param microRNA_ID The microRNA ID.
#' @param microRNA_seed_sequence The microRNA seed sequence.
#' @usage MicroRNAGene(ID, HUGO_symbol, gene_name, description, gene_structure, microRNA_ID, microRNA_seed_sequence)
#' @return An instance of the MicroRNAGene class.
#' @author Niccolò Bianchi\cr UniMi + PoliMI\cr Maintainer: Niccolò
#' Bianchi\cr E-Mail: <niccolo.bianchi2@@studenti.unimi.it>
#' @references \url{https://www.ensembl.org/index.html?redirect=no}\cr
#' @examples
#' gene_structure <- GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges::IRanges(5000, 6000), strand = "+")
#' MicroRNAGene(ID = "ENSG000003", HUGO_symbol = "MIR00001",
#'              gene_name = "MicroRNA 1",
#'              description = "MicroRNA",
#'              gene_structure = gene_structure,
#'              microRNA_ID = "ENST000005", microRNA_seed_sequence = "UGAGGUAGUAGGUUGUAUGGUAG")
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom methods is
#' @importFrom methods new
#' @export
MicroRNAGene <- function(ID, HUGO_symbol, gene_name, description, gene_structure,
                         microRNA_ID, microRNA_seed_sequence) {
  if (missing(ID) || !is(ID, "character")) stop("'ID' must be of class 'character'")
  if (missing(HUGO_symbol) || !is(HUGO_symbol, "character")) stop("'HUGO_symbol' must be of class 'character'")
  if (missing(gene_name) || !is(gene_name, "character")) stop("'gene_name' must be of class 'character'")
  if (missing(description) || !is(description, "character")) stop("'description' must be of class 'character'")
  if (missing(gene_structure) || !is(gene_structure, "GRanges")) stop("'gene_structure' must be of class 'GRanges'")
  if (missing(microRNA_ID) || !is(microRNA_ID, "character")) stop("'microRNA_ID' must be of class 'character'")
  if (missing(microRNA_seed_sequence) || !is(microRNA_seed_sequence, "character")) stop("'microRNA_seed_sequence' must be of class 'character'")
  
  new("MicroRNAGene", ID = ID, HUGO_symbol = HUGO_symbol, gene_name = gene_name,
      description = description, gene_structure = gene_structure,
      microRNA_ID = microRNA_ID, microRNA_seed_sequence = microRNA_seed_sequence)
}
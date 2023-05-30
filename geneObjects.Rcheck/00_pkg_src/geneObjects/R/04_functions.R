#' Generic function for length of gene product
#' 
#' @description This method gets the length of the product of the gene.
#' @param gene An object of class Gene.
#' @return The length of the gene product.
#' @rdname lengthProduct-methods
#' @export
setGeneric("lengthProduct", function(gene) standardGeneric("lengthProduct"))

#' Length of protein-coding gene product
#' 
#' @description Method to get the length of the product for a ProteinCodingGene object.
#' @param gene An object of class ProteinCodingGene.
#' @return The length of the protein sequence.
#' @rdname lengthProduct-methods
#' @importFrom methods is
#' @export
setMethod("lengthProduct", "ProteinCodingGene", function(gene) {
  if (!is(gene@protein_sequence, "character")) stop("'protein_sequence' must be of class 'character'")
  nchar(gene@protein_sequence)
})

#' Length of long non-coding RNA gene product
#' 
#' @description Method to get the length of the product for a LncRNAGene object.
#' @param gene An object of class LncRNAGene.
#' @return The length of the RNA sequence.
#' @rdname lengthProduct-methods
#' @importFrom methods is
#' @export
setMethod("lengthProduct", "LncRNAGene", function(gene) {
  if (!is(gene@RNA_sequence, "character")) stop("'RNA_sequence' must be of class 'character'")
  nchar(gene@RNA_sequence)
})

#' Length of microRNA gene product
#' 
#' @description Method to get the length of the product for a MicroRNAGene object.
#' @param gene An object of class MicroRNAGene.
#' @return The length of the microRNA seed sequence.
#' @rdname lengthProduct-methods
#' @importFrom methods is
#' @export
setMethod("lengthProduct", "MicroRNAGene", function(gene) {
  if (!is(gene@microRNA_seed_sequence, "character")) stop("'microRNA_seed_sequence' must be of class 'character'")
  nchar(gene@microRNA_seed_sequence)
})
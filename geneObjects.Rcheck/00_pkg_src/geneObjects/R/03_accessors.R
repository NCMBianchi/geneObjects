#' Generic getter for gene ID
#'
#' @description This method gets the ID for the gene.
#' @param gene An object of class Gene.
#' @return A character string representing the ID of the gene.
#' @rdname ID-methods
#' @export
setGeneric("getID", function(gene) standardGeneric("getID"))

#' Gene-specific getter for gene ID
#' 
#' @description Method to get the ID for a Gene object.
#' @param gene An object of class Gene.
#' @return A character string representing the ID of the gene.
#' @rdname ID-methods
#' @importFrom methods is
#' @export
setMethod("getID", "Gene", function(gene) {
  if (!is(gene, "Gene")) stop("'gene' must be of class 'Gene'")
  return(gene@ID)
})

#' ProteinCodingGene-specific getter for gene ID
#' 
#' @description Method to get the ID for a ProteinCodingGene object.
#' @param gene An object of class ProteinCodingGene.
#' @return A character string representing the ID of the gene.
#' @rdname ID-methods
#' @importFrom methods is
#' @export
setMethod("getID", "ProteinCodingGene", function(gene) {
  if (!is(gene, "ProteinCodingGene")) stop("'gene' must be of class 'ProteinCodingGene'")
  return(gene@ID)
})

#' LncRNAGene-specific getter for gene ID
#' 
#' @description Method to get the ID for a LncRNAGene object.
#' @param gene An object of class LncRNAGene.
#' @return A character string representing the ID of the gene.
#' @rdname ID-methods
#' @importFrom methods is
#' @export
setMethod("getID", "LncRNAGene", function(gene) {
  if (!is(gene, "LncRNAGene")) stop("'gene' must be of class 'LncRNAGene'")
  return(gene@ID)
})

#' MicroRNAGene-specific getter for gene ID
#' 
#' @description Method to get the ID for a MicroRNAGene object.
#' @param gene An object of class MicroRNAGene.
#' @return A character string representing the ID of the gene.
#' @rdname ID-methods
#' @importFrom methods is
#' @export
setMethod("getID", "MicroRNAGene", function(gene) {
  if (!is(gene, "MicroRNAGene")) stop("'gene' must be of class 'MicroRNAGene'")
  return(gene@ID)
})

#' Generic setter for gene ID
#'
#' @description This method sets the ID for the gene.
#' @param gene An object of class Gene.
#' @param value A character string representing the new ID for the gene.
#' @return An updated Gene object with the new ID.
#' @rdname ID-methods
#' @export
setGeneric("setID<-", function(gene, value) standardGeneric("setID<-"))

#' Gene-specific setter for gene ID
#' 
#' @description Method to set the ID for a Gene object.
#' @param gene An object of class Gene.
#' @param value A character string representing the new ID for the gene.
#' @return An updated Gene object with the new ID.
#' @rdname ID-methods
#' @importFrom methods is
#' @export
setMethod("setID<-", "Gene", function(gene, value) {
  if (!is(gene, "Gene")) stop("'gene' must be of class 'Gene'")
  if (!is(value, "character")) stop("'value' must be of class 'character'")
  gene@ID <- value
  gene
})

#' ProteinCodingGene-specific setter for gene ID
#' 
#' @description Method to set the ID for a ProteinCodingGene object.
#' @param gene An object of class ProteinCodingGene.
#' @param value A character string representing the new ID for the gene.
#' @return An updated ProteinCodingGene object with the new ID.
#' @rdname ID-methods
#' @importFrom methods is
#' @export
setMethod("setID<-", "ProteinCodingGene", function(gene, value) {
  if (!is(gene, "ProteinCodingGene")) stop("'gene' must be of class 'ProteinCodingGene'")
  if (!is(value, "character")) stop("'value' must be of class 'character'")
  gene@ID <- value
  gene
})

#' LncRNAGene-specific setter for gene ID
#' 
#' @description Method to set the ID for a LncRNAGene object.
#' @param gene An object of class LncRNAGene.
#' @param value A character string representing the new ID for the gene.
#' @return An updated LncRNAGene object with the new ID.
#' @rdname ID-methods
#' @importFrom methods is
#' @export
setMethod("setID<-", "LncRNAGene", function(gene, value) {
  if (!is(gene, "LncRNAGene")) stop("'gene' must be of class 'LncRNAGene'")
  if (!is(value, "character")) stop("'value' must be of class 'character'")
  gene@ID <- value
  gene
})

#' MicroRNAGene-specific setter for gene ID
#' 
#' @description Method to set the ID for a MicroRNAGene object.
#' @param gene An object of class MicroRNAGene.
#' @param value A character string representing the new ID for the gene.
#' @return An updated MicroRNAGene object with the new ID.
#' @rdname ID-methods
#' @importFrom methods is
#' @export
setMethod("setID<-", "MicroRNAGene", function(gene, value) {
  if (!is(gene, "MicroRNAGene")) stop("'gene' must be of class 'MicroRNAGene'")
  if (!is(value, "character")) stop("'value' must be of class 'character'")
  gene@ID <- value
  gene
})
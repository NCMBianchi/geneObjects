#' Gene class
#'
#' This is a virtual S4 class to represent a general gene with associated information.
#'
#' @slot ID The gene ID (Ensembl ID or NCBI gene ID).
#' @slot HUGO_symbol The HUGO gene symbol.
#' @slot gene_name The full gene name.
#' @slot description A brief description of the gene.
#' @slot gene_structure The gene structure represented as a GRanges object.
#' @author Niccolò Bianchi\cr UniMi + PoliMI\cr Maintainer: Niccolò
#' Bianchi\cr E-Mail: <niccolo.bianchi2@@studenti.unimi.it>
#' @references \url{https://www.ensembl.org/index.html?redirect=no}\cr
#' @seealso \code{\link{ProteinCodingGene}}\cr
#' \code{\link{LncRNAGene}}\cr
#' \code{\link{MicroRNAGene}}\cr
#' @importFrom GenomicRanges GRanges
#' @examples
#' # Not usually instantiated directly, used as a superclass for other gene classes.
#'
#' @export
setClass("Gene", slots = list(ID = "character", HUGO_symbol = "character",
                              gene_name = "character", description = "character",
                              gene_structure = "GRanges"), contains = "VIRTUAL")


#' Protein-coding gene class
#'
#' This is an S4 class to represent a protein-coding gene, inheriting from the Gene class.
#'
#' @slot protein_ID The protein ID.
#' @slot protein_sequence The protein sequence.
#' @author Niccolò Bianchi\cr UniMi + PoliMI\cr Maintainer: Niccolò
#' Bianchi\cr E-Mail: <niccolo.bianchi2@@studenti.unimi.it>
#' @references \url{https://www.ensembl.org/index.html?redirect=no}\cr
#' @seealso \code{\link{LncRNAGene}}\cr
#' \code{\link{MicroRNAGene}}\cr
#' @importFrom GenomicRanges GRanges
#' @examples
#' gene_structure <- GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges::IRanges(start = 1000, end = 2000), strand = "+")
#' ProteinCodingGene(ID = "ENSG000001", HUGO_symbol = "BRCA1",
#'                   gene_name = "Breast cancer 1",
#'                   description = "Breast cancer type 1 susceptibility protein",
#'                   gene_structure = gene_structure,
#'                   protein_ID = "ENSP000003", protein_sequence = "MENSDRNSIKVAV...")
#'
#' @export
setClass("ProteinCodingGene", slots = list(protein_ID = "character",
                                           protein_sequence = "character"),
         contains = "Gene")


#' Long non-coding RNA gene class
#'
#' This is an S4 class to represent a long non-coding RNA gene, inheriting from the Gene class.
#'
#' @slot lncRNA_ID The lncRNA ID.
#' @slot RNA_sequence The RNA sequence.
#' @author Niccolò Bianchi\cr UniMi + PoliMI\cr Maintainer: Niccolò
#' Bianchi\cr E-Mail: <niccolo.bianchi2@@studenti.unimi.it>
#' @references \url{https://www.ensembl.org/index.html?redirect=no}\cr
#' @seealso \code{\link{ProteinCodingGene}}\cr
#' \code{\link{MicroRNAGene}}\cr
#' @importFrom GenomicRanges GRanges
#' @examples
#' gene_structure <- GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges::IRanges(start = 3000, end = 4000), strand = "-")
#' LncRNAGene(ID = "ENSG000002", HUGO_symbol = "LINC00001",
#'            gene_name = "Long intergenic non-protein coding RNA 1",
#'            description = "Long non-coding RNA",
#'            gene_structure = gene_structure,
#'            lncRNA_ID = "ENST000004", RNA_sequence = "ACUGCUAGCUAGUCA...")
#'
#' @export
setClass("LncRNAGene", slots = list(lncRNA_ID = "character",
                                    RNA_sequence = "character"), contains = "Gene")


#' MicroRNA gene class
#'
#' This is an S4 class to represent a microRNA gene, inheriting from the Gene class.
#'
#' @slot microRNA_ID The microRNA ID.
#' @slot microRNA_seed_sequence The microRNA seed sequence.
#' @author Niccolò Bianchi\cr UniMi + PoliMI\cr Maintainer: Niccolò
#' Bianchi\cr E-Mail: <niccolo.bianchi2@@studenti.unimi.it>
#' @references \url{https://www.ensembl.org/index.html?redirect=no}\cr
#' @seealso \code{\link{ProteinCodingGene}}\cr
#' \code{\link{LncRNAGene}}\cr
#' @importFrom GenomicRanges GRanges
#' @examples
#' gene_structure <- GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges::IRanges(start = 5000, end = 6000), strand = "+")
#' MicroRNAGene(ID = "ENSG000003", HUGO_symbol = "MIR00001",
#'              gene_name = "MicroRNA 1",
#'              description = "MicroRNA",
#'              gene_structure = gene_structure,
#'              microRNA_ID = "ENST000005", microRNA_seed_sequence = "UGAGGUAGUAGGUUGUAUGGUAG")
#'
#' @export
setClass("MicroRNAGene", slots = list(microRNA_ID = "character",
                                      microRNA_seed_sequence = "character"),
         contains = "Gene")
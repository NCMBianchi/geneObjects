% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/02_constructors.R
\name{LncRNAGene}
\alias{LncRNAGene}
\title{Constructor function for LncRNAGene class}
\usage{
LncRNAGene(ID, HUGO_symbol, gene_name, description, gene_structure, lncRNA_ID, RNA_sequence)
}
\arguments{
\item{ID}{The gene ID (Ensembl ID or NCBI gene ID).}

\item{HUGO_symbol}{The HUGO gene symbol.}

\item{gene_name}{The full gene name.}

\item{description}{A brief description of the gene.}

\item{gene_structure}{The gene structure represented as a GRanges object.}

\item{lncRNA_ID}{The lncRNA ID.}

\item{RNA_sequence}{The RNA sequence.}
}
\value{
An instance of the LncRNAGene class.
}
\description{
This function creates a new instance of the LncRNAGene class.
}
\examples{
gene_structure <- GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges::IRanges(3000, 4000), strand = "-")
LncRNAGene(ID = "ENSG000002", HUGO_symbol = "LINC00001",
           gene_name = "Long intergenic non-protein coding RNA 1",
           description = "Long non-coding RNA",
           gene_structure = gene_structure,
           lncRNA_ID = "ENST000004", RNA_sequence = "ACUGCUAGCUAGUCA...")
}
\references{
\url{https://www.ensembl.org/index.html?redirect=no}\cr
}
\author{
Niccolò Bianchi\cr UniMi + PoliMI\cr Maintainer: Niccolò
Bianchi\cr E-Mail: \href{mailto:niccolo.bianchi2@studenti.unimi.it}{niccolo.bianchi2@studenti.unimi.it}
}

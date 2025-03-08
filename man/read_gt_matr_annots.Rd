% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vcf_parsing_helpers.R
\name{read_gt_matr_annots}
\alias{read_gt_matr_annots}
\title{Get variant position and SnpEff annotations from Genotype Matrix}
\usage{
read_gt_matr_annots(gt_file, genes_bed_file = NULL)
}
\arguments{
\item{gt_file}{Unmodified output of `variantMatrixToPlink.pl.`}

\item{genes_bed_file}{necessary to get the BP-coordinate
at which collapsed variants of each gene are shown.
Format is tab-delimited chr, start, end, gene_ID. Do not use a header row.}
}
\description{
Get variant position and SnpEff annotations from Genotype Matrix
}

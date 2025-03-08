% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vcf_parsing_helpers.R
\name{parse_HGVS.c}
\alias{parse_HGVS.c}
\title{Parse HGVS codon code to a simpler format}
\usage{
parse_HGVS.c(df, col = "HGVS.c", coalesce_HGVS.a = TRUE)
}
\description{
Mostly intended for interpreting variants with "modifier" effect,
e.g. promoter regions
}

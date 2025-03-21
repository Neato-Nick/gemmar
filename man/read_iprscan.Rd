% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vcf_parsing_helpers.R
\name{read_iprscan}
\alias{read_iprscan}
\title{Read interproscan output}
\usage{
read_iprscan(
  filename,
  trim_cols = TRUE,
  run_GO = TRUE,
  run_pathways = FALSE,
  split_GO = FALSE,
  convert_hyphens = FALSE
)
}
\arguments{
\item{filename}{TSV output from InterProScan}
}
\description{
Read interproscan output
}

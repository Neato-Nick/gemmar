% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vcf_parsing_helpers.R
\name{generate_dilution_scheme}
\alias{generate_dilution_scheme}
\title{Pair of functions to map 2-fold dilution schemes to observed MIC values.
Generated in cooperation with chatGPT.
Dynamically generate the 2-fold step scheme.}
\usage{
generate_dilution_scheme(min_value, max_value)
}
\description{
Pair of functions to map 2-fold dilution schemes to observed MIC values.
Generated in cooperation with chatGPT.
Dynamically generate the 2-fold step scheme.
}
\examples{
dilution_scheme <- generate_dilution_scheme(min_value, max_value)
}

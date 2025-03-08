% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vcf_parsing_helpers.R
\name{contingency_fisher}
\alias{contingency_fisher}
\title{Fisher's exact test on any variant within a summarized count matrix}
\usage{
contingency_fisher(check_var, x = gt_polymorphic_interest_phen_counts)
}
\description{
Fisher's exact test on any variant within a summarized count matrix
}
\details{
Input format example:
\tabular{lrll}{
  \strong{Variant} \tab \strong{count} \tab \strong{FLU_32_res_named} \tab \strong{gt_num_named} \cr
  CP060341.1_215724 \tab 135 \tab Flu32.S \tab Ref.allele\cr
  CP060341.1_215724 \tab 491 \tab Flu32.R \tab Ref.allele\cr
  CP060341.1_215724 \tab   8 \tab Flu32.S \tab Alt.allele\cr
  CP060341.1_215724 \tab 228 \tab Flu32.R \tab Alt.allele
}
}
\examples{
fisher_outs <- map(check_specific_vars, contingency_fisher) \%>\% bind_rows()
}

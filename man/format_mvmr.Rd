% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/format_function.R
\name{format_mvmr}
\alias{format_mvmr}
\title{format_mvmr}
\usage{
format_mvmr(BXGs, BYG, seBXGs, seBYG, RSID)
}
\arguments{
\item{BXGs}{A matrix containing beta-coefficient values for genetic associations with the each exposure. Columns should indicate exposure number, with rows representing estimates for a given genetic variant.}

\item{BYG}{A numeric vector of beta-coefficient values for genetic associations with the outcome.}

\item{seBXGs}{A matrix containing standard errors corresponding to the matrix of beta-coefficients \code{BXGs}.}

\item{seBYG}{A numeric vector of standard errors corresponding to the beta-coefficients \code{BYG}.}

\item{RSID}{A vector of names for genetic variants included in the analysis. If variant IDs are not provided a vector of ID numbers will be generated.}
}
\value{
A formatted data frame of class \code{mvmr_format}.
}
\description{
Reads in summary data. Checks and organises columns for use in calculating multivariable Mendelian Randomization analyses. Where variant IDs are not provided, a vector is generated for variant identification.
}
\examples{
r_input <- format_mvmr(
  BXGs = rawdat_mvmr[, c("LDL_beta", "HDL_beta")],
  BYG = rawdat_mvmr$SBP_beta,
  seBXGs = rawdat_mvmr[, c("LDL_se", "HDL_se")],
  seBYG = rawdat_mvmr$SBP_se,
  RSID = rawdat_mvmr$SNP
)
names(r_input)
class(r_input)
}
\references{
Sanderson, E., et al., An examination of multivariable Mendelian randomization in the single-sample and two-sample summary data settings. International Journal of Epidemiology, 2019, 48, 3, 713--727. \doi{10.1093/ije/dyy262}
}
\author{
Wes Spiller; Eleanor Sanderson; Jack Bowden.
}

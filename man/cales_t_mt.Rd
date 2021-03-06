% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cal_nes_mt.R
\name{cales_t_mt}
\alias{cales_t_mt}
\title{calculate neoantigens NES for mutation type of a sample}
\usage{
cales_t_mt(
  data,
  barcode,
  calp = FALSE,
  cal_type = "exp",
  mhc_type = "I",
  IC50_threshold = 500,
  Rank_threshold = 10,
  type = "I",
  trim,
  DAI,
  DAI_threshold,
  sample_counts
)
}
\arguments{
\item{barcode}{Tumor sample barcode.}

\item{calp}{Bool value, calculate p value or not.}

\item{cal_type}{Should be "ccf" or "exp", calculate enrichment score in a aspect of Cancer Cell Fraction or Expression.}

\item{mhc_type}{Should be "I" or "II" ,calculate enrichmrnt score for MHC-I or MHC-II.}

\item{IC50_threshold}{Mutations below this threshold are considered as MHC-I neoantigens.}

\item{Rank_threshold}{Mutations below this threshold are considered as MHC-II neoantigens.}

\item{file}{Mutation file,folowing column are possible needed: MT_mean,exp,sample,chromosome,position,\%Rank_best_perL,ccf_cn_assume.}
}
\value{
es,nes,p
}
\description{
calculate neoantigens NES for mutation type of a sample
}

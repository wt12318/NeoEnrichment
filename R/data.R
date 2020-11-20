#' Cancer type code and TCGA Study Abbreviations
#'
#' A dataset containing the relationship between Cancer type code and TCGA Study Abbreviations
#'
#' @format A data frame with 830 rows and 3 variables:
#' \describe{
#'   \item{cancer_code}{two characters, representing 6 and 7 position of TCGA sample barcode}
#'   \item{cancer_type_full_name}{full name of TCGA cancer type}
#'   \item{cancer_type}{TCGA cancer type abbreviations}
#' }
#' @source \url{https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tissue-source-site-codes}
#' @source \url{https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations}
#' @export
"cancer_type_code"

#' Test data
#'
#' A dataset containing mutations (SNV) of 1061 samples from TCGA
#'
#' @format A data frame with 61414 rows and 9 variables:
#' \describe{
#'   \item{sample}{TCGA sample barcode}
#'   \item{chromosome}{chromosome, chr}
#'   \item{position}{position for mutation}
#'   \item{MT_mean}{aggregated IC50}
#'   \item{exp}{expression of gene that has this mutation, TPM}
#'   \item{ccf_cn_assume}{CCF of mutation}
#'   \item{cancer_type}{TCGA cancer type}
#'   \item{ref}{ref allele}
#'   \item{alt}{alt allele}
#' }
"mt_test"

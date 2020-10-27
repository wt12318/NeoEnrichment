#' get cancer type from tumor sample barcode
#'
#' @param tumor_barcode TCGA tumor sample barcode
#' @param cores How many cores to use if parallelization is used
#' @param par bool value, use parallelization or not
#'
#' @return cancer type
#' @export
#' @examples function("TCGA-OR-A5J1",par=FALSE)

getcancer_type <- function(tumor_barcode,cores,par=FALSE){

  cancer_type_code <- NeoEnrichment:::cancer_type_code
  if(!par){
    code <- substr(tumor_barcode,6,7)
    cancer_type <- sapply(code,
                          function(x){
                            cancer_type_code$cancer_type[grep(x,cancer_type_code$cancer_code)]
                          })
    return(cancer_type)
  }else{
    cat(paste0("you have ",detectCores(logical = F)," cores"))
    cl <- makeCluster(getOption("cl.cores", cores),type="FORK")
    code <- substr(tumor_barcode,6,7)
    cancer_type <- parSapply(cl=cl,code,
                             function(x){
                               cancer_type_code$cancer_type[grep(x,cancer_type_code$cancer_code)]
                             })
    stopCluster(cl)
    return(cancer_type)
  }

}


#' Get TCGA cancer type abbreviation
#'
#' Return TCGA cancer type abbreviation for TCGA sample barcode, this function use the 6 and 7 position
#' of sample barcode to extract cancer type based on the table of Tissue Source Site Codes of GDC, the description
#' of this data can be seen at document by `?cancer_type_code`
#'
#' @param input TCGA tumor sample barcode or TCGA cancer type full name
#' @param parallel Bool value, use parallelization or not, default FALSE
#' @param cores How many cores to use if parallelization is used, only can ce used when parallel is TRUE
#' @param need_full_name Bool value, return TCGA cancer type full name or not, default FALSE
#' @param input_type "code" or "full_name" for input of TCGA tumor sample barcode or TCGA cancer type full name, respectively.
#'
#' @return cancer type
#' @export
#' @examples
#' get_cancer_type("TCGA-OR-A5J1",par=FALSE)
#' get_cancer_type("Liver hepatocellular carcinoma",input_type="full_name")

get_cancer_type <- function(input,cores,parallel=FALSE,need_full_name=FALSE,input_type="code"){

  cancer_type_code <- NeoEnrichment::cancer_type_code
  if (input_type=="code"){
    if (!need_full_name){
      if(!parallel){
        code <- substr(input,6,7)
        cancer_type <- sapply(code,
                              function(x){
                                cancer_type_code$cancer_type[grep(x,cancer_type_code$cancer_code)]
                              })
        return(cancer_type)
      }else{
        cat(paste0("you have ",detectCores(logical = F)," cores"))
        cl <- makeCluster(getOption("cl.cores", cores),type="FORK")
        code <- substr(input,6,7)
        cancer_type <- parSapply(cl=cl,code,
                                 function(x){
                                   cancer_type_code$cancer_type[grep(x,cancer_type_code$cancer_code)]
                                 })
        stopCluster(cl)
        return(cancer_type)
      }
    }else{
      if(!parallel){
        code <- substr(input,6,7)
        cancer_type <- sapply(code,
                              function(x){
                                cancer_type_code$cancer_type_full_name[grep(x,cancer_type_code$cancer_code)]
                              })
        return(cancer_type)
      }else{
        cat(paste0("you have ",detectCores(logical = F)," cores"))
        cl <- makeCluster(getOption("cl.cores", cores),type="FORK")
        code <- substr(input,6,7)
        cancer_type <- parSapply(cl=cl,code,
                                 function(x){
                                   cancer_type_code$cancer_type_full_name[grep(x,cancer_type_code$cancer_code)]
                                 })
        stopCluster(cl)
        return(cancer_type)
      }
    }
  }else if (input_type=="full_name"){
    if(!parallel){
      cancer_type <- sapply(input,
                            function(x){
                              unique(cancer_type_code$cancer_type[grep(x,cancer_type_code$cancer_type_full_name)])
                            })
      return(cancer_type)
    }else{
      cat(paste0("you have ",detectCores(logical = F)," cores"))
      cl <- makeCluster(getOption("cl.cores", cores),type="FORK")
      cancer_type <- parSapply(cl=cl,input,
                               function(x){
                                 cancer_type_code$cancer_type[grep(x,cancer_type_code$cancer_type_full_name)]
                               })
      stopCluster(cl)
      return(cancer_type)
    }
  }else{
    cat("input_type must be 'code' or 'full_name' ")
  }
}


#' calculate neoantigens NES for a sample
#'
#' @param data Mutation file,folowing column are possible needed: MT_mean,exp,sample,chromosome,position,\%Rank_best_perL,ccf_cn_assume.
#' @param barcode Tumor sample barcode.
#' @param calp Bool value, calculate p value or not.
#' @param cal_type Should be "ccf" or "exp", calculate enrichment score in a aspect of Cancer Cell Fraction or Expression.

#' @return es,nes,p
#' @export
#' @importFrom rlang .data
#'
cales_t <- function(data,barcode,calp=FALSE,cal_type="exp",type="I",
                    sample_counts){
  if(cal_type=="exp"){

    file <- data %>%
      filter(!is.na(exp))

    test <- file %>% filter(sample==barcode)%>%
      dplyr::arrange(desc(.data$exp)) %>% dplyr::mutate(index=row_number())

    test$rank <- abs(nrow(test)/2 - seq(nrow(test):1))+1
  }else{

    file <- data %>%
      filter(!is.na(vaf)) %>%
      filter(!is.na(ccf))

    test <- file %>% dplyr::filter(sample==barcode)%>%
      dplyr::arrange(desc(.data$ccf),desc(.data$vaf)) %>%
      dplyr::mutate(index=row_number())
    test$rank <- abs(nrow(test)/2 - seq(nrow(test):1))+1
  }
  neo_list <- test %>% filter(.data$neo == "neo") %>% dplyr::select(.data$index)
  neo_list <- neo_list[[1]]


  if(length(neo_list)==0){
    return(paste(barcode,"no neoantigen"))
  }else{
    es <- cales(test,neo_list,type=type)

    if(calp==T){
      r <- cal_p_and_normalized(es,neo_list,test,type=type,sample_counts)
    }else{r <- es}
  }
}

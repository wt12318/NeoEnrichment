#' calculate neoantigens NES for a sample, combing ccf and exp
#'
#' @param data Mutation file,folowing column are possible needed: MT_mean,exp,sample,chromosome,position,\%Rank_best_perL,ccf_cn_assume.
#' @param barcode Tumor sample barcode.
#' @param calp Bool value, calculate p value or not.
#' @param mhc_type Should be "I" or "II" ,calculate enrichmrnt score for MHC-I or MHC-II.
#' @param IC50_threshold Mutations below this threshold are considered as MHC-I neoantigens.
#' @param Rank_threshold Mutations below this threshold are considered as MHC-II neoantigens.
#'
#' @return es,nes,p
#' @export
#' @importFrom rlang .data
#'
cales_t_combine <- function(data,barcode,calp=FALSE,mhc_type="I",IC50_threshold=500,Rank_threshold=10,type="I",trim,DAI,DAI_threshold,
                    sample_counts){
  file <- data %>%
    filter(!is.na(exp)) %>%
    filter(!is.na(ccf_cn_assume))

  test <- file %>% filter(sample==barcode)%>%
    dplyr::arrange(desc(.data$ccf_cn_assume),desc(.data$exp)) %>%
    dplyr::mutate(index=row_number())
  test$index <- rev(test$index)
  a <- nrow(test)
  test <- test %>%
    dplyr::mutate(rank=abs((a/2)-index)+1)


  if(DAI==TRUE){
    neo_list <- ifelse(mhc_type=="I",test %>% filter(.data$MT_mean<IC50_threshold & .data$DAI>DAI_threshold) %>% dplyr::select(.data$index),
                       test %>% filter(.data$`%Rank_best_perL`<Rank_threshold & .data$DAI>DAI_threshold) %>% dplyr::select(.data$index))
    neo_list <- neo_list[[1]]
  }else{
    neo_list <- ifelse(mhc_type=="I",test %>% filter(.data$MT_mean<IC50_threshold) %>% dplyr::select(.data$index),
                       test %>% filter(.data$`%Rank_best_perL`<Rank_threshold) %>% dplyr::select(.data$index))
    neo_list <- neo_list[[1]]
  }


  if(length(neo_list)==0){
    return(paste(barcode,"no neoantigen"))
  }else{
    es <- cales(test,neo_list,cal_type=NULL,type=type,trim=trim)

    if(calp==T){
      r <- cal_p_and_normalized(es,neo_list,test,cal_type=NULL,type=type,trim=trim,sample_counts)
    }else{r <- es}
  }
}

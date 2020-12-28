#' calculate neoantigens NES for mutation type of a sample
#'
#' @param file Mutation file,folowing column are possible needed: MT_mean,exp,sample,chromosome,position,\%Rank_best_perL,ccf_cn_assume.
#' @param barcode Tumor sample barcode.
#' @param calp Bool value, calculate p value or not.
#' @param cal_type Should be "ccf" or "exp", calculate enrichment score in a aspect of Cancer Cell Fraction or Expression.
#' @param mhc_type Should be "I" or "II" ,calculate enrichmrnt score for MHC-I or MHC-II.
#' @param IC50_threshold Mutations below this threshold are considered as MHC-I neoantigens.
#' @param Rank_threshold Mutations below this threshold are considered as MHC-II neoantigens.
#'
#' @return es,nes,p
#' @export
#' @importFrom rlang .data
#'
cales_t_mt <- function(data,barcode,calp=FALSE,cal_type="exp",mhc_type="I",IC50_threshold=500,Rank_threshold=10,type="I",trim,DAI,DAI_threshold,
                    sample_counts){
  file <- data %>%
    filter(!is.na(ccf_cn_assume))

  if(mhc_type=="I"){
    test <- file %>% dplyr::filter(sample==barcode)%>%
      dplyr::arrange(desc(.data$ccf_cn_assume)) %>% dplyr::mutate(index=row_number())
  }else{
    test <- file %>% dplyr::filter(sample==barcode)%>%
      dplyr::arrange(desc(.data$ccf_cn_assume)) %>% dplyr::mutate(index=row_number())
  }
  test <- test %>%
    dplyr::mutate(rank = as.numeric(factor(rank(.data$ccf_cn_assume))))
  if(trim==TRUE){
    a <- max(test$rank)
    test <- test %>%
      dplyr::mutate(rank=abs((a/2)-rank)+1)
    test$rank <- ifelse(test$ccf_cn_assume==1,(test$rank)/2,test$rank)
  }else{
    a <- nrow(test)
    test <- test %>%
      dplyr::mutate(rank=abs((a/2)-rank)+1)
  }

  r <- enframe(table(test$mutation_type))
  r_es <- vector("list",length = nrow(r))
  names(r_es) <- r$name

  for(i in seq_along(r_es)){
    if(DAI==TRUE){
      neo_list <- ifelse(mhc_type=="I",test %>% filter(.data$MT_mean<IC50_threshold & .data$DAI>DAI_threshold & .data$mutation_type==names(r_es)[i]) %>% dplyr::select(.data$index),
                         test %>% filter(.data$`%Rank_best_perL`<Rank_threshold & .data$DAI>DAI_threshold & .data$mutation_type==names(r_es)[i]) %>% dplyr::select(.data$index))
      neo_list <- neo_list[[1]]
    }else{
      neo_list <- ifelse(mhc_type=="I",test %>% filter(.data$MT_mean<IC50_threshold & .data$mutation_type==names(r_es)[i]) %>% dplyr::select(.data$index),
                         test %>% filter(.data$`%Rank_best_perL`<Rank_threshold & .data$mutation_type==names(r_es)[i]) %>% dplyr::select(.data$index))
      neo_list <- neo_list[[1]]
    }
    if(length(neo_list)==0){
      r_es[[i]] <- "no neoantigen"
    }else{
      es <- cales(test,neo_list,cal_type=cal_type,type=type,trim=trim)

      if(calp==T){
        nes <- cal_p_and_normalized(es,neo_list,test,cal_type=cal_type,type=type,trim=trim,sample_counts)
      }else{nes <- es}
      r_es[[i]] <- nes
    }
  }
  r_es <- Filter(function(x){length(x)>1},r_es)
  r_es <- bind_rows(r_es,.id = "mutation_type")
  return(r_es)
}

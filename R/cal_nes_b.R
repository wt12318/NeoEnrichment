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
cales_tb <- function(data,calp=FALSE,cal_type="exp",type="I",
                    sample_counts){
  if(cal_type=="exp"){

    file <- data %>%
      filter(!is.na(exp))

    test <- file %>% filter(sample==barcode)%>%
      dplyr::arrange(desc(.data$exp)) %>% dplyr::mutate(index=row_number())

    test <- test %>%
      dplyr::mutate(rank = as.numeric(factor(rank(.data$exp))))
    a <- nrow(test)
    test <- test %>%
      dplyr::mutate(rank=abs((a/2)-rank)+1)
  }else{

    # test_neo <- data %>%
    #   filter(neo %in% c("neo","not_neo")) %>%
    #   dplyr::arrange(desc(.data$ccf),desc(.data$dna_vaf)) %>%
    #   dplyr::mutate(index=row_number())
    # l_test_neo <- nrow(test_neo)
    # test_neo <- test_neo %>%
    #   dplyr::mutate(rank = abs((l_test_neo/2)-c(l_test_neo:1))+1)
    #
    # test_sys <- data %>%
    #   filter(neo %in% c("sys","not_neo")) %>%
    #   dplyr::arrange(desc(.data$ccf),desc(.data$dna_vaf)) %>%
    #   dplyr::mutate(index=row_number())
    # l_test_sys <- nrow(test_sys)
    # test_sys <- test_sys %>%
    #   dplyr::mutate(rank = abs((l_test_sys/2)-c(l_test_sys:1))+1)
    ll <- nrow(data)
    dt <- data %>%
      dplyr::arrange(desc(.data$ccf),desc(.data$dna_vaf)) %>%
      dplyr::mutate(index=row_number()) %>%
      dplyr::mutate(rank = abs((ll/2)-c(ll:1))+1)

  }
  neo_list <- dt %>% filter(.data$neo == "neo") %>% dplyr::select(.data$index)
  neo_list <- neo_list[[1]]

  sys_list <- dt %>% filter(.data$neo == "sys") %>% dplyr::select(.data$index)
  sys_list <- sys_list[[1]]

  if(length(neo_list)==0 | length(sys_list)==0){
    return(paste(unique(data$sample),"no neoantigen/sys"))
  }else{
    es_neo <- cales(dt,neo_list,type=type)
    es_sys <- cales(dt,sys_list,type=type)

    diff <- es_neo - es_sys
    r <- cal_p_and_normalized_b(diff,type=type,sample_counts,data=dt)
    r <- cbind(
      data.frame(es_neo=es_neo,es_sys=es_sys,diff=diff,sample=unique(data$sample)),
      r
    )
    return(r)
  }
}

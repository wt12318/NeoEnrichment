#' calculate neoantigens NES for a sample
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
#'
#'
cales_t <- function(file,barcode,calp=FALSE,cal_type="exp",mhc_type="I",IC50_threshold=500,Rank_threshold=10){
  if(cal_type=="exp"){

    file <- file %>%
      filter(!is.na(exp))

    if(mhc_type=="I"){
      test <- file %>% filter(sample==barcode)%>%
        dplyr::select(MT_mean,exp,sample,chromosome,position) %>%
        dplyr::mutate(index=paste(sample,chromosome,position,sep = ",")) %>%
        dplyr::distinct(index,.keep_all=TRUE) %>%
        dplyr::arrange(desc(exp)) %>% dplyr::mutate(index=row_number())
    }else{
      test <- file %>% filter(sample==barcode)%>%
        dplyr::select(`%Rank_best_perL`,exp,sample,chromosome,position) %>%
        dplyr::mutate(index=paste(sample,chromosome,position,sep = ",")) %>%
        dplyr::distinct(index,.keep_all=TRUE) %>%
        dplyr::arrange(desc(exp)) %>% dplyr::mutate(index=row_number())
    }
    a <- nrow(test)
    test <- test %>%
      dplyr::mutate(rank = rank(exp)) %>%
      dplyr::mutate(rank=abs((a/2)-rank)+1)
  }else{

    file <- file %>%
      filter(!is.na(ccf_cn_assume))

    if(mhc_type=="I"){
      test <- file %>% dplyr::filter(sample==barcode)%>%
        dplyr::select(MT_mean,sample,chromosome,position,ccf_cn_assume) %>%
        dplyr::mutate(index=paste(sample,chromosome,position,sep = ",")) %>%
        dplyr::distinct(index,.keep_all=T) %>%
        dplyr::arrange(desc(ccf_cn_assume)) %>% dplyr::mutate(index=row_number())
    }else{
      test <- file %>% dplyr::filter(sample==barcode)%>%
        dplyr::select(`%Rank_best_perL`,sample,chromosome,position,ccf_cn_assume) %>%
        dplyr::mutate(index=paste(sample,chromosome,position,sep = ",")) %>%
        dplyr::distinct(index,.keep_all=T) %>%
        dplyr::arrange(desc(ccf_cn_assume)) %>% dplyr::mutate(index=row_number())
    }
    a <- nrow(test)
    test <- test %>%
      dplyr::mutate(rank = rank(ccf_cn_assume)) %>%
      dplyr::mutate(rank=abs((a/2)-rank)+1)
  }
  neo_list <- ifelse(mhc_type=="I",test %>% filter(MT_mean<IC50_threshold) %>% dplyr::select(index),
                     test %>% filter(`%Rank_best_perL`<Rank_threshold) %>% dplyr::select(index))
  neo_list <- neo_list[[1]]
  if(length(neo_list)==0){
    return(paste(barcode,"no neoantigen"))
  }else{
    es <- cales(test,neo_list)
    if(calp==T){
      r <- cal_p_and_normalized(es,neo_list,test)
    }else{r <- es}
  }
}

#' calculate enrichment score for neoantigens
#'
#' @param mutation_dt Mutations dataframe
#' @param neoantigen_list Neoantigens list
#'
#' @return enrichment score
#' @export
#'
#'
#'
cales <- function(mutation_dt,neoantigen_list,cal_type,type){
  a <- sum(mutation_dt[mutation_dt$index %in% neoantigen_list,"rank"])
  if(type=="I"){
    b <- (nrow(mutation_dt)-length(neoantigen_list))
    mutation_dt$re <- ifelse(mutation_dt$index %in% neoantigen_list,((mutation_dt$rank)/a),-(1/b))

  }else{
    b <- sum(mutation_dt[!(mutation_dt$index %in% neoantigen_list),"rank"])
    mutation_dt$re <- ifelse(mutation_dt$index %in% neoantigen_list,((mutation_dt$rank)/a),-((mutation_dt$rank)/b))

  }
  mutation_dt$cum_re <- cumsum(mutation_dt$re)
  if(cal_type=="ccf"){
    mutation_dt_1 <- mutation_dt %>% dplyr::filter(ccf_cn_assume==1)
    if(nrow(mutation_dt_1)==0){
      es <- max(0,mutation_dt$cum_re)-abs(min(mutation_dt$cum_re,0))
    }else{
      mutation_dt_no_1 <- mutation_dt %>% dplyr::filter(ccf_cn_assume!=1)
      test <- rbind(mutation_dt_1[nrow(mutation_dt_1),],mutation_dt_no_1)
      es <- max(0,test$cum_re)-abs(min(test$cum_re,0))
      return(es)
    }
  }else{
    es <- max(0,mutation_dt$cum_re)-abs(min(mutation_dt$cum_re,0))
    return(es)
  }
}


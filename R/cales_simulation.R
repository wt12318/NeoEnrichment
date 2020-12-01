#' calculate es for simulation
#'
#' @param test Mutation file.
#' @param neoantigen_list Neoantigen list.
#' @export
#'
#' @return es
cales_simulation <- function(test,neoantigen_list,cal_type,type){
  a <- sum(test[test$mutation_id %in% neoantigen_list,"rank"])
  if(type=="I"){
    b <- (nrow(test)-length(neoantigen_list))
    test$re <- ifelse(test$mutation_id %in% neoantigen_list,((test$rank)/a),-(1/b))

  }else{
    b <- sum(test[!(test$mutation_id %in% neoantigen_list),"rank"])
    test$re <- ifelse(test$mutation_id %in% neoantigen_list,((test$rank)/a),-((test$rank)/b))

  }
  test$cum_re <- cumsum(test$re)
  if(cal_type=="ccf"){
    test_1 <- test %>% filter(ccf_cn_assume==1)
    if(nrow(test_1)==0){
      es <- max(0,test$cum_re)-abs(min(test$cum_re,0))
    }else{
      test_no1 <- test %>% filter(ccf_cn_assume!=1)
      test <- rbind(test_1[nrow(test_1),],test_no1)
      es <- max(0,test$cum_re)-abs(min(test$cum_re,0))
    }
    return(es)
  }else{
    es <- max(0,test$cum_re)-abs(min(test$cum_re,0))
    return(es)
  }
}

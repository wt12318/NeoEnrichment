#' calculate es for simulation
#'
#' @param test Mutation file.
#' @param neoantigen_list Neoantigen list.
#' @export
#'
#' @return es
cales_simulation <- function(test,neoantigen_list){
  a <- sum(test[test$mutation_id %in% neoantigen_list,"rank"])
  b <- (nrow(test)-length(neoantigen_list))

  test$re <- ifelse(test$mutation_id %in% neoantigen_list,((test$rank)/a),-(1/b))
  # test$re <- case_when(
  #   test$mutation_id %in% neoantigen_list & test$ccf_cn_assume==1 ~ ((test$rank)/a)*(1/2),
  #   test$mutation_id %in% neoantigen_list & test$ccf_cn_assume!=1 ~ ((test$rank)/a),
  #   !(test$mutation_id %in% neoantigen_list) & test$ccf_cn_assume==1 ~ -(1/b)*(1/2),
  #   !(test$mutation_id %in% neoantigen_list) & test$ccf_cn_assume!=1 ~ -(1/b)
  # )
  test$cum_re <- cumsum(test$re)
  test_1 <- test %>% filter(ccf_cn_assume==1)
  if(nrow(test_1)==0){
    es <- max(0,test$cum_re)-abs(min(test$cum_re,0))
  }else{
    test_no1 <- test %>% filter(ccf_cn_assume!=1)
    test <- rbind(test_1[nrow(test_1),],test_no1)
    es <- max(0,test$cum_re)-abs(min(test$cum_re,0))
  }
  return(es)
}

#' calculate es for simulation
#'
#' @param test Mutation file.
#' @param neoantigen_list Neoantigen list.
#'
#' @return es
cales_simulation <- function(test,neoantigen_list){
  a <- sum(test[test$mutation_id %in% neoantigen_list,"rank"])
  b <- (nrow(test)-length(neoantigen_list))

  test$re <- ifelse(test$mutation_id %in% neoantigen_list,((test$rank)/a),-(1/b))
  test$cum_re <- cumsum(test$re)
  es <- max(0,test$cum_re)-abs(min(test$cum_re,0))
}

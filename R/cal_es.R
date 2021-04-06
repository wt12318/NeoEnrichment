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
cales <- function(mutation_dt,neoantigen_list,type){
  a <- sum(mutation_dt[mutation_dt$index %in% neoantigen_list,"rank_ccf"])
  if(type=="I"){
    b <- (nrow(mutation_dt)-length(neoantigen_list))
    mutation_dt$re <- ifelse(mutation_dt$index %in% neoantigen_list,((mutation_dt$rank_ccf)/a),-(1/b))

  }else{
    b <- sum(mutation_dt[!(mutation_dt$index %in% neoantigen_list),"rank_ccf"])
    mutation_dt$re <- ifelse(mutation_dt$index %in% neoantigen_list,((mutation_dt$rank_ccf)/a),-((mutation_dt$rank)/b))
  }
  mutation_dt$cum_re <- cumsum(mutation_dt$re)
  es <- max(0,mutation_dt$cum_re)-abs(min(mutation_dt$cum_re,0))
}


#' calculate p value and normalized enrichment score for neoantigens
#'
#' @param es Enrichment score
#' @param neo_list Neoantigens list
#' @param mutation_dt Mutations dataframe
#'
#' @return enrichment score, enrichment score, p value
#' @export
#'
#'
#'
cal_p_and_normalized <- function(es,neo_list,mutation_dt,cal_type,type){
  sample_res <- data.frame(res=rep(1,1000))
  for (i in 1:1000) {
    neoantigen_list <- sample(mutation_dt$index,length(neo_list),replace = F)
    sample_res$res[i] <- cales(mutation_dt,neoantigen_list,cal_type=cal_type,type=type)
  }
  p <- ifelse(es<0,mean(sample_res$res<es),mean(sample_res$res>es))
  nes <- ifelse(es<0,es/abs(mean(sample_res$res[sample_res$res<0])),
                es/mean(sample_res$res[sample_res$res>0]))
  es_p <- data.frame(es=es,nes=nes,p_value=p)
}

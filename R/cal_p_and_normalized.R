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
cal_p_and_normalized <- function(es,neo_list,mutation_dt,cal_type,type,trim,sample_counts,nes_type="I"){
  sample_res <- data.frame(res=rep(1,sample_counts))
  for (i in 1:sample_counts) {
    neoantigen_list <- sample(mutation_dt$index,length(neo_list),replace = F)
    sample_res$res[i] <- cales(mutation_dt,neoantigen_list,cal_type=cal_type,type=type,trim = trim)
  }
  #p <- ifelse(es<0,mean(sample_res$res<es),mean(sample_res$res>es))
  p <- ifelse(es<0,(sum(sample_res$res<es)+1)/(sample_counts+1),
              (sum(sample_res$res>es)+1)/(sample_counts+1))
  if(nes_type=="I"){
    nes <- ifelse(es<0,es/abs(mean(sample_res$res[sample_res$res<0])),
                  es/mean(sample_res$res[sample_res$res>0]))
  }else{
    nes <- (es-mean(sample_res$res))/sd(sample_res$res)
  }

  es_p <- data.frame(es=es,nes=nes,p_value=p)
}

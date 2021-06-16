#' cal nes test
#'
#' @param dt dataframe
#'
#' @return es_p
#' @export
#'
#' @examples
cal_nes_new_test <- function(dt,sample_counts){

  ##neo nes
  es_neo <- cal_es_new_test(dt,type = "neo")
  ##not neo nes
  es_not_neo <- cal_es_new_test(dt,type = "not_neo")
  if (is.character(es_neo) | is.character(es_not_neo)){
    return(es_neo)
  }else{

    diff <- es_neo - es_not_neo

    sample_res <- data.frame(tmp_es_neo=rep(1,sample_counts),
                             tmp_es_not_neo=rep(1,sample_counts),
                             tmp_diff=rep(1,sample_counts))
    sample_res$tmp_es_not_neo <- es_not_neo
    neo_counts <- sum(dt$neo == "yes")
    sys_counts <- sum(dt$neo == "no")
    interval_mid <- seq(0,1,length.out=max(neo_counts,sys_counts)) %>%
      round(.,digits = 2) + 0.01
    interval_mid <- interval_mid[-length(interval_mid)]
    for (i in 1:sample_counts) {

      tmp_dt <- dt

      tmp_dt$ccf <- sample(interval_mid,nrow(dt),replace = TRUE)
      sample_res$tmp_es_neo[i] <-  cal_es_new_test(tmp_dt,type = "neo")
      sample_res$tmp_es_not_neo[i] <-  cal_es_new_test(tmp_dt,type = "not_neo")



    }
    sample_res$tmp_diff <- sample_res$tmp_es_neo - sample_res$tmp_es_not_neo

    p <- ifelse(diff<0,(sum(sample_res$tmp_diff <= diff)+1)/(sample_counts+1),
                (sum(sample_res$tmp_diff >= diff)+1)/(sample_counts+1))



    es_p <- data.frame(es_neo=es_neo,es_not_neo=es_not_neo,
                       diff=diff,diff_p=p,sample=unique(dt$sample))

    return(es_p)
  }
}

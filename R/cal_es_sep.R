#' cal nes test
#'
#' @param dt dataframe
#'
#' @return es_p
#' @export
#'
#' @examples
cal_es_sep <- function(dt,sample_counts,range){

  ##neo nes
  es_neo <- cal_es_new(dt,type = "neo")
  ##sys nes
  es_sys <- cal_es_new(dt,type = "not_neo")
  if (is.character(es_neo) | is.character(es_sys)){
    return(es_neo)
  }else{

    neo <- dt %>%
      dplyr::filter(neo == "yes")
    sys <- dt %>%
      dplyr::filter(neo == "no")

    neo_counts <- nrow(neo)
    sys_counts <- nrow(sys)

    interval_mid <- seq(range[1],range[2],0.01)

    res_neo <- vector("numeric",1000)
    res_sys <- vector("numeric",1000)

    for (i in 1:1000){

      tmp_neo <- neo
      tmp_neo$ccf <- sample(interval_mid,neo_counts,replace = TRUE)
      tmp_res_neo <- cal_es_new(tmp_neo,type = "neo")
      res_neo[i] <- tmp_res_neo

      tmp_sys <- sys
      tmp_sys$ccf <- sample(interval_mid,sys_counts,replace = TRUE)
      tmp_res_sys <- cal_es_new(tmp_sys,type = "not_neo")
      res_sys[i] <- tmp_res_sys

    }

    p_neo <-  ifelse(es_neo<0,(sum(res_neo <= es_neo)+1)/(sample_counts+1),
                     (sum(res_neo >= es_neo)+1)/(sample_counts+1))
    p_sys <-  ifelse(es_sys<0,(sum(res_sys <= es_sys)+1)/(sample_counts+1),
                     (sum(res_sys >= es_sys)+1)/(sample_counts+1))

    es_p <- data.frame(es_neo=es_neo,es_sys=es_sys,
                       p_neo=p_neo,p_sys=p_sys,sample=unique(dt$sample))

    return(es_p)
  }
}

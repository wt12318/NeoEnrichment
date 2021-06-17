#' cal nes test
#'
#' @param dt dataframe
#'
#' @return es_p
#' @export
#'
#' @examples
cal_nes_new_test <- function(dt,sample_counts){

  es <- cal_es_new_test(dt)

  sample_res <- vector(length = sample_counts)

  for (i in 1:sample_counts) {

    tmp_dt <- dt

    tmp_dt$neo <- sample(dt$neo,nrow(dt),replace = FALSE)
    sample_res[i] <- cal_es_new_test(tmp_dt)

  }

  p <- ifelse(es<0,(sum(sample_res <= es)+1)/(sample_counts+1),
              (sum(sample_res >= es)+1)/(sample_counts+1))

  es_p <- data.frame(es=es,p=p,sample=unique(dt$sample))

  return(es_p)
}

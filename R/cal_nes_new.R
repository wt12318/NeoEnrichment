#' cal nes test
#'
#' @param dt dataframe
#'
#' @return es_p
#' @export
#'
#' @examples
cal_nes_new <- function(dt,sample_counts){

  ##neo nes
  es_neo <- cal_es_new(dt,type = "neo")
  ##not neo nes
  es_not_neo <- cal_es_new(dt,type = "not_neo")
  if (is.character(es_neo) | is.character(es_not_neo)){
    return(es_neo)
  }else{

    diff <- es_neo - es_not_neo

    sample_res <- data.frame(tmp_es_neo=rep(1,sample_counts),
                             tmp_es_not_neo=rep(1,sample_counts),
                             tmp_diff=rep(1,sample_counts))
    #sample_res$tmp_es_not_neo <- es_not_neo
    #neo_counts <- sum(dt$neo=="yes")
    interval_mid <- seq(0.005,1,0.01)
    for (i in 1:sample_counts) {
      # tmp_dt <- dt %>%
      #   arrange(neo)
      # inter <- seq(range[1]+0.005,range[2],0.01)
      # tmp_dt$ccf <- c(sample(inter,sum(tmp_dt$neo=="no"),replace = T),
      #                 sample(inter,sum(tmp_dt$neo=="yes"),replace = T))

      tmp_dt <- dt
      #tmp_dt$ccf <- runif(nrow(dt))
      #tmp_dt$ccf[tmp_dt$neo=="yes"] <- runif(neo_counts)
      #tmp_dt$ccf[tmp_dt$neo=="yes"] <- sample(interval_mid,neo_counts,replace = TRUE)

      #tmp_dt$ccf[tmp_dt$neo=="yes"] <- sample(interval_mid,neo_counts,replace = TRUE)
      tmp_dt$ccf <- sample(interval_mid,nrow(tmp_dt),replace = TRUE)

      sample_res$tmp_es_not_neo[i] <- cal_es_new(tmp_dt,type = "not_neo")
      sample_res$tmp_es_neo[i] <- cal_es_new(tmp_dt,type = "neo")
      #sample_res$tmp_es_not_neo[i] <-  cal_es_new(tmp_dt,type = "not_neo")
    }
    sample_res$tmp_diff <- sample_res$tmp_es_neo - sample_res$tmp_es_not_neo

    p <- ifelse(diff<0,(sum(sample_res$tmp_diff <= diff)+1)/(sample_counts+1),
                (sum(sample_res$tmp_diff >= diff)+1)/(sample_counts+1))

    # diff_norm <- case_when(
    #   diff<0 ~ diff/abs(mean(sample_res[sample_res$tmp_diff < 0,"tmp_diff"])),
    #   diff>0 ~ diff/abs(mean(sample_res[sample_res$tmp_diff > 0,"tmp_diff"])),
    #   diff==0 ~ 0
    # )

    #diff_norm1 <- (diff-mean(sample_res$tmp_diff))/sd(sample_res$tmp_diff)

    # nes_neo <- case_when(
    #   es_neo<0 ~ es_neo/abs(mean(sample_res[sample_res$tmp_es_neo < 0,"tmp_es_neo"])),
    #   es_neo>0 ~ es_neo/abs(mean(sample_res[sample_res$tmp_es_neo > 0,"tmp_es_neo"])),
    #   es_neo==0 ~ 0
    # )
    #
    # nes_not_neo <- case_when(
    #   es_not_neo<0 ~ es_not_neo/abs(mean(sample_res[sample_res$tmp_es_not_neo < 0,"tmp_es_not_neo"])),
    #   es_not_neo>0 ~ es_not_neo/abs(mean(sample_res[sample_res$tmp_es_not_neo > 0,"tmp_es_not_neo"])),
    #   es_not_neo==0 ~ 0
    # )

    es_p <- data.frame(es_neo=es_neo,es_not_neo=es_not_neo,
                       diff=diff,diff_p=p,sample=unique(dt$sample))

    return(es_p)
  }
}

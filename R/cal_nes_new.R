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
    # sample_res$tmp_es_not_neo <- es_not_neo
    # neo_counts <- sum(dt$neo == "yes")
    # #sys_counts <- sum(dt$neo != "yes")
    # interval_mid <- seq(0.005,1,0.01)
    for (i in 1:sample_counts) {
      # tmp_dt <- dt %>%
      #   arrange(neo)
      # inter <- seq(range[1]+0.005,range[2],0.01)
      # tmp_dt$ccf <- c(sample(inter,sum(tmp_dt$neo=="no"),replace = T),
      #                 sample(inter,sum(tmp_dt$neo=="yes"),replace = T))

      tmp_dt <- dt
      #tmp_dt$ccf <- sample(dt$ccf,nrow(dt),replace = F)
      # sample_res$tmp_es_not_neo[i] <-  cal_es_new(tmp_dt,type = "not_neo")
      # sample_res$tmp_es_neo[i] <-  cal_es_new(tmp_dt,type = "neo")
      #tmp_sys <- dt[sample(1:nrow(dt),sys_counts),] %>% mutate(neo="no")
      #sample_res$tmp_es_not_neo[i] <- cal_es_new(tmp_sys,type = "not_neo")
      #tmp_dt$ccf <- runif(nrow(dt))
      #tmp_dt$ccf[tmp_dt$neo=="yes"] <- runif(neo_counts)
      #tmp_dt$ccf[tmp_dt$neo=="yes"] <- sample(interval_mid,neo_counts,replace = TRUE)
      #tmp_dt$ccf <- sample(interval_mid,nrow(dt),replace = TRUE)
      tmp_dt$neo <- sample(dt$neo,nrow(dt),replace = FALSE)
      sample_res$tmp_es_neo[i] <-  cal_es_new(tmp_dt,type = "neo")
      sample_res$tmp_es_not_neo[i] <-  cal_es_new(tmp_dt,type = "not_neo")


      #tpm_dt$ccf[tpm_dt$neo=="yes"] <- sample(interval_mid,neo_counts,replace = TRUE)
      #tpm_dt$ccf[tpm_dt$neo=="yes"] <- sample(dt$ccf,neo_counts)
      #tpm_dt$ccf <- sample(interval_mid,nrow(tpm_dt),replace = TRUE)

      #tpm_dt$neo <- sample(dt$neo,nrow(dt))
      #sample_res$tmp_es_not_neo[i] <- cal_es_new(tpm_dt,type = "not_neo")
      # tmp_neo <- dt[sample(1:nrow(dt),neo_counts),] %>% mutate(neo="yes")
      # sample_res$tmp_es_neo[i] <- cal_es_new(tmp_neo,type = "neo")
      #sample_res$tmp_es_not_neo[i] <-  cal_es_new(tmp_dt,type = "not_neo")
    }
    sample_res$tmp_diff <- sample_res$tmp_es_neo - sample_res$tmp_es_not_neo

    p <- ifelse(diff<0,(sum(sample_res$tmp_diff <= diff)+1)/(sample_counts+1),
                (sum(sample_res$tmp_diff >= diff)+1)/(sample_counts+1))

    #diff_norm <- diff/(max(sample_res$tmp_diff)-min(sample_res$tmp_diff))
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

#' cal nes test
#'
#' @param dt dataframe
#'
#' @return
#' @export
#'
#' @examples
cal_nes_new <- function(dt){

  interval <- seq(1,0,length.out=100)

  neo <- dt %>%
    filter(MT_mean < 500)

  tmp <- table(cut(neo$ccf,interval)) %>% as.data.frame()
  tmp$Freq[1] <- tmp$Freq[1] + sum(neo$ccf==0)

  tmp <- tmp %>%
    mutate(index = c(1:nrow(tmp))) %>%
    arrange(desc(index)) %>%
    mutate(rank=abs((100/2)-index)+1)

  add <- 1/sum(tmp$Freq*tmp$rank)
  minus <- 1/(100-sum(tmp$Freq!=0))

  tmp$re <- ifelse(tmp$Freq == 0,-minus,tmp$Freq*tmp$rank * add)
  tmp$cum_re <- cumsum(tmp$re)

  es <- max(0,tmp$cum_re)-abs(min(tmp$cum_re,0))

  neo_counts <- sum(dt$MT_mean<500)

  re <- vector("numeric",1000)
  # for (i in 1:1000){
  #   tmp_dt <- sample(tmp$Var1,neo_counts,
  #                    replace = TRUE) %>% table %>% as.data.frame(stringsAsFactors=FALSE)
  #   tmp_dt <- tmp_dt %>%
  #     mutate(index = c(1:nrow(tmp_dt))) %>%
  #     arrange(desc(index)) %>%
  #     mutate(rank=abs((100/2)-index))
  #
  #   add1 <- 1/sum(tmp_dt$Freq*tmp_dt$rank)
  #   minus1 <- 1/(100-sum(tmp_dt$Freq!=0))
  #
  #   tmp_dt$re <- ifelse(tmp_dt$Freq == 0,-minus1,tmp_dt$Freq*tmp_dt$rank * add1)
  #   tmp_dt$cum_re <- cumsum(tmp_dt$re)
  #
  #   re[[i]] <- max(0,tmp_dt$cum_re)-abs(min(tmp_dt$cum_re,0))
  # }

  re <- mapply(
    function(x,...){
      tmp_dt <- sample(tmp$Var1,neo_counts,
                       replace = TRUE) %>% table %>% as.data.frame(stringsAsFactors=FALSE)
      tmp_dt <- tmp_dt %>%
        mutate(index = c(1:nrow(tmp_dt))) %>%
        arrange(desc(index)) %>%
        mutate(rank=abs((100/2)-index)+1)

      add1 <- 1/sum(tmp_dt$Freq*tmp_dt$rank)
      minus1 <- 1/(100-sum(tmp_dt$Freq!=0))

      tmp_dt$re <- ifelse(tmp_dt$Freq == 0,-minus1,tmp_dt$Freq*tmp_dt$rank * add1)
      tmp_dt$cum_re <- cumsum(tmp_dt$re)

      return(max(0,tmp_dt$cum_re)-abs(min(tmp_dt$cum_re,0)))
    },c(1:1000),MoreArgs=list(tmp,neo_counts)
  )

  p <- ifelse(es<0,(sum(re<es)+1)/(1000+1),
              (sum(re>es)+1)/(1000+1))
  nes <- ifelse(es<0,es/abs(mean(re[re<0])),
                es/mean(re[re>=0]))

  es_p <- data.frame(es=es,nes=nes,p_value=p)

  return(es_p)
}

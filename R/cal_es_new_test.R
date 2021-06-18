#' cal nes test
#'
#' @param dt dataframe
#' @param type "neo" or "not neo"
#'
#' @return
#' @export
#'
#' @examples
cal_es_new_test <- function(dt){


  interval <- seq(0,1,length.out=100)
  obj_neo <- dt %>%
    dplyr::filter(neo == "yes")
  obj_sys <- dt %>%
    dplyr::filter(neo == "no")

  tmp_obj_neo <- table(cut(obj_neo$ccf,interval)) %>% as.data.frame()
  tmp_obj_sys <- table(cut(obj_sys$ccf,interval)) %>% as.data.frame()
  if (nrow(tmp_obj_neo)==0){
    return("there are no neoantigen/no_neoantigen in this sample")
  }else{
    ll <- nrow(tmp_obj_neo)
    tmp_obj_neo$Freq[1] <- tmp_obj_neo$Freq[1] + sum(obj_neo$ccf==0)
    tmp_obj_neo <- tmp_obj_neo %>%
      mutate(index = c(1:nrow(tmp_obj_neo))) %>%
      arrange(desc(index))%>%
      mutate(rank=abs((ll/2)-index)+1)

    ll2 <- nrow(tmp_obj_sys)
    tmp_obj_sys$Freq[1] <- tmp_obj_sys$Freq[1] + sum(obj_sys$ccf==0)
    tmp_obj_sys <- tmp_obj_sys %>%
      mutate(index = c(1:nrow(tmp_obj_sys))) %>%
      arrange(desc(index)) %>%
      mutate(rank=abs((ll2/2)-index)+1)

    add <- 1/sum(tmp_obj_neo[tmp_obj_neo$Freq!=0,"rank"]*tmp_obj_neo[tmp_obj_neo$Freq!=0,"Freq"])
    add2 <- 1/sum(tmp_obj_sys[tmp_obj_sys$Freq!=0,"rank"]*tmp_obj_sys[tmp_obj_sys$Freq!=0,"Freq"])

    tmp_obj_neo$re <- ifelse(tmp_obj_neo$Freq == 0,0,
                         ((tmp_obj_neo$rank)*tmp_obj_neo$Freq)*add)
    tmp_obj_neo$cum_re <- cumsum(tmp_obj_neo$re)

    tmp_obj_sys$re <- ifelse(tmp_obj_sys$Freq == 0,0,
                             ((tmp_obj_sys$rank)*tmp_obj_sys$Freq)*add2)
    tmp_obj_sys$cum_re <- cumsum(tmp_obj_sys$re)

    diff <- tmp_obj_neo$cum_re - tmp_obj_sys$cum_re

    es <- max(0,diff)-abs(min(diff,0))
    return(es)
  }
}

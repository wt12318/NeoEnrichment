#' cal nes test
#'
#' @param dt dataframe
#' @param type "neo" or "not neo"
#'
#' @return
#' @export
#'
#' @examples
cal_es_new_test <- function(dt,type){

  n_neo <- sum(dt$neo=="yes")
  n_sys <- sum(dt$neo=="no")
  interval <- seq(0,1,length.out=max(n_neo,n_sys)) %>%
    round(.,digits = 2)

  if (type == "neo"){
    obj <- dt %>%
      dplyr::filter(neo == "yes")
  }else{
    obj <- dt %>%
      dplyr::filter(neo == "no")
  }

  tmp_obj <- table(cut(obj$ccf,interval)) %>% as.data.frame()

  if (nrow(obj)==0){
    return("there are no neoantigen/no_neoantigen in this sample")
  }else{
    ll <- nrow(tmp_obj)
    tmp_obj$Freq[1] <- tmp_obj$Freq[1] + sum(obj$ccf==0)
    tmp_obj <- tmp_obj %>%
      mutate(index = c(1:nrow(tmp_obj))) %>%
      arrange(desc(index))%>%
      mutate(rank=abs((ll/2)-index)+1)

    add <- 1/sum(tmp_obj[tmp_obj$Freq!=0,"rank"]*tmp_obj[tmp_obj$Freq!=0,"Freq"])
    minus <- 1/sum(tmp_obj$Freq==0)

    tmp_obj$re <- ifelse(tmp_obj$Freq == 0,-minus,
                         ((tmp_obj$rank)*tmp_obj$Freq)*add)
    tmp_obj$cum_re <- cumsum(tmp_obj$re)
    es <- max(0,tmp_obj$cum_re)-abs(min(tmp_obj$cum_re,0))
    return(es)
  }
}

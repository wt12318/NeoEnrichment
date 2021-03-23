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
cal_p_and_normalized_b <- function(diff,type,sample_counts,data){
  sample_res <- data.frame(res_neo=rep(1,sample_counts),res_sys=rep(1,sample_counts),
                           res_diff=rep(1,sample_counts))
  for (i in 1:sample_counts) {
    tmp_data <- data
    tmp_data$neo <- sample(data$neo,nrow(data))

    tmp_neo_list <- tmp_data %>% filter(.data$neo == "neo") %>% dplyr::select(.data$index)
    tmp_neo_list <- tmp_neo_list[[1]]

    tmp_sys_list <- tmp_data %>% filter(.data$neo == "sys") %>% dplyr::select(.data$index)
    tmp_sys_list <- tmp_sys_list[[1]]

    sample_res$res_neo[i] <- cales(tmp_data,tmp_neo_list,type=type)
    sample_res$res_sys[i] <- cales(tmp_data,tmp_sys_list,type=type)

    sample_res$res_diff[i] <- sample_res$res_neo[i] - sample_res$res_sys[i]
  }
  #p <- ifelse(es<0,mean(sample_res$res<es),mean(sample_res$res>es))
  p <- ifelse(diff<0,(sum(sample_res$res_diff<=diff)+1)/(sample_counts+1),
              (sum(sample_res$res_diff>=diff)+1)/(sample_counts+1))
  # diff_norm <- ifelse(diff<0,diff/abs(mean(sample_res$res_diff[sample_res$res_diff<0])),
  #               diff/mean(sample_res$res_diff[sample_res$res_diff>0]))
  diff_norm <- case_when(
    diff<0 ~ diff/abs(mean(sample_res$res_diff[sample_res$res_diff<0])),
    diff>0 ~ diff/abs(mean(sample_res$res_diff[sample_res$res_diff>0])),
    diff==0 ~ 0
  )

  es_p <- data.frame(diff_norm=diff_norm,p=p)
}

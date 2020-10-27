#' run the function of calculating NES
#'
#' @param cancer mutation file
#' @param calp bool value ,calculate p value or not
#' @param cal_type "exp" or "ccf"
#' @param mhc_type "I" or "II"
#' @param IC50_threshold threshold for MHC-I neoantigens
#' @param Rank_threshold threshold for MHC-II neoantigens
#'
#' @return a dataframe with es, nes, p_value
#'
#'
run_cal_nes <- function(cancer,calp=F,cal_type="exp",mhc_type="I",IC50_threshold,Rank_threshold){
  if(cal_type=="exp"){
    cancer <- cancer %>% filter(!is.na(exp))
  }else{cancer <- cancer %>% filter(!is.na(ccf_cn_assume))}

  sample <- unique(cancer$sample)
  results <- data.frame(es=c(1:length(sample)))
  for (i in 1:length(sample)){
    results$es[i] <- cales_t(cancer,sample[i],calp = calp,cal_type=cal_type,mhc_type=mhc_type,IC50_threshold = IC50_threshold,Rank_threshold=Rank_threshold)
  }
  results <- results %>% separate(es,into = c("es","nes","p_value"),sep = ",")
  results$sample <- sample
  results <- na.omit(results)
  results$es <- as.numeric(results$es)
  results$nes <- as.numeric(results$nes)
  results$p_value <- as.numeric(results$p_value)
  return(results)
}

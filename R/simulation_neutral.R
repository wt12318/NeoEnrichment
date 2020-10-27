#' simulate neutral distribution of NES by cancer type
#'
#' @param mutation_file  At least contains columns of MT_mean(IC50) and cancer type for sampling and other columns for calculating NES
#'
#' @return cancer type neutral simulation
#' @export
#'
#' @examples neutral_simulation <- simulation_neutral(mt_test)
simulation_neutral <- function(mutation_file){
  if(!("cancer_type" %in% colnames(mutation_file))){
    stop("there is no column named \"cancer_type\",please provide it (maybe you need use get_cancer_type() function from this package)")
  }
  cancer_type <- unique(mutation_file$cancer_type)

  sample_base <- data.frame(sample_id=NA,sample_mt=NA,cancer_type=NA,es=NA,nes=NA)
  for (i in 1:length(cancer_type)){
    cancer <- mutation_file %>% filter(cancer_type==cancer_type[i])
    cancer %>% group_by(sample) %>%
      summarise(mutation_counts=n()) -> mutation_counts
    ##sampling mutation counts
    sample_mt <- sample(mutation_counts$mutation_counts,10000,replace = T)

    sample <- data.frame(sample_id=c(1:10000),sample_mt=sample_mt,cancer_type=cancer_type[i],es_nes=NA)

    for(j in 1:nrow(sample)){
      tmp_sample <- data.frame(mutation_id=c(1:sample$sample_mt[j]),
                               IC50=sample(cancer$MT_mean,sample$sample_mt[j]),
                               rank=c(sample$sample_mt[j]:1))
      tmp_sample <- tmp_sample %>% dplyr::mutate(rank=abs((nrow(tmp_sample)/2)-rank)+1)
      neoantigen_list <- tmp_sample %>% filter(IC50<500) %>% select(mutation_id)
      neoantigen_list <- neoantigen_list[[1]]
      if(length(neoantigen_list)==0){
        cat("no neoantigen for sample",j,"in",cancer_type[i],"\n",sep = " ")
        sample$es_nes[j] <- NA
      }else{
        es <- cales_simulation(tmp_sample,neoantigen_list)
        nes <- cal_p_and_normalized_simulation(es,neoantigen_list,tmp_sample)
        sample$es_nes[j] <- nes
      }
    }

    sample <- sample %>% filter(!is.na(es_nes)) %>%
      separate(es_nes,into = c("es","nes"),sep = ",")
    sample$es <- as.numeric(sample$es)
    sample$nes <- as.numeric(sample$nes)
    sample_base <- rbind(sample_base,sample)
  }

  sample_base <- na.omit(sample_base)
  return(sample_base)
}

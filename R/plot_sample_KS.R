#' plot K-S like graph for a single sample
#'
#' @param mutation_file Mutations with neoantigen prediction (IC50).
#' @param barcode TCGA sample barcode.
#' @param mhc_type MHC type ,can be I or II.
#' @param IC50_threshold Threshold which IC50 below it can be thought as neoantigen
#'
#' @return NULL
#' @export
#'
plot_KS <- function(mutation_file,barcode,mhc_type="I",IC50_threshold=500){

  test <- mutation_file %>% dplyr::filter(sample==barcode)%>%
    dplyr::filter(!is.na(ccf_cn_assume)) %>%
    dplyr::select(MT_mean,sample,chromosome,position,ccf_cn_assume) %>%
    dplyr::mutate(index=paste(sample,chromosome,position,sep = ",")) %>%
    dplyr::distinct(index,.keep_all=T) %>%
    dplyr::arrange(desc(ccf_cn_assume)) %>% dplyr::mutate(index=row_number()) %>%
    dplyr::mutate(rank = rank(ccf_cn_assume))
  test <- test %>% dplyr::mutate(rank=abs((nrow(test)/2)-rank)+1)
  neo_list <- ifelse(mhc_type=="I",test %>% filter(MT_mean<IC50_threshold) %>% dplyr::select(index),
                     test %>% filter(`%Rank_best_perL`<Rank_threshold) %>% dplyr::select(index))
  neo_list <- neo_list[[1]]

  mutation_dt <- test
  a <- sum(mutation_dt[mutation_dt$index %in% neo_list,"rank"])
  b <- (nrow(mutation_dt)-length(neo_list))
  mutation_dt <- mutation_dt %>%
    mutate(`Neoantigen Set`=NA,`Non Neoantigen Set`=NA)
  mutation_dt$`Neoantigen Set`[1] <- ifelse(mutation_dt$index[1] %in% neo_list,(mutation_dt$rank[1])/a,0)
  mutation_dt$`Non Neoantigen Set`[1] <- ifelse(mutation_dt$index[1] %in% neo_list,0,1/b)

  for (i in 2:nrow(mutation_dt)){
    mutation_dt$`Neoantigen Set`[i] <- ifelse(mutation_dt$index[i] %in% neo_list,
                                              mutation_dt$`Neoantigen Set`[i-1]+(mutation_dt$rank[i])/a,
                                              mutation_dt$`Neoantigen Set`[i-1])
    mutation_dt$`Non Neoantigen Set`[i] <- ifelse(mutation_dt$index[i] %in% neo_list,mutation_dt$`Non Neoantigen Set`[i-1],
                                                  mutation_dt$`Non Neoantigen Set`[i-1]+1/b)
  }
  mutation_dt <- mutation_dt %>%
    pivot_longer(cols = c("Neoantigen Set","Non Neoantigen Set"),names_to="type",values_to="value")
  es <- NeoEnrichment:::cales(test,neo_list)
  nes <- NeoEnrichment:::cal_p_and_normalized(es = es,neo_list = neo_list,mutation_dt = test) %>%
    strsplit(.,",") %>% unlist %>% as.numeric() %>% format(., scientific = TRUE,digits = 3)
  p <- ggplot(mutation_dt) +
    geom_step( aes(x=index, y=value,color=type),size=1.5)+
    theme_classic()+
    annotate("text", x=quantile(mutation_dt$index)[2], y=0.75, label=paste0("ES=",nes[1],"\nNES=",nes[2],"\np=",nes[3]), size=6)
  return(p)
}

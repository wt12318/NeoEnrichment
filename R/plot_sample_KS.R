#' Plot K-S like graph
#'
#' Plot K-S like graph of neoantigen enrichment for a single sample. The x axis is mutation index ordered
#' by descending mutation CCF. The y axis is the cumulative value for mutations. Walking through the mutation
#' list (descend by ccf), when encounter a mutation is naoantigen, add a value (the rank of this mutation devided by sum of rank of all mutations)
#' to the distribution of Neoantigen Set, otherwise add a value (number of non neoantigen mutations) to the distribution of
#' Non Neoantigen Set. \cr
#' This plot can be used to compare two distributions and calculate Neoantigen Enrichment score.
#'
#' @param mutation_file Mutations with neoantigen prediction (IC50).
#' @param barcode TCGA sample barcode.
#' @param mhc_type MHC type ,can be I or II.
#' @param IC50_threshold Threshold which IC50 below it can be thought as neoantigen
#'
#' @return NULL
#' @export
#'
#' @examples
#' plot_KS(mutation_file = mt_test,barcode = "TCGA-OR-A5J1-01A")
#'
#' @importFrom rlang .data
#'
plot_KS <- function(mutation_file,barcode,mhc_type="I",IC50_threshold=500){

  test <- mutation_file %>% dplyr::filter(sample==barcode)%>%
    dplyr::filter(!is.na(.data$ccf_cn_assume)) %>%
    dplyr::select(.data$MT_mean,.data$sample,.data$chromosome,.data$position,.data$ccf_cn_assume) %>%
    dplyr::mutate(index=paste(.data$sample,.data$chromosome,.data$position,sep = ",")) %>%
    dplyr::distinct(.data$index,.keep_all=T) %>%
    dplyr::arrange(desc(.data$ccf_cn_assume)) %>% dplyr::mutate(index=row_number()) %>%
    dplyr::mutate(rank = rank(.data$ccf_cn_assume))
  test <- test %>% dplyr::mutate(rank=abs((nrow(test)/2)-rank)+1)
  neo_list <- ifelse(mhc_type=="I",test %>% filter(.data$MT_mean<IC50_threshold) %>% dplyr::select(.data$index),
                     test %>% filter(.data$`%Rank_best_perL`<Rank_threshold) %>% dplyr::select(.data$index))
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
  es <- NeoEnrichment::cales(test,neo_list)
  nes <- NeoEnrichment::cal_p_and_normalized(es = es,neo_list = neo_list,mutation_dt = test)
  nes[,1:3] <- apply(nes[,1:3],2,function(x){format(x, scientific = TRUE,digits = 3)})
  p <- ggplot(mutation_dt) +
    geom_step( aes(x=index, y=value,color=type),size=1.5)+
    theme_classic()+
    annotate("text", x=quantile(mutation_dt$index)[2], y=0.75, label=paste0("ES=",nes$es,"\nNES=",nes$nes,"\np=",nes$p_value), size=6)
  return(p)
}

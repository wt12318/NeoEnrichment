#' get mutation type (SBS96) for mutation data
#'
#' @param mutation mutation data, at least contains following columns: position, chromosome
#'
#' @return mutation data with mutation type
#' @export
#' @import data.table
#'
get_mutation_type <- function(mutation){

  hsgs.installed = BSgenome::installed.genomes(splitNameParts = TRUE)
  data.table::setDT(x = hsgs.installed)
  ref_genome <- "hg38"

  ref_genome = BSgenome::getBSgenome(genome = ref_genome)
  mutation$position <- as.numeric(mutation$position)
  mutation <- mutation %>%
    mutate(Start = position-1,End = position+1)
  extract.tbl <- data.table::setDT(mutation)

  ss = BSgenome::getSeq(x = ref_genome, names = extract.tbl[,chromosome], start = extract.tbl[,Start] , end = extract.tbl[,End], as.character = TRUE)
  extract.tbl[,trinucleotide:= as.character(ss)]

  sub <- data.frame(x=c("T","C","G","A"),y=c("A","G","C","T"),stringsAsFactors = F)

  extract.tbl[,`:=`(a=sapply(trinucleotide,
                             function(z){sub$y[which(sub$x==substr(z,1,1))]}),
                    b=sapply(trinucleotide,
                             function(z){sub$y[which(sub$x==substr(z,2,2))]}),
                    c=sapply(trinucleotide,
                             function(z){sub$y[which(sub$x==substr(z,3,3))]}))]
  extract.tbl[,`:=`(reverse=stringi::stri_reverse(paste(a,b,c,sep = "")))]
  extract.tbl <- extract.tbl[,-c("a","b","c")]
  extract.tbl[,alt_reverse:=sapply(alt,function(z){sub$y[which(sub$x==z)]})]

  ##get mutation type
  extract.tbl[,mutation_type:=ifelse(ref %in% c("C","T"),
                                     paste(substr(trinucleotide,1,1),"[",
                                           substr(trinucleotide,2,2),">",
                                           alt,"]",substr(trinucleotide,3,3),sep = ""),
                                     paste(substr(reverse,1,1),"[",
                                           substr(reverse,2,2),">",
                                           alt_reverse,"]",substr(reverse,3,3),sep = ""))]
  extract.tbl <- extract.tbl %>%
    select(-c("Start","End","trinucleotide","reverse","alt_reverse"))
  return(extract.tbl)
}


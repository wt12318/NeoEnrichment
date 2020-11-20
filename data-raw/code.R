## code to prepare `code` dataset goes here

usethis::use_data(code, overwrite = TRUE)

##The TSS code and study name data is from https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tissue-source-site-codes
##The Study Abbreviations data is from https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations

code <- read.csv("~/NeoEnrichment/data-raw/TCGA_code.csv",stringsAsFactors = FALSE) %>%
  mutate(TSS.Code=ifelse(is.na(TSS.Code),"NA",TSS.Code))
study <- read.csv("~/NeoEnrichment/data-raw/study.csv")
code <- left_join(code,
                  study,
                  by="Study.Name")
cancer_type_code <- code
cancer_type_code <- cancer_type_code %>%
  rename(cancer_code=TSS.Code,cancer_type=Study.Abbreviation,
         cancer_type_full_name=Study.Name) %>%
  select(-Source.Site,-BCR)
cancer_type_code[,1:3] <- apply(cancer_type_code[,1:3],2,as.character)
usethis::use_data(cancer_type_code,internal = TRUE,overwrite = TRUE)

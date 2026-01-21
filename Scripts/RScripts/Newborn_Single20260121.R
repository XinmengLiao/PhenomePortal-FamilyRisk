## ---- Load libraries ----
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(xlsx)

for (i in c("select","filter", "mutate","rename", "left_join", "slice")){
  conflicted::conflict_prefer(i, "dplyr")
}
rm(i)
conflicted::conflicts_prefer(stats::sd)
conflicted::conflicts_prefer(httr::content)
conflicted::conflicts_prefer(plotly::layout)

## ---- Define file paths ----
location = "server" # server or local 

args <- commandArgs(trailingOnly = TRUE)
sampleID <- args[1]
result_file <- args[2]
gender <- args[3]
output_dir <- args[4]
genelist <- args[5]

if (location == "local"){
  compare_file <- '/Users/xinmengliao/Documents/Project/20250710_NewbornRisk/Datasets/NBSeq_Results.xlsx'
  TR_removed_variant <- '/Users/xinmengliao/Documents/Project/20250710_NewbornRisk/Datasets/TRpipelineRemovedVariants.txt'
  genedb_file <- "/Users/xinmengliao/Documents/Project/20250710_NewbornRisk/Datasets/genelists/Preset_screening_list_GenCC_20251125.txt"
}

if (location == "server"){
  compare_file <- '/mnt/nas/Genomics/Genome/FamilyRisk/Datasets/NBSeq_Results.xlsx'
  TR_removed_variant <- '/mnt/nas/Genomics/Genome/FamilyRisk/Datasets/TRpipelineRemovedVariants.txt'
  genedb_file <- "/mnt/nas/Genomics/Genome/FamilyRisk/Datasets/Preset_screening_list_GenCC_20251125.txt"
}


## ---- Load files ----
result <- read.csv(result_file,header = T,sep = "\t") 
colnames(result) <- gsub("am_","AlphaMissense_", colnames(result))
genelist <- unlist(strsplit(genelist,","))

genedb <- read.csv(genedb_file,header = T,sep = "\t") %>% filter(Project %in% genelist)
genedb1 <- genedb %>% select(Genes, MIM, Inheritance, Disorder_Group,GenCC_Classification, GenCC_Submitter ) %>% unique()
result1 <- result %>% 
  mutate(MIM = as.character(MIM)) %>% 
  left_join(., genedb1, by = c("Genes","MIM", "Inheritance")) %>% unique()
  
if ("NBScreening" %in% genelist){
  remove.variant <- read.csv(TR_removed_variant, header = T, sep = "\t")
  result <- result %>% filter(!(variant_info %in% remove.variant$variant_info))
}

compare_positive <- read.xlsx(compare_file,sheetName = "All-Positive") %>% filter(!is.na(Reported.Disease)) %>% select(1:9) %>% 
  mutate(Allele.Count = round(Count),
         new_count = paste0(Count.Type, Count)) %>% 
  mutate(result_tmp = paste0(Project, " ", new_count, " Allele Frequency:",Allele.Freuqency)) %>% 
  group_by(Reported.Disease, Genes, Transcript, Protein) %>% 
  mutate(`Results from Other NBSeq Projects` = paste(result_tmp, collapse = "; ")) %>% ungroup()

compare_carrier <- read.xlsx(compare_file,sheetName = "All-Carrier") %>% filter(!is.na(Reported.Disease)) %>% select(1:8) %>% 
  mutate(Allele.Count = round(Allele.Count)) %>% 
  mutate(result_tmp = paste0(Project, " Allele Count:", Allele.Count, " Allele Frequency:", Allele.Frequency)) %>%
  group_by(Reported.Disease, Genes, Transcript, Protein) %>% 
  mutate(`Results from Other NBSeq Projects` = paste(result_tmp, collapse = "; ")) %>% ungroup()

## ---- General Summary Table ----
total_variantnum <- length(unique(result$variant_info))
total_plp <- result %>% filter(grepl("Pathogenic|Likely_pathogneic",ClinVar_CLNSIG) |
                                 grepl("Pathogenic|Likely_pathogneic",acmg_classification)) %>% 
  pull(variant_info) %>% unique() %>% length()
general.info = data.frame(
  Content = c("Sample ID", "Sex", "Filtered variants", "Clinical P/LP variants"),
  Count = c(sampleID, gender, total_variantnum, total_plp)
)


## ---- Functions ----
genotype.count.single <- function(data, sampleid, gender){
  data <- data %>%
    mutate(Gender = gender) %>% 
    mutate(Genotype = sapply(strsplit(.data[[sampleid]], ":"), `[`, 1)) %>% 
    mutate(Zygosity = case_when(Genotype %in% c("1/0", "0/1", "0|1", "1|0") ~ "Heterozygous",
                                Genotype %in% c("1/1", "1|1") ~ "Homozygous",
                                 TRUE ~ "Other")) %>%
    mutate(Zygosity = case_when(
      (grepl("X", X.CHROM) & Zygosity == "Heterozygous" & gender == "Male") ~ "Hemizygous",
      TRUE ~ Zygosity
      )) %>% 
    mutate(HGVSc = gsub("^.*:","",HGVSc),
           HGVSp = gsub("^.*:","",HGVSp),
           ClinVar_CLNSIG = sapply(strsplit(split = "&", ClinVar_CLNSIG),`[`,1))
  
  return(data)
}

positive.monogenic.plp.single <- function(data){
  plp <- data %>% filter(grepl("Pathogenic|Likely_pathogneic",ClinVar_CLNSIG) &
                           grepl("Pathogenic|Likely_pathogneic",acmg_classification))
  
  # positive monogenic disease
  plp.monogenic <- plp %>% filter(
    (Inheritance %in% c("AD","Multiple and/or complex pattern","AR; AD","AD; AR","X-linked multiple and/or complex pattern","XLD") & 
       Zygosity %in% c("Heterozygous", "Hemizygous","Homozygous")) |
    (Inheritance %in% c("AR","XLR") & Zygosity %in% c("Hemizygous","Homozygous"))) %>% unique()
  
  plp.monogenic.final <- plp.monogenic %>% 
    mutate(Disease = paste(unique(Disease), collapse = "; ")) %>% 
    select(Disease, Inheritance, Genes, X.CHROM, POS, REF, ALT, Zygosity, HGVSc, HGVSp, ClinVar_CLNSIG, acmg_classification, 
           IMPACT, MAX_AF) %>% unique() %>% 
    left_join(., compare_positive %>% 
                select(Reported.Disease, `Results from Other NBSeq Projects`, Genes, Transcript) %>% unique(),
              by = c("Genes", "HGVSc" = "Transcript")) %>% 
    rename(`Reported Disease from Other NBSeq Projects` = Reported.Disease) %>% unique()
  
  return(plp.monogenic.final)
}

carrier.plp <- function(data){
  carrier <- data %>% filter(grepl("Pathogenic|Likely_pathogneic",ClinVar_CLNSIG) &
                             grepl("Pathogenic|Likely_pathogneic",acmg_classification)) %>% 
    filter(Inheritance %in% c("AR","XLR")) 
  
  carrier.disease <- carrier %>% 
    filter((Inheritance %in% c("AR","XLR") & Zygosity == "Heterozygous")) %>% unique() %>% 
    select(Disease, Inheritance, Genes, X.CHROM, POS, REF, ALT, Zygosity, HGVSc, HGVSp, ClinVar_CLNSIG, acmg_classification, 
           IMPACT, MAX_AF) %>% unique() %>% 
    left_join(., compare_carrier %>% 
                select(Reported.Disease, `Results from Other NBSeq Projects`, Genes, Transcript) %>% unique(),
              by = c("Genes", "HGVSc" = "Transcript")) %>% 
    rename(`Reported Disease from Other NBSeq Projects` = Reported.Disease) %>% unique()
  
  return(carrier.disease)
}

## ---- Generate screening results ----
result1 <- genotype.count.single(result, sampleid = sampleID, gender = gender)
plp.positive <- positive.monogenic.plp.single(result1)

# carrier status 
if("NBScreening" %in% genelist){
  result_carrier <- result1 %>% 
    filter(Genes != "MEFV" & HGVSp != "p.Met694Val") %>% 
    filter(Genes != "MEFV" & HGVSp != "p.Met680Ile")
}else{
  result_carrier <- result1
}
carrier.status <- carrier.plp(result_carrier)

## ---- Writing output tables ----
# Sample and variant summary 
write.table(general.info, paste0(output_dir, "/Results/general_summary.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

# scren-positive monogenic disease and carrier status 
if(nrow(plp.positive) == 0){
  plp.positive <- data.frame(Message = "No positive monogenic disease variant detected based on the current gene panel and filtering criteria.")
}
if(nrow(carrier.status) == 0){
  carrier.status <- data.frame(Message = "No positive monogenic disease carrier status variant detected based on the current gene panel and filtering criteria.")
}
write.table(plp.positive, paste0(output_dir, "/Results/monogenic_positive_results.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(carrier.status, paste0(output_dir, "/Results/carrier_status_results.txt"), sep = "\t", row.names = FALSE, quote = FALSE)


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
args <- commandArgs(trailingOnly = TRUE)
sampleID <- args[1]
result_file <- args[2]
metadata_file <- args[3]
output_dir <- args[4]
genelist <- args[5]

result_file <- "/Users/xinmengliao/Documents/Project/20250710_NewbornRisk/Batch/newborn103_vep_merged_rmmissingalt_biallelic2.txt"
genedb_file <- "/Users/xinmengliao/Documents/Project/20250710_NewbornRisk/Datasets/genelists/Preset_screening_list_GenCC_20251125.txt"
gender_file <- "/Users/xinmengliao/Documents/Project/20250710_NewbornRisk/Batch/newborn103gender.txt"
output_dir <- '/Users/xinmengliao/Documents/Project/20250710_NewbornRisk/Results20251124/'
genelist <- "TR"

compare_file <- '/Users/xinmengliao/Documents/Project/20250710_NewbornRisk/Datasets/NBSeq_Results.xlsx'
TR_removed_variant <- '/Users/xinmengliao/Documents/Project/20250710_NewbornRisk/Scripts/MainFunctions/batch/TRpipelineRemovedVariants.txt'

## ---- Load files ----
result <- read.csv(result_file,header = T,sep = "\t") 
colnames(result) <- gsub("am_","AlphaMissense_", colnames(result))
genelist <- unlist(strsplit(genelist,","))

if ("TR" %in% genelist){
  remove.variant <- read.csv(TR_removed_variant, header = T, sep = "\t")
  result <- result %>% filter(!(variant_info %in% remove.variant$variant_info))
}

# genedb <- read.csv(genedb_file,header = T,sep = "\t") %>% filter(Project == "NBScreening")
# genedb1 <- genedb %>% select(Genes, MIM, Inheritance, Disorder_Group,GenCC_Classification, GenCC_Submitter ) %>% unique()
# result <- result %>% left_join(., genedb1, by = c("Genes","MIM", "Inheritance")) %>% unique()
# which(is.na(result$Disorder_Group))

gender <- read.csv(gender_file,header = F,sep = "\t") # include sample ID and gender 
colnames(gender) <- c("Sample","Gender")

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
total_samplenum <- nrow(gender)
total_variantnum <- length(unique(result$variant_info))
total_plp <- result %>% filter(grepl("Pathogenic|Likely_pathogneic",ClinVar_CLNSIG) |
                                 grepl("Pathogenic|Likely_pathogneic",acmg_classification)) %>% 
  pull(variant_info) %>% unique() %>% length()

## ---- Functions ----
genotype.count <- function(data, sampleid, gender = F){
  #data = plp
  # clean the genotype, only keeps the genotyps information
  data[, sampleid] <- lapply(data[, sampleid], function(col) {
    sapply(strsplit(col, ":"), `[`, 1)
  })
  
  # dominant count 
  data$homo_count <- apply(data[, sampleid], 1, function(row) {
    sum(row %in% c("1/1", "1|1"), na.rm = TRUE)
  })
  
  data$homo_sample <- apply(data[, sampleid], 1, function(row) {
    paste(names(row)[row %in% c("1/1", "1|1")], collapse = ",")
  })
  
  data <- data %>% separate_rows(homo_sample, sep = ",")
  
  # recessive count 
  data$het_count <- apply(data[, sampleid], 1, function(row) {
    sum(row %in% c("1/0", "1|0","0/1","0|1"), na.rm = TRUE)
  })
  
  data$het_sample <- apply(data[, sampleid], 1, function(row) {
    paste(names(row)[row %in% c("1/0", "1|0", "0/1", "0|1")], collapse = ",")
  })
  
  data <- data %>% separate_rows(het_sample, sep = ",")
  
  # Extract homozygous genotype per sample
  data$homo_genotype <- mapply(function(row, sample) {
    if (sample %in% colnames(data)) {
      data[row, sample]
    } else {
      NA
    }
  }, seq_len(nrow(data)), data$homo_sample)
  
  # Extract heterozygous genotype per sample
  data$het_genotype <- mapply(function(row, sample) {
    if (sample %in% colnames(data)) {
      data[row, sample]
    } else {
      NA
    }
  }, seq_len(nrow(data)), data$het_sample)
  
  
  # add zygosity with matching genders 
  if (!identical(gender, FALSE)){
    data <- data %>% left_join(., gender, by = c("homo_sample" = "Sample")) %>% 
      rename(homo_gender = Gender) %>% 
      left_join(., gender, by = c("het_sample" = "Sample")) %>% 
      rename(het_gender = Gender) %>% unique() 
  }
  
  data <- data %>% 
    mutate(het_zygosity = if_else(!is.na(het_gender), case_when((het_gender == "Male" & X.CHROM == "chrX" & 
                                                                   Inheritance %in% c("AR", "XLR")) ~ "Hemizygous",T ~ "Heterozygous"), NA),
           homo_zygosity = if_else(!is.na(homo_gender), "Homozygous", NA) ) %>% unique()
  
  index <- which(data$het_zygosity == "Hemizygous")
  data$homo_zygosity[index] = "Hemizygous"
  data$homo_sample[index] = data$het_sample[index]
  data$het_sample[index] = ""
  data$homo_count[index] = data$homo_count[index] + 1
  data$het_count[index] = data$het_count[index] - 1
  
  return(data)
}

positive.monogenic.plp <- function(data){
  plp <- data %>% filter(grepl("Pathogenic|Likely_pathogneic",ClinVar_CLNSIG) &
                           grepl("Pathogenic|Likely_pathogneic",acmg_classification))
  
  plp <- genotype.count(plp, sampleid = gender$Sample, gender = gender)
  
  # manage plp.df colnames 
  plp <- plp %>% mutate(HGVSc = gsub("^.*:","",HGVSc),
                        HGVSp = gsub("^.*:","",HGVSp),
                        ClinVar_CLNSIG = sapply(strsplit(split = "&", ClinVar_CLNSIG),`[`,1))
  
  
  # positive monogenic disease
  plp.monogenic <- plp %>% filter((Inheritance %in% c("AD","Multiple and/or complex pattern","AR; AD","AD; AR",
                                                      "X-linked multiple and/or complex pattern","XLD") & het_count > 0) |
                                    (Inheritance %in% c("AD","Multiple and/or complex pattern","AR; AD","AD; AR",
                                                        "X-linked multiple and/or complex pattern","XLD") & homo_count > 0) | 
                                    (homo_zygosity == "Hemizygous" | het_zygosity == "Hemizygous") |
                                    (Inheritance %in% c("AR","XLR") & homo_count > 0)) %>% 
    mutate(across(c(het_sample, homo_sample), as.character)) 
  
  plp.disease <- plp.monogenic %>% 
    group_by(Genes) %>% 
    mutate(Disease = paste(unique(Disease), collapse = "; ")) %>% ungroup() %>% unique() %>% 
    pivot_longer(cols = c(het_sample, homo_sample), values_to = "sample_id") %>% 
    filter(sample_id != "" & !is.na(sample_id)) %>% 
    filter(!is.na(Disease) & Disease != "") %>% 
    group_by(Disease) %>% summarise(disease_count = n_distinct(sample_id), .groups = "drop") %>% 
    rename(`No. of sample for disease` = disease_count) %>% 
    mutate(`Disease carreir frequency` = `No. of sample for disease` / total_samplenum)
  
  variant.freq <- plp.monogenic %>% 
    group_by(Genes) %>% 
    mutate(Disease = paste(unique(Disease), collapse = "; ")) %>% ungroup() %>% unique() %>% 
    pivot_longer(cols = c(het_sample, homo_sample), values_to = "sample_id") %>% 
    filter(sample_id != "" & !is.na(sample_id)) %>% 
    filter(!is.na(Disease) & Disease != "") %>% 
    select(variant_info, Genes, Disease,sample_id) %>% unique() %>% 
    group_by(across(-sample_id)) %>% 
    summarise(`No. of sample for variant` = n(), .groups = "drop") %>% 
    mutate(`Variant carrier frequency` = `No. of sample for variant` / total_samplenum )
  
  plp.monogenic.gene <- plp.monogenic %>%
    pivot_longer(cols = c(het_sample, homo_sample), names_to = "zygosity_type", values_to = "sample_id") %>% 
    filter(sample_id != "" & !is.na(sample_id)) %>% 
    filter(!is.na(Disease) & Disease != "") %>% 
    select(variant_info, Genes, Disease,sample_id) %>% unique() %>% 
    group_by(Genes) %>%
    summarise(gene_count = n_distinct(sample_id), .groups = "drop") %>% 
    arrange(desc(gene_count)) %>% 
    mutate(proportion = gene_count / total_samplenum) %>% unique()
  
  plp.monogenic.gene1 <- plp.monogenic.gene %>% 
    left_join(., plp.monogenic, by = "Genes") %>% 
    pivot_longer(cols = c(homo_zygosity, het_zygosity),values_to = "Zygosity") %>% 
    filter(!is.na(Zygosity) & Zygosity != "") %>% 
    select(Disease, Inheritance, Genes, X.CHROM, POS, REF, ALT, Zygosity, HGVSc, HGVSp, ClinVar_CLNSIG, acmg_classification, 
           IMPACT, MAX_AF, gene_count, proportion) %>% 
    rename(Inh = Inheritance, Chr = X.CHROM, Pos = POS, Ref = REF, Alt = ALT, `Global Max AF` = MAX_AF,  
           `No. of sample for gene` = gene_count, `Gene carreir frequency` = proportion) %>% unique() %>% 
    group_by(Genes) %>% 
    mutate(Disease = paste(unique(Disease), collapse = "; ")) %>% ungroup() %>% unique() 
  
  # merge other frequencies together 
  plp.monogenic.gene.final <- plp.monogenic.gene1 %>% 
    left_join(., plp.disease, by = "Disease") %>% 
    mutate(variant_info = paste(Chr, Pos, Ref, Alt, sep = "_")) %>% 
    left_join(., variant.freq) %>% 
    select(Disease, Inh, Genes, Chr, Pos, Ref, Alt, Zygosity, HGVSc, HGVSp, ClinVar_CLNSIG, acmg_classification, 
           IMPACT, `Global Max AF`, `No. of sample for variant`, `Variant carrier frequency`,
           `No. of sample for gene`, `Gene carreir frequency`, 
           `No. of sample for disease`, `Disease carreir frequency`, ) %>% 
    left_join(., compare_positive %>% 
                select(Reported.Disease, `Results from Other NBSeq Projects`, Genes, Transcript) %>% 
                unique(), by = c("Genes", "HGVSc" = "Transcript")) %>% unique() %>% 
    rename(`Reported Disease from Other NBSeq Projects` = Reported.Disease) %>% 
    arrange(
      desc(`No. of sample for disease`),
      desc(`No. of sample for gene`),
      desc(`No. of sample for variant`)
    ) 
  
  # statistics of the variant per newborn
  plp.monogenic.stat <- plp.monogenic %>% 
    pivot_longer(cols = c(het_sample, homo_sample), names_to = "zygosity_type", values_to = "sample_id") %>% 
    filter(sample_id != "" & !is.na(sample_id)) %>% 
    filter(!is.na(Disease) & Disease != "") %>% unique() %>% 
    select(variant_info, sample_id) %>% unique() %>% 
    group_by(sample_id) %>% summarise(count = n()) %>% ungroup() %>% 
    group_by(count) %>% summarise(nsample = n_distinct(sample_id))
    
  
  return(list(monogenic.positive = plp.monogenic.gene.final, monogenic.variant = plp.monogenic$variant_info, 
              stat_data = plp.monogenic.stat))
}

carrier.plp <- function(data){
  carrier <- data %>% filter(grepl("Pathogenic|Likely_pathogneic",ClinVar_CLNSIG) &
                               grepl("Pathogenic|Likely_pathogneic",acmg_classification)) %>% 
    filter(Inheritance %in% c("AR","XLR")) 
  
  carrier <- genotype.count(carrier, sampleid = gender$Sample, gender = gender) %>% 
    filter(!(X.CHROM == "chrX" & het_zygosity == "Hemizygous")) %>% 
    mutate(HGVSc = gsub("^.*:","",HGVSc),
           HGVSp = gsub("^.*:","",HGVSp),
           ClinVar_CLNSIG = sapply(strsplit(split = "&", ClinVar_CLNSIG),`[`,1))
  
  carrier.status <- carrier %>% filter((Inheritance %in% c("AR","XLR") & het_count > 0)) %>% 
    mutate(across(c(het_sample, homo_sample), as.character)) 
  
  carrier.disease <- carrier.status %>% 
    group_by(Genes) %>% 
    mutate(Disease = paste(unique(Disease), collapse = "; ")) %>% ungroup() %>% unique() %>% 
    pivot_longer(cols = c(het_sample, homo_sample), values_to = "sample_id") %>% 
    filter(sample_id != "" & !is.na(sample_id)) %>% 
    filter(!is.na(Disease) & Disease != "") %>% 
    group_by(Disease) %>% summarise(disease_count = n_distinct(sample_id), .groups = "drop") %>% 
    rename(`No. of sample for disease` = disease_count) %>% 
    mutate(`Disease carreir frequency` = `No. of sample for disease` / total_samplenum)
  
  variant.freq <- carrier.status %>% 
    group_by(Genes) %>% 
    mutate(Disease = paste(unique(Disease), collapse = "; ")) %>% ungroup() %>% unique() %>% 
    pivot_longer(cols = c(het_sample, homo_sample), values_to = "sample_id") %>% 
    filter(sample_id != "" & !is.na(sample_id)) %>% 
    filter(!is.na(Disease) & Disease != "") %>% 
    select(variant_info, Genes, Disease,sample_id) %>% unique() %>% 
    group_by(variant_info) %>% 
    summarise(`No. of sample for variant` = n(), .groups = "drop") %>% 
    mutate(`Variant carrier frequency` = `No. of sample for variant` / total_samplenum )
  
  carrier.status.gene <- carrier.status %>%
    pivot_longer(cols = c(het_sample, homo_sample), names_to = "zygosity_type", values_to = "sample_id") %>% 
    filter(sample_id != "" & !is.na(sample_id)) %>% 
    group_by(Genes) %>%
    summarise(gene_count = n_distinct(sample_id), .groups = "drop") %>% 
    arrange(desc(gene_count)) %>% 
    mutate(proportion = gene_count / total_samplenum) %>% unique()
  
  carrier.status.gene1 <- carrier.status.gene %>% 
    left_join(., carrier.status, by = "Genes") %>% 
    pivot_longer(cols = c(homo_zygosity, het_zygosity),values_to = "Zygosity") %>% 
    filter(!is.na(Zygosity) & Zygosity != "") %>% 
    select(Disease, Inheritance, Genes, X.CHROM, POS, REF, ALT, Zygosity, HGVSc, HGVSp, ClinVar_CLNSIG, acmg_classification, 
           IMPACT, MAX_AF, gene_count, proportion) %>% 
    rename(Inh = Inheritance, Chr = X.CHROM, Pos = POS, Ref = REF, Alt = ALT, `Global Max AF` = MAX_AF,  
           `No. of sample for gene` = gene_count, `Gene carreir frequency` = proportion) %>% unique() %>% 
    group_by(Genes) %>% 
    mutate(Disease = paste(unique(Disease), collapse = "; ")) %>% ungroup() %>% unique() 
  
  # merge other frequencies together 
  carrier.status.gene.final <- carrier.status.gene1 %>% 
    left_join(., carrier.disease, by = "Disease") %>% 
    mutate(variant_info = paste(Chr, Pos, Ref, Alt, sep = "_")) %>% 
    left_join(., variant.freq) %>% 
    select(Disease, Inh, Genes, Chr, Pos, Ref, Alt, Zygosity, HGVSc, HGVSp, ClinVar_CLNSIG, acmg_classification, 
           IMPACT, `Global Max AF`, `No. of sample for variant`, `Variant carrier frequency`,
           `No. of sample for gene`, `Gene carreir frequency`, 
           `No. of sample for disease`, `Disease carreir frequency`, ) %>% 
    left_join(., compare_carrier %>% 
                select(Reported.Disease, `Results from Other NBSeq Projects`, Genes, Transcript) %>% 
                unique(), by = c("Genes", "HGVSc" = "Transcript")) %>% unique() %>% 
    rename(`Reported Disease from Other NBSeq Projects` = Reported.Disease) %>% 
    arrange(
      desc(`No. of sample for disease`),
      desc(`No. of sample for gene`),
      desc(`No. of sample for variant`)
    ) 
  
  # statistics of the variant per newborn
  carrier.status.stat <- carrier.status %>% 
    pivot_longer(cols = c(het_sample, homo_sample), names_to = "zygosity_type", values_to = "sample_id") %>% 
    filter(sample_id != "" & !is.na(sample_id)) %>% 
    filter(!is.na(Disease) & Disease != "") %>% unique() %>% 
    select(variant_info, sample_id) %>% unique() %>% 
    group_by(sample_id) %>% summarise(count = n()) %>% ungroup() %>% 
    group_by(count) %>% summarise(nsample = n_distinct(sample_id))
  
  return(list(carrier_final = carrier.status.gene.final, carrier_stat = carrier.status.stat))
}

## ---- Generate screening results ----
plp.positive <- positive.monogenic.plp(result)$monogenic.positive
plp.variant <- positive.monogenic.plp(result)$monogenic.variant

# carrier status 
if("TR" %in% genelist){
  result_carrier <- result %>% 
    filter(Genes != "MEFV" & HGVSp != "p.Met694Val") %>% 
    filter(Genes != "MEFV" & HGVSp != "p.Met680Ile")
}
carrier.status <- carrier.plp(result_carrier)$carrier_final


## ---- Statistics of positive variant each newborn carry (Barplot)----
stat.carrier <- carrier.plp(result_carrier)$carrier_stat %>% mutate(Type = "Carrier status")
stat.monogenic =  positive.monogenic.plp(result)$stat_data %>% mutate(Type = "Monogenic Diseases")
stat.all <- rbind(stat.carrier,stat.monogenic) %>% 
  rename(Var_num = count, Sample_num = nsample)

## ---- Statistics of Disease Ontology ----
on.child <- plp.positive %>% 
  mutate(variant_info = paste(Chr, Pos, Ref, Alt, sep = "_")) %>% 
  separate_rows(sep = "; ",Disease) %>% 
  left_join(., genedb %>% select(Genes, Disease, Inheritance, Disorder_Group), by = c("Genes","Inh" = "Inheritance","Disease")) %>% 
  mutate(Disorder_Group = ifelse(grepl(";", Disorder_Group), "Complex genetic syndromes", Disorder_Group)) %>% 
  select(Disease, Disorder_Group) %>% unique() %>% 
  group_by(Disorder_Group) %>% 
  summarise(count = n()) %>% ungroup() %>% mutate(Type = "Positives")

on.carrier<- carrier.status %>% 
  mutate(variant_info = paste(Chr, Pos, Ref, Alt, sep = "_")) %>% 
  separate_rows(sep = "; ",Disease) %>% 
  left_join(., genedb %>% select(Genes, Disease, Inheritance, Disorder_Group), by = c("Genes","Inh" = "Inheritance","Disease")) %>% 
  select(Disease, Disorder_Group) %>% unique() %>% 
  mutate(Disorder_Group = if_else(grepl(";",Disorder_Group),"Complex genetic syndromes",Disorder_Group)) %>% unique() %>% 
  group_by(Disorder_Group) %>% 
  summarise(count = n()) %>% ungroup() %>% mutate(Type = "Carrier") 

on.all <- rbind(on.child, on.carrier)

## ---- Inhosue allele frequency ----
sample.col <- colnames(result)[which(colnames(result) %in% gender$Sample)]
female.num <- table(gender$Gender) %>% as.data.frame() %>% filter(Var1 == "Female") %>% pull(Freq) %>% as.numeric()
male.num <- table(gender$Gender) %>% as.data.frame() %>% filter(Var1 == "Male") %>% pull(Freq) %>% as.numeric()
normal.count <- (female.num + male.num) * 2
x.count <- female.num * 2 + male.num
y.count <- male.num

inhouse.af <- result %>% 
  pivot_longer(cols = all_of(sample.col), values_to = "genotype",names_to = "Sample") %>% 
  filter(!grepl("\\./\\.", genotype)) %>% 
  #left_join(., gender) %>% 
  mutate(genotype = substr(genotype, 1, 3),
         a1 = as.numeric(sapply(strsplit(split = "/|\\|",genotype), `[`,1)),
         a2 = as.numeric(sapply(strsplit(split = "/|\\|",genotype), `[`,2)),
         cohort_AC = a1+a2) %>% 
  select(variant_info, genotype, a1, a2, cohort_AC, Sample) %>% unique() %>% 
  select(variant_info,cohort_AC, Sample) %>% 
  group_by(variant_info) %>% 
  summarise(sum_cohort_AC = sum(cohort_AC)) %>% ungroup() %>% 
  mutate(cohort_AF = case_when(grepl("X", variant_info) ~ sum_cohort_AC / x.count,
                               grepl("Y", variant_info) ~ sum_cohort_AC / y.count,
                                TRUE ~ sum_cohort_AC / normal.count)) %>% 
  rename(cohort_AC = sum_cohort_AC)

## ---- Writing output tables ----
# Cohort and variant summary 
general.info = data.frame(
  Content = c("Total samples", "Filtered variants", "Clinical P/LP variants"),
  Count = c(total_samplenum, total_variantnum, total_plp)
)
write.table(general.info, paste0(output_dir, "/general_summary.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

# scren-positive monogenic disease and carrier status 
if(nrow(plp.positive) == 0){
  plp.positive <- data.frame(Message = "No positive monogenic disease variant detected based on the current gene panel and filtering criteria.")
}
if(nrow(carrier.status) == 0){
  plp.positive <- data.frame(Message = "No positive monogenic disease variant detected based on the current gene panel and filtering criteria.")
}
write.table(plp.positive, paste0(output_dir, "/monogenic_positive_results.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(carrier.status, paste0(output_dir, "/carrier_status_results.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

# Statistics of positive variant each newborn carry (Barplot)
write.table(stat.all, paste0(output_dir, "/plp_var_stat.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

# Statistics of disease ontology (Barplot)
write.table(on.all, paste0(output_dir, "/ontology_stat.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

# Inhouse allele frequency
write.table(inhouse.af, paste0(output_dir, "/inhouse_allele_frequency.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

## ---- Load libraries ----
library(kinship2)
library(tidyr)
library(dplyr)
for (i in c("select","filter", "mutate","rename", "left_join", "slice")){
  conflicted::conflict_prefer(i, "dplyr")
}
rm(i)
conflicted::conflicts_prefer(stats::sd)
conflicted::conflicts_prefer(httr::content)
conflicted::conflicts_prefer(plotly::layout)

## ---- Define file paths ----
location = "local" # server or local 

args <- commandArgs(trailingOnly = TRUE)
sampleid <- args[1]
input_file <- args[2]
ped_file <- args[3]
output_dir <-  args[4]
genelist <- args[5]


sampleid = "rwgsF1"
input_file = "/Users/xinmengliao/Desktop/rwgsF1.txt"
ped_file = "/Users/xinmengliao/Desktop/rwgs_F1.ped"
output_dir = "/Users/xinmengliao/Desktop"
genelist <- c("NBScreening")

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

print(paste0("Now pocessing ", sampleid))
print(paste0("Input file: ", input_file))
print(paste0("Ped file: ", ped_file))
print(paste0("Output directory: ", output_dir))

## ---- Load files ----
result <- read.csv(input_file,header = T,sep = "\t") %>% unique() %>% 
  mutate(ClinVar_CLNSIG = gsub("_", " ", ClinVar_CLNSIG),
         acmg_classification = gsub("_", " ", acmg_classification))
colnames(result) <- gsub("am_","AlphaMissense_", colnames(result))
if ("NBScreening" %in% genelist){
  remove.variant <- read.csv(TR_removed_variant, header = T, sep = "\t")
  result <- result %>% filter(!(variant_info %in% remove.variant$variant_info))
}

ped <- read.csv(ped_file, header = F,sep = "\t") %>% unique()
colnames(ped) <- c("FID","Individual","Father","Mother","Sex","Genotype")
familyid <- ped$FID %>% unique()

## ---- Family and Variant Summary ----
family_info <- ped %>% 
  mutate(Relationship = case_when(Sex == 1 ~ "Father", Sex == 2 ~ "Mother", TRUE ~ NA)) %>% 
  select(Individual, Relationship) %>% rename(SampleID = Individual)

plp.count <- result %>% 
  filter(grepl("Pathogenic|Likely_pathogenic",ClinVar_CLNSIG) | grepl("Pathogenic|Likely_pathogenic",acmg_classification)) %>% 
  select(variant_info) %>% unique()

var_summary <- data.frame(Content = c("Filtered variants", "Clinical P/LP variants"),
                          Count = c(length(unique(result$variant_info)), 
                                    nrow(plp.count)))

## ---- 2. Functions ----

# Function 1: filter out the compound heterozygous for recessive disease and add back to the df.
compound.heterozygous.res <- function(df){
    
	brief.tmp <- df %>%
		filter(Inheritance %in% c("AR", "XLR")) %>%
    mutate(ClinVar_CLNSIG = gsub("_", " ", ClinVar_CLNSIG),
        acmg_classification = gsub("_", " ", acmg_classification),
        ClinVar_CLNSIG = sapply(strsplit(split = "&", ClinVar_CLNSIG), `[`, 1),
        ClinVar_CLNSIG = if_else(ClinVar_CLNSIG == "Pathogenic/Likely pathogenic/Pathogenic",
                                 "Pathogenic/Likely pathogenic",
                                 as.character(ClinVar_CLNSIG)),
        Existing_variation = sapply(strsplit(split = "&", Existing_variation), `[`, 1)) %>% 
	  select(Disease, Genes, Inheritance, variant_info, Existing_variation, ClinVar_CLNSIG, acmg_classification, FatherGenotype, MotherGenotype) %>%
	  rename(`Existing variant` = Existing_variation) %>% distinct()
	
	
	# function for infer the compound heterozygous origin
	infer_origin <- function(father, mother) {
	  # if (length(father) == 0 || length(mother) == 0) {
	  #   return("Error, please check again")
	  # }
	  # if (is.na(father) || is.na(mother)) {
	  #   return("Error, please check again")
	  # }
	  if (father %in% c("0/1", "1/0", "1|1", "1/1", "1|0", "0|1") & mother %in% c("0/0","0|0")){
	    return("Father")
	  }else if (father %in% c("1/1","1|1") & mother %in% c("0/1", "1/0")){
	    return("Father")
	  }else if (mother %in% c("0/1", "1/0", "1|1", "1/1", "1|0", "0|1") & father %in% c("0/0","0|0")){
	    return("Mother")
	  }else if (mother %in% c("1/1","1|1") & father %in% c("0/1", "1/0")){
	    return("Mother")
	  }else if (father == "1|0" & mother  == "1|0"){
	    return("Father")
	  }else if (father == "0|1" & mother  == "0|1"){
	    return("Mother")
	  }else if (father %in% c("./.",".|.")){
	    return("Mother")
	  }else if (mother %in% c("./.",".|.")) {
	    return("Father")
	  }else if (father == "1" & mother == "."){
	    return("Uncertained chrX variant")
	  }else{
	    return("Uncertained")
	  }
	}
	
  brief.com.het1 <- brief.tmp %>%
    select(Genes, variant_info, ClinVar_CLNSIG, acmg_classification,FatherGenotype, MotherGenotype) %>% 
    rowwise() %>%
    mutate(origin = infer_origin(FatherGenotype, MotherGenotype)) %>%
    ungroup() %>% 
    filter(ClinVar_CLNSIG %in% c("Pathogenic/Likely pathogenic","Pathogenic","Likely pathogenic","Uncertain significance") &
             acmg_classification %in% c("Pathogenic","Likely pathogenic","Uncertain significance")) %>% unique()
    
  ch.df <- brief.com.het1 %>% group_by(Genes) %>%
    filter(n() >= 2) %>%   # must have more than two variant
    summarise(
      origins = list(unique(origin)),
      Compound_heterozygous = case_when(
        all(c("Father","Mother") %in% origins[[1]]) ~ "Yes",
        all(origins[[1]] == "Father") ~ "No",
        all(origins[[1]] == "Mother") ~ "No",
        TRUE ~ "Uncertained"
      ),
      .groups = "drop"
    ) 
  
  brief.com.het2 <- brief.com.het1 %>% 
    left_join(., ch.df) %>% select(-origins) %>% filter(Compound_heterozygous %in% c("Yes","Uncertained")) 
  
  df1 <- df %>% 
    left_join(., brief.com.het2 %>% select(-ClinVar_CLNSIG, -acmg_classification)) %>% unique() 

	return(df1)
}

# Function 2: Potential recessive disease on the offspring 
carrier.couple <- function(df){
  df.ch <- compound.heterozygous.res(df)
  plp.carrier <- df.ch %>% 
    filter(grepl("Pathogenic|Likely pathogenic",ClinVar_CLNSIG) & 
           grepl("Pathogenic|Likely pathogenic",acmg_classification) |
           !is.na(Compound_heterozygous)) %>%  
    filter(Inheritance %in% c("AR", "XLR")) %>% unique()
  
  if(nrow(plp.carrier) == 0){
    return(plp.carrier)
  }else{
    plp.carrier.final <- plp.carrier %>% 
      mutate(HGVSc = gsub("^.*:","",HGVSc), HGVSp = gsub("^.*:","",HGVSp)) %>% 
      select(Disease, Inheritance, Genes, X.CHROM, POS, REF, ALT, HGVSc, HGVSp, ClinVar_CLNSIG, acmg_classification, 
             IMPACT, MAX_AF, FatherGenotype, MotherGenotype, Compound_heterozygous) %>% unique() 
    return(plp.carrier.final)
  }

}


## ---- Generate screening results ----
carrier.res.couple <-  carrier.couple(df = result)

## ---- Writing output tables ----
# Family and variant summary 
write.table(family_info, paste0(output_dir, "/Results/family_information.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(var_summary, paste0(output_dir, "/Results/variant_summary.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

# scren-positive monogenic disease and carrier status 
if(nrow(carrier.res.couple) == 0){
  carrier.res.couple <- data.frame(Message = "No potential recessive disease caused by clinical pathogenic or compound heterozygous variants to offspring based on the carrier screening.")
}
write.table(carrier.res.couple, paste0(output_dir, "/Results/carrier_screening_result_couple.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

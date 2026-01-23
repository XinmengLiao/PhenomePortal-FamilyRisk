## ---- Load libraries ----
library(kinship2)
library(tidyr)
library(dplyr)
librring(stringr)
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
sampleid <- args[1]
input_file <- args[2]
gender <- args[3]
output_dir <-  args[4]
genelist <- args[5]

print(paste0("Now pocessing ", sampleid))
print(paste0("Input file: ", input_file))
print(paste0("Gender: ", gender))
print(paste0("Output directory: ", output_dir))

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
result <- read.csv(input_file,header = T,sep = "\t") %>% unique()
colnames(result) <- gsub("am_","AlphaMissense_", colnames(result))
if ("NBScreening" %in% genelist){
  remove.variant <- read.csv(TR_removed_variant, header = T, sep = "\t")
  result <- result %>% filter(!(variant_info %in% remove.variant$variant_info)) %>%
    filter(Genes != "MEFV" & HGVSp != "p.Met694Val") %>% 
    filter(Genes != "MEFV" & HGVSp != "p.Met680Ile")
}

## ---- 2. Functions ----
# Function 1: Potential recessive disease on the offspring 
carrier.single <- function(df, gender){
  plp.carrier <- df %>% 
    mutate(Gender = gender) %>% 
    filter(grepl("Pathogenic|Likely_pathogenic",ClinVar_CLNSIG) & 
           grepl("Pathogenic|Likely_pathogenic",acmg_classification)) %>%  
    filter(Inheritance %in% c("AR", "XLR")) %>% unique()

  if (nrow(plp.carrier) == 0){
    return(plp.carrier)
  }else{
    genotype_col <- sapply(df, function(x) any(grepl("^[0-1]/[0-1]:", x)))
    genotype_col_name <- names(df)[genotype_col]
    
    # Zygosity function 
    zygosity.fun <- function(genotype, gender, chrom){
      if(gender == "Male" & chrom == "X" & genotype %in% c("0/1","1/0","1|0","0|1")){
        return("Hemizygous")
      } else if(genotype %in% c("0/1","1/0","1|0","0|1")){
        return("Heterozygous")
      } else if(genotype %in% c("1/1","1|1")){
        return("Homozygous")
      } else {
        return("Wrong Genotype")
      }
    }
    
    plp.carrier.final <- plp.carrier %>% 
      mutate(Genotype = substr(.data[[genotype_col_name]], 1, 3)) %>% 
      rowwise() %>%
      mutate(Zygosity = zygosity.fun(Genotype, .data[["Gender"]], X.CHROM)) %>%
      ungroup() %>%
      mutate(HGVSc = gsub("^.*:","",HGVSc), HGVSp = gsub("^.*:","",HGVSp)) %>% 
      select(Disease, Inheritance, Genes, X.CHROM, POS, REF, ALT, HGVSc, HGVSp,Zygosity, ClinVar_CLNSIG, acmg_classification, 
             IMPACT, MAX_AF) %>% unique() 
    
    return(plp.carrier.final)
  }
  
}

## ---- Generate screening results ----
carrier.res.single <-  carrier.single(df = result, gender = gender)

## ---- Writing output tables ----
# scren-positive monogenic disease and carrier status 
if(nrow(carrier.res.single) == 0){
  carrier.res.single <- data.frame(Message = "No potential recessive disease caused by clinical pathogenic or compound heterozygous variants to offspring based on the carrier screening.")
}
write.table(carrier.res.single, paste0(output_dir, "/Results/carrier_screening_result_single.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

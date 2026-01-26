## ---- Load libraries ----
library(kinship2)
library(tidyr)
library(dplyr)
library(ggplot2)
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
sampleid <- args[1]
input_file <- args[2]
ped_file <- args[3]
output_dir <-  args[4]
genelist <- args[5]

results_dir <- file.path(output_dir, "Results")
dir.create(results_dir, showWarnings = FALSE)


if (location == "local"){
  compare_file <- '/Users/xinmengliao/Documents/Project/20250710_FamilyRisk/Datasets/NBSeq_Results.xlsx'
  TR_removed_variant <- '/Users/xinmengliao/Documents/Project/20250710_FamilyRisk/Datasets/TRpipelineRemovedVariants.txt'
  genedb_file <- "/Users/xinmengliao/Documents/Project/20250710_FamilyRisk/Datasets/genelists/Preset_screening_list_GenCC_20251125.txt"
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
result <- read.csv(input_file,header = T,sep = "\t") %>% unique()
colnames(result) <- gsub("am_","AlphaMissense_", colnames(result))
if ("NBScreening" %in% genelist){
  remove.variant <- read.csv(TR_removed_variant, header = T, sep = "\t")
  result <- result %>% filter(!(variant_info %in% remove.variant$variant_info))
}

ped <- read.csv(ped_file, header = F,sep = "\t") %>% unique()
colnames(ped) <- c("FID","Individual","Father","Mother","Sex","Genotype")
familyid <- ped$FID %>% unique()

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

## ---- Family and Variant Summary ----
family_info <- ped %>% 
  mutate(Relationship = case_when((Father == 0 & Mother == 0 & Sex == 2) ~ "Mother",
                                  (Father == 0 & Mother == 0 & Sex == 1) ~ "Father", 
                                  (Father != 0 & Mother != 0) ~ paste0("Kid", cumsum(Father != 0 & Mother != 0)), T ~ NA),
         Sex = if_else(Sex == 2, "Female", "Male")) %>% select(-Genotype) %>% 
  select(Relationship, Individual, Sex) %>% rename(sampleID = Individual)


plp.count <- result %>% 
  filter(grepl("Pathogenic|Likely_pathogenic",ClinVar_CLNSIG) | grepl("Pathogenic|Likely_pathogenic",acmg_classification)) %>% 
  select(variant_info) %>% unique()

var_summary <- data.frame(Content = c("Filtered variants", "Clinical P/LP variants"),
                          Count = c(length(unique(result$variant_info)), 
                                    nrow(plp.count)))

## ---- 2. Functions ----

# Function 1: filter out the compound heterozygous for recessive disease and add back to the df.
# This results will be merged with the screen-positive results
compound.heterozygous.res <- function(df){
	index <- grep("Kid", colnames(df))
	kid_pattern_col <- intersect(grep("Pattern",colnames(df)), index)
	kid_genotype_col <- intersect(grep("Genotype",colnames(df)), index)
	kid_cols_to_keep <- c(kid_pattern_col, kid_genotype_col)
    
	brief.tmp <- df %>%
		filter(Inheritance %in% c("AR", "XLR")) %>%
    filter(rowSums(!is.na(select(., all_of(kid_genotype_col)))) > 0) %>% 
    mutate(ClinVar_CLNSIG = gsub("_", " ", ClinVar_CLNSIG),
        acmg_classification = gsub("_", " ", acmg_classification),
        ClinVar_CLNSIG = sapply(strsplit(split = "&", ClinVar_CLNSIG), `[`, 1),
        ClinVar_CLNSIG = if_else(ClinVar_CLNSIG == "Pathogenic/Likely pathogenic/Pathogenic",
                                 "Pathogenic/Likely pathogenic",
                                 as.character(ClinVar_CLNSIG)),
        Existing_variation = sapply(strsplit(split = "&", Existing_variation), `[`, 1)) %>% 
	  select(Disease, Genes, Inheritance, variant_info, Existing_variation, ClinVar_CLNSIG, acmg_classification, FatherGenotype, MotherGenotype,all_of(kid_cols_to_keep)) %>%
	  rename(`Existing variant` = Existing_variation, ClinVar = ClinVar_CLNSIG, ACMG = acmg_classification) %>% distinct()
	
	
	# function for infer the compound heterozygous origin
	infer_origin <- function(kid, father, mother) {
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
	  }else if (father %in% c("0/0","0|0") & mother %in% c("0/0","0|0")){
	    return("de novo")
	  }else if (father %in% c("./.",".|.")){
	    return("Mother")
	  }else if (mother %in% c("./.",".|.")) {
	    return("Father")
	  }else if (father == "1" & mother == "."){
	    return("Uncertained chrX phenotype")
	  }else{
	    return("Uncertained")
	  }
	}
	
	ch.df.all <-data.frame(Genes = character(), variant_info = character(), stringsAsFactors = FALSE)

	for ( i in kid_genotype_col){
	  kid_col <- colnames(df)[i]
	  kid_origin_col <- gsub("Genotype","Origin",kid_col)
	  kid_ch_col <- gsub("Genotype","Compound_heterozygous",kid_col)
	  CH_colname <- gsub("Genotype","Compound_heterozygous",kid_col)
	  
	  brief.com.het1 <- brief.tmp %>%
	    select(Genes, variant_info, ClinVar, ACMG,FatherGenotype, MotherGenotype, all_of(kid_col)) %>% 
	    rowwise() %>%
	    mutate(origin = infer_origin(.data[[kid_col]], FatherGenotype, MotherGenotype)) %>%
	    ungroup() %>%
	    filter(.data[[kid_col]] %in% c("0/1","1/0","1|0","0|1")) %>% unique() %>% 
	    filter(ClinVar %in% c("Pathogenic/Likely pathogenic","Pathogenic","Likely pathogenic","Uncertain significance") &
	             ACMG %in% c("Pathogenic","Likely pathogenic","Uncertain significance")) %>% unique() %>% 
	  filter(origin != "de novo" )
	  
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
	  
	  if(nrow(ch.df) == 0){
	    ch.df = data.frame(origin = NA, Compound_heterozygous = NA) %>% 
	      rename(!!kid_origin_col := origin, !!kid_ch_col := Compound_heterozygous)
	    brief.com.het2 <- cbind(brief.com.het1, ch.df)
	    ch.df.all <- full_join(ch.df.all, brief.com.het2)
	    kid_col1 = colnames(ch.df.all)[grep("Origin|Compound", colnames(ch.df.all))]
	    df <- df %>% left_join(., ch.df.all %>% select(Genes, variant_info, all_of(kid_col1)) %>% unique(), by = c("Genes","variant_info")) 
	  }else{
	    brief.com.het2 <- brief.com.het1 %>% 
	      left_join(., ch.df) %>% select(-origins) %>% 
	      filter(Compound_heterozygous %in% c("Yes","Uncertained","Uncertained chrX phenotype")) %>% 
	      rename(!!kid_origin_col := origin, !!kid_ch_col := Compound_heterozygous)
	    ch.df.all <- full_join(ch.df.all, brief.com.het2, by = c("Genes","variant_info"))
	    kid_col1 = colnames(ch.df.all)[grep("Origin|Compound", colnames(ch.df.all))]
	    df <- df %>% left_join(., ch.df.all %>% select(Genes, variant_info, all_of(kid_col1)) %>% unique(), by = c("Genes","variant_info")) 
	  }
	}
	
	index1 <- grep("Kid", colnames(df))
	kid_compoundhet_col <- intersect(grep("Compound_heterozygous",colnames(df)), index1)
	df <- df %>% 
	  filter(if_any(all_of(kid_compoundhet_col),
	                ~ .x %in% c("Yes", "Uncertained","Uncertained chrX phenotype")))
	
	return(df)
}

# Function 2: Generate the pedigree plot for each clinical P/LP variant
plot_pedig.res <- function(df, variant_info, ped){
  
  for ( i in variant_info){
    variant = i
    print(paste0("Now ploting pedigree plot for ", i))
    
    row <- df %>% filter(variant_info == variant)
    ped_data <- ped
    transcript <- if_else(is.na(unlist(strsplit(row$HGVSc,split = ":"))[2]), "", unlist(strsplit(row$HGVSc,split = ":"))[2])
    protein <- if_else(is.na(unlist(strsplit(row$HGVSp,split = ":"))[2]), "", paste0("(",unlist(strsplit(row$HGVSp,split = ":"))[2],")"))
    gene <- row$SYMBOL %>% unique()
    kid_gender_col <- colnames(result)[grep("Kid.*Gender",colnames(result))]
    kid_genotype_col <- colnames(df)[grep("Kid.*Genotype",colnames(df))]
    
    for (k in 1:length(kid_gender_col)) {
      kidid <- unlist(strsplit(kid_gender_col[k], split = "Gender_"))[2]
      genotype = row[1,kid_genotype_col[k]]
      ped_data$Genotype[which(ped_data$Individual == kidid)] = genotype
    }
    
    # for father and mother genotype 
    f_genotype = row[1,"FatherGenotype"]
    ped_data$Genotype[which(ped_data$Father ==0 & ped_data$Mother == 0 & ped_data$Sex == 1)] = f_genotype
    m_genotype = row[1,"MotherGenotype"]
    ped_data$Genotype[which(ped_data$Father ==0 & ped_data$Mother == 0 & ped_data$Sex == 2)] = m_genotype
    ped_data$Father[ped_data$Father == "0"] <- NA
    ped_data$Mother[ped_data$Mother == "0"] <- NA
    ped_data$Sex <- as.integer(ped_data$Sex)
    ped_data <- ped_data %>% tidyr::separate(Genotype, into = c("affected","avail"),sep = "[/|]")
    ped_data$Mother[ped_data$Mother == "0"] <- NA
    ped_data$affected <- as.numeric(ped_data$affected)
    ped_data$avail <- as.numeric(ped_data$avail)
    
    ped_res <- pedigree(id = ped_data$Individual, dadid = ped_data$Father, momid = ped_data$Mother, sex = ped_data$Sex,
                        affected = cbind(ped_data$affected, ped_data$avail), famid = ped_data$FID)
    ped1 <- ped_res[familyid]
    
    png(paste0(output_dir,"/Results/", familyid, "_", gene,"_",transcript, "_pedigree_plot.png"), 
        width=800, height=900,res = 300)
    plot(ped1, col=ifelse(ped_data$avail, 2, 1), cex=0.6)
    #title(main="Pedigree analysis",cex.main = 0.7)
    mtext(paste(gene, transcript, protein), side=3, line=0.5, cex=0.5)
    dev.off()
  }
}

# Function 3: Manage the screen-positive diseases with compound heterozygous
positive.monogenic.plp.family <- function(df, ped){
  # analyse compound heterozygous results
  result.ch <- compound.heterozygous.res(df)
  
  # analyse monogenic positive results
  kid_pattern_cols <- colnames(df)[grep("Kid.*Pattern", colnames(df))]
  kid_genotype_cols <- colnames(df)[grep("Kid.*Genotype", colnames(df))]
  plp.monogenic <- df %>% 
    filter(
      ((grepl("Pathogenic|Likely_pathogenic",ClinVar_CLNSIG) & grepl("Pathogenic|Likely_pathogenic",acmg_classification)) & 
         if_any(all_of(kid_pattern_cols), ~ grepl("dominant", .x, ignore.case = TRUE)))) %>% unique()
  
  if(nrow(plp.monogenic) > 0){
    ## generate the pedigree for clinical P/LP variants by default, not for the compound heterozygous ones
    positive.res.pedigree <- plp.monogenic %>% 
      filter(grepl("Pathogenic|Likely_pathogenic",ClinVar_CLNSIG) & grepl("Pathogenic|Likely_pathogenic",acmg_classification))
    plot_pedig.res(df = positive.res.pedigree, ped = ped, variant_info = positive.res.pedigree$variant_info)
  }
  
  # merge two results together
  plp.monogenic.final <- rbind(plp.monogenic, result.ch) %>% unique()
  
  if(nrow(plp.monogenic.final) > 0 ){
    plp.monogenic.final <- plp.monogenic.final %>% 
      mutate(HGVSc = gsub("^.*:","",HGVSc), HGVSp = gsub("^.*:","",HGVSp)) %>% 
      select(Disease, Inheritance, Genes, X.CHROM, POS, REF, ALT, HGVSc, HGVSp, ClinVar_CLNSIG, acmg_classification, 
             IMPACT, MAX_AF, FatherGenotype, MotherGenotype, all_of(kid_genotype_cols), all_of(kid_ch_cols)) %>% unique() %>% 
      left_join(., compare_positive %>% 
                  select(Reported.Disease, `Results from Other NBSeq Projects`, Genes, Transcript) %>% unique(),
                by = c("Genes", "HGVSc" = "Transcript")) %>% 
      rename(`Reported Disease from Other NBSeq Projects` = Reported.Disease) %>% unique()
  }
  
  return(plp.monogenic.final)
}

# Function 4: Manage the carrier status for autosomal recessive and X-linked recessive diseases
carrier.plp.family <- function(data){
  
  kid_pattern_cols <- colnames(data)[grep("Kid.*Pattern", colnames(data))]
  kid_genotype_cols <- colnames(data)[grep("Kid.*Genotype", colnames(data))]
  
  carrier <- data %>% 
    filter(grepl("Pathogenic|Likely_pathogneic",ClinVar_CLNSIG) &
           grepl("Pathogenic|Likely_pathogneic",acmg_classification)) %>% 
    filter(Inheritance %in% c("AR","XLR")) %>% 
    filter(if_any(all_of(kid_pattern_cols), ~ grepl("recessive", .x, ignore.case = TRUE)))
  
  carrier.disease <- carrier %>% 
    select(Disease, Inheritance, Genes, X.CHROM, POS, REF, ALT, HGVSc, HGVSp, ClinVar_CLNSIG, acmg_classification, 
           IMPACT, MAX_AF, FatherGenotype, MotherGenotype, all_of(kid_genotype_cols)) %>% unique() %>% 
    left_join(., compare_carrier %>% 
                select(Reported.Disease, `Results from Other NBSeq Projects`, Genes, Transcript) %>% unique(),
              by = c("Genes", "HGVSc" = "Transcript")) %>% 
    rename(`Reported Disease from Other NBSeq Projects` = Reported.Disease) %>% unique()
  
  return(carrier.disease)
}


## ---- Generate screening results ----
positive.res <-  positive.monogenic.plp.family(df = result, ped = ped)

# carrier status 
if("NBScreening" %in% genelist){
  result_carrier <- result %>% 
    filter(Genes != "MEFV" & HGVSp != "p.Met694Val") %>% 
    filter(Genes != "MEFV" & HGVSp != "p.Met680Ile")
}
carrier.res <- carrier.plp.family(result_carrier)

## ---- Writing output tables ----
# Family and variant summary 
write.table(family_info, paste0(output_dir, "/Results/family_information.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(var_summary, paste0(output_dir, "/Results/variant_summary.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

# scren-positive monogenic disease and carrier status 
if(nrow(positive.res) == 0){
  positive.res <- data.frame(Message = "No positive monogenic disease variant detected based on the current gene panel and filtering criteria.")
}
if(nrow(carrier.res) == 0){
  carrier.res <- data.frame(Message = "No positive monogenic disease carrier status variant detected based on the current gene panel and filtering criteria.")
}
write.table(positive.res, paste0(output_dir, "/Results/monogenic_positive_results.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(carrier.res, paste0(output_dir, "/Results/carrier_status_results.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

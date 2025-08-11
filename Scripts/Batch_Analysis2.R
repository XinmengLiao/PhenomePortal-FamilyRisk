library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(viridis)
library(biomaRt)
library(ggthemes)
library(patchwork)


for (i in c("select","filter", "mutate","rename", "left_join", "slice")){
  conflicted::conflict_prefer(i, "dplyr")
}
rm(i)
conflicted::conflicts_prefer(stats::sd)
conflicted::conflicts_prefer(httr::content)
conflicted::conflicts_prefer(plotly::layout)

sampleID = "newborn103"
result_file = "/Users/xinmengliao/Documents/Project/20250710_NewbornRisk/Batch/newborn103_vep_merged_rmmissingalt_biallelic2.txt"
metadata_file = "/Users/xinmengliao/Documents/Project/20250710_NewbornRisk/Batch/newborn103gender.txt"
output_dir = "/Users/xinmengliao/Desktop/VarXOmics server scripts/NewbornRisk"
genelist = "TR"
pgx_list <- list.files(paste0(output_dir, "PGX"), full.names = TRUE, pattern = "PGx")

args <- commandArgs(trailingOnly = TRUE)
sampleID <- args[1]
result_file <- args[2]
metadata_file <- args[3]
output_dir <- args[4]
output_dir <- paste0(output_dir, "/")
genelist <- args[5] # TR, babyseq, earlycheck, babydetect, guardian, babyscreen, ACMG, customized gene list
pgx_list <- list.files(paste0(output_dir, "PGX"), full.names = TRUE, pattern = "PGx")


## ---- Load files  ----
result <- read.csv(result_file,header = T,sep = "\t") 
metadata <- read.csv(metadata_file,header = F,sep = "\t") # include sample ID and gender 
colnames(metadata) <- c("Sample","Gender")
total_samplenum <- nrow(metadata)

# GeneDB list 
path = "sysmed" # local or sysmed
if (path == "local"){
  genelist_paths <- list(
    'babyseq' = '/Users/xinmengliao/Documents/Project/20250710_NewbornRisk/Datasets/genelists/BabySeq.txt',
    'earlycheck' = '/Users/xinmengliao/Documents/Project/20250710_NewbornRisk/Datasets/genelists/EarlyCheck.txt',
    'babydetect' = '/Users/xinmengliao/Documents/Project/20250710_NewbornRisk/Datasets/genelists/BabyDetect.txt',
    'babyscreen' = '/Users/xinmengliao/Documents/Project/20250710_NewbornRisk/Datasets/genelists/BabyScreen.txt',
    'gemonic101' = '/Users/xinmengliao/Documents/Project/20250710_NewbornRisk/Datasets/genelists/Genomic101.txt',
    'TR' = '/Users/xinmengliao/Documents/Project/20250710_NewbornRisk/Datasets/genelists/NBScreening.txt',
    'ACMG' = '/Users/xinmengliao/Documents/Project/20250710_NewbornRisk/Datasets/genelists/ACMGv3.3.txt',
    "TR_removed_variant" = "/Users/xinmengliao/Documents/Project/20250710_NewbornRisk/Scripts/batch/TRpipelineRemovedVariants.txt",
    "compare_file" = '/Users/xinmengliao/Documents/Project/20250710_NewbornRisk/Datasets/NBSeq_Results.xlsx'
  )
}else if(path == "sysmed"){
  genelist_paths <- list(
    'babyseq' = '/mnt/storage_pool/Genomics/Genome/NewbornRisk/Datasets/genelists/BabySeq.txt',
    'earlycheck' = '/mnt/storage_pool/Genomics/Genome/NewbornRisk/Datasets/genelists/EarlyCheck.txt',
    'babydetect' = '/mnt/storage_pool/Genomics/Genome/NewbornRisk/Datasets/genelists/BabyDetect.txt',
    'babyscreen' = '/mnt/storage_pool/Genomics/Genome/NewbornRisk/Datasets/genelists/BabyScreen.txt',
    'gemonic101' = '/mnt/storage_pool/Genomics/Genome/NewbornRisk/Datasets/genelists/Genomic101.txt',
    "guardian" = '/mnt/storage_pool/Genomics/Genome/NewbornRisk/Datasets/genelists/Guarduan.txt',
    'TR' = '/mnt/storage_pool/Genomics/Genome/NewbornRisk/Datasets/genelists/NBScreening.txt',
    'ACMG' = '/mnt/storage_pool/Genomics/Genome/NewbornRisk/Datasets/genelists/ACMGv3.3.txt',
    "TR_removed_variant" = "/mnt/storage_pool/Genomics/Genome/NewbornRisk/Datasets/TRpipelineRemovedVariants.txt",
    "compare_file" = '/mnt/storage_pool/Genomics/Genome/NewbornRisk/Datasets/NBSeq_Results.xlsx'
  )
}

if (genelist == "babyseq") {
  genedb <- read.delim(genelist_paths[["babyseq"]])
} else if (genelist == "earlycheck") {
  genedb <- read.delim(genelist_paths[["earlycheck"]])
} else if (genelist == "babydetect") {
  genedb <- read.delim(genelist_paths[["babydetect"]])
} else if (genelist == "babyscreen") {
  genedb <- read.delim(genelist_paths[["babyscreen"]])
} else if (genelist == "guardian") {
  genedb <- read.delim(genelist_paths[["guardian"]])
} else if (genelist == "TR") {
  genedb <- read.delim(genelist_paths[["TR"]])
} else if (genelist == "ACMG") {
  genedb <- read.delim(genelist_paths[["ACMG"]])
}

# remove the variants deleted by root cause analysis 
if (genelist == "TR"){
  remove.variant <- read.csv(genelist_paths[["TR_removed_variant"]],header = T,sep = "\t")
  result <- result %>% filter(!(variant_info %in% remove.variant$variant_info))
}


## ---- Functions ----
genotype.count <- function(data, sampleid, metadata){
  sampleid <- as.character(sampleid)
  data$X.CHROM <- as.character(data$X.CHROM)
  data$Inheritance <- as.character(data$Inheritance)
  
  # Clean the genotype
  data[, sampleid] <- lapply(data[, sampleid], function(col) {
    col <- as.character(col)
    col[is.na(col)] <- ""
    sapply(strsplit(col, ":"), `[`, 1)
  })
  
  # Convert back to data frame if needed
  data[, sampleid] <- as.data.frame(data[, sampleid])
  
  # Rest of your existing code unchanged...
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
  if (!identical(metadata$Gender, FALSE)){
    data <- data %>% left_join(., metadata, by = c("homo_sample" = "Sample")) %>% 
      rename(homo_gender = Gender) %>% 
      left_join(., metadata, by = c("het_sample" = "Sample")) %>% 
      rename(het_gender = Gender) %>% unique() 
  }
  
  data <- data %>% 
    mutate(het_zygosity = if_else(!is.na(het_gender), case_when((het_gender == "Male" & X.CHROM == "chrX" & 
                                                                   Inheritance %in% c("AR", "XLR")) ~ "Hemizygous",T ~ "Heterozygous"), NA_character_),
           homo_zygosity = if_else(!is.na(homo_gender), "Homozygous", NA_character_) ) %>% unique()
  
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
  
  plp <- genotype.count(plp, sampleid = metadata$Sample, metadata = metadata)
  
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
    group_by(variant_info) %>% 
    summarise(`No. of sample for variant` = n(), .groups = "drop") %>% 
    mutate(`Variant carrier frequency` = `No. of sample for variant` / total_samplenum )
  
  plp.monogenic.gene <- plp.monogenic %>%
    pivot_longer(cols = c(het_sample, homo_sample), names_to = "zygosity_type", values_to = "sample_id") %>% 
    filter(sample_id != "" & !is.na(sample_id)) %>% 
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
    mutate(HGVSc = gsub("^.*:","",HGVSc),
           HGVSp = gsub("^.*:","",HGVSp),
           ClinVar_CLNSIG = sapply(strsplit(split = "&", as.character(ClinVar_CLNSIG)),`[`,1),
    ) %>% 
    rename(Inh = Inheritance, Chr = X.CHROM, Pos = POS, Ref = REF, Alt = ALT, `Global Max AF` = MAX_AF,  
           `No. of sample for gene` = gene_count, `Gene carreir frequency` = proportion) %>% unique() %>% 
    group_by(Genes) %>% 
    mutate(Disease = paste(unique(Disease), collapse = "; ")) %>% ungroup() %>% unique() 
  
  # merge other frequencies together 
  plp.monogenic.gene.final <- plp.monogenic.gene1 %>% 
    left_join(., plp.disease, by = "Disease") %>% 
    mutate(variant_info = paste(Chr, Pos, Ref, Alt, sep = "_")) %>% 
    left_join(., variant.freq, by = "variant_info") %>% select(-variant_info) %>% 
    select(Disease, Inh, Genes, Chr, Pos, Ref, Alt, Zygosity, HGVSc, HGVSp, ClinVar_CLNSIG, acmg_classification, 
           IMPACT, `Global Max AF`, `No. of sample for variant`, `Variant carrier frequency`,
           `No. of sample for gene`, `Gene carreir frequency`, 
           `No. of sample for disease`, `Disease carreir frequency`, )
  
  return(list(monogenic.positive = plp.monogenic.gene.final, monogenic.variant = plp.monogenic$variant_info, 
              raw_data = plp.monogenic))
}

carrier.plp <- function(data){
  carrier <- data %>% filter(Inheritance %in% c("AR","XLR")) 
  carrier <- genotype.count(carrier, sampleid = metadata$Sample, metadata = metadata) %>% 
    filter(!(X.CHROM == "chrX" & het_zygosity == "Hemizygous"))
  
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
    mutate(HGVSc = gsub("^.*:","",HGVSc),
           HGVSp = gsub("^.*:","",HGVSp),
           ClinVar_CLNSIG = sapply(strsplit(split = "&", as.character(ClinVar_CLNSIG)),`[`,1),
    ) %>% 
    rename(Inh = Inheritance, Chr = X.CHROM, Pos = POS, Ref = REF, Alt = ALT, `Global Max AF` = MAX_AF,  
           `No. of sample for gene` = gene_count, `Gene carreir frequency` = proportion) %>% unique() %>% 
    group_by(Genes) %>% 
    mutate(Disease = paste(unique(Disease), collapse = "; ")) %>% ungroup() %>% unique() 
  
  # merge other frequencies together 
  carrier.status.gene.final <- carrier.status.gene1 %>% 
    left_join(., carrier.disease, by = "Disease") %>% 
    mutate(variant_info = paste(Chr, Pos, Ref, Alt, sep = "_")) %>% 
    left_join(., variant.freq, by = "variant_info") %>% select(-variant_info) %>% 
    select(Disease, Inh, Genes, Chr, Pos, Ref, Alt, Zygosity, HGVSc, HGVSp, ClinVar_CLNSIG, acmg_classification, 
           IMPACT, `Global Max AF`, `No. of sample for variant`, `Variant carrier frequency`,
           `No. of sample for gene`, `Gene carreir frequency`, 
           `No. of sample for disease`, `Disease carreir frequency`, )
  
  return(list(carrier_final = carrier.status.gene.final, carrier_raw = carrier.status))
}


## ---- Monogenic positive cases and carrier status cases (Tables) ----
#### ---- monogenic positive ----
plp.positive <- positive.monogenic.plp(result)$monogenic.positive
plp.variant <- positive.monogenic.plp(result)$monogenic.variant

#### ---- carrier status  ----
if(genelist == "TR"){
  result_carrier <- result %>% filter(!(Genes == "MEFV" & HGVSp == "NP_000234.1:p.Met694Val")) %>% 
    filter(!(Genes == "MEFV" & HGVSp == "NP_000234.1:p.Met680Ile")) %>% 
    filter(!(variant_info %in% plp.variant))
}
carrier.status <- carrier.plp(result_carrier)$carrier_final

write.table(plp.positive, paste0(output_dir, "/monogenic_positive_results.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(carrier.status, paste0(output_dir, "/carrier_status_results.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

monogenic.result <- positive.monogenic.plp(result)$raw_data
carrier.result <- carrier.plp(result_carrier)$carrier_raw


## ---- bar plots for monogenic disease frequency and variant number ----
# count gene length
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_lengths <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", 
                                     "chromosome_name", "start_position", 
                                     "end_position"),
                      filters = "biotype",
                      values = "protein_coding",
                      mart = ensembl)
gene_lengths$gene_length <- gene_lengths$end_position - gene_lengths$start_position + 1

number <- plp.positive %>% arrange(desc(`No. of sample for gene`)) %>% 
  select(Genes, `No. of sample for gene`) %>% unique()
plp_number <- plp.positive %>% 
  arrange(desc(`No. of sample for gene`)) %>% 
  rename(Pathogenicity = ClinVar_CLNSIG) %>% 
  mutate(Pathogenicity = if_else(Pathogenicity == "Pathogenic/Likely_pathogenic" , "Pathogenic", Pathogenicity)) %>% 
  unique() %>% 
  group_by(Genes, Pathogenicity) %>% 
  summarise(count = n(), .groups = 'drop') %>% 
  filter(Pathogenicity %in% c("Pathogenic", "Likely_pathogenic")) %>%
  filter(Genes %in% number$Genes) %>% 
  left_join(gene_lengths[, c("hgnc_symbol", "gene_length")], by = c("Genes" = "hgnc_symbol")) %>% unique() 

if("Pathogenic" %in% unique(plp_number$Pathogenicity) & 
   !("Likely_pathogenic" %in% unique(plp_number$Pathogenicity))){
  plp_number <- rbind(plp_number, plp_number %>% mutate(Pathogenicity = "Likely_pathogenic", count = 0))
}else if ("Likely_pathogenic" %in% unique(plp_number$Pathogenicity) & 
          !("Pathogenic" %in% unique(plp_number$Pathogenicity))){
  plp_number <- rbind(plp_number, plp_number %>% mutate(Pathogenicity = "Pathogenic", count = 0))
}

plp_number <- plp_number %>% 
  pivot_wider(names_from = Pathogenicity, values_from = count, values_fill = list(count = 0)) %>% 
  pivot_longer(cols = c(Pathogenic, Likely_pathogenic), names_to = "Pathogenicity", values_to = "count") %>%
  mutate(P.adjust = (count / gene_length) * 100000)

plp_number$Genes <- factor(plp_number$Genes, 
                           levels = rev(number$Genes[order(number$`No. of sample for gene`, decreasing = TRUE)]))

number.disease <- plp.positive %>% select(Genes, Disease) %>% unique() %>% 
  group_by(Genes) %>% summarise(Disease = paste(unique(Disease), collapse = ", "), .groups = 'drop') %>% 
  mutate(Disease = paste(Genes, Disease, sep = " "))
number <- number %>% left_join(., number.disease) %>% 
  select(-Genes) %>% rename(Genes = Disease)
number$Genes <- factor(number$Genes, 
                       levels = rev(number$Genes[order(number$`No. of sample for gene`, decreasing = TRUE)]))

p1 <- ggplot(plp_number, aes(y = Genes, x = P.adjust, fill = Genes,alpha = Pathogenicity)) +
  geom_bar(stat = "identity") +
  labs( x = "Adjusted Variant Number", y = NULL, fill = NULL)+
  theme_classic() +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
  )+
  scale_x_reverse()+
  scale_y_discrete(expand = c(0, 0), position = "right") +  
  #scale_fill_viridis(option = "C", discrete = TRUE)+
  #scale_fill_tableau(palette = "Tableau 20",direction = -1) + 
  scale_alpha_manual(values = c("P" = 1, "LP" = 0.6))

if (length(unique(plp_number$Genes)) <= 20) {
  p1 <- p1 + scale_fill_tableau(palette = "Tableau 20", direction = -1)
} else {
  p1 <- p1 + scale_fill_viridis(option = "C", discrete = TRUE)
}

p2 <- ggplot(number, aes(y = Genes, x = "")) +
  geom_tile(aes(fill =`No. of sample for gene`)) +
  theme_classic()+
  scale_fill_gradient(low="#dadaeb", high="darkblue") +
  labs(
    x = NULL,
    y = NULL, fill = str_wrap("Mutated Gene Frequency", width = 10)
  )+
  scale_x_reverse()+
  scale_y_discrete(position = "right") +
  scale_x_discrete(position = "bottom")+
  labs(x = str_wrap("Gene mutated frequency", width = 10))+
  theme(#axis.title.x = element_text(size = 10,family = "Helvetica",color = "black"), 
    axis.title.y = element_blank(), axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(family = "Helvetica",color = "black", angle = 0 ,hjust =1,vjust = -3),
    axis.text.y = element_text(family = "Helvetica",color = "black"), 
    legend.title = element_text(size = 10, family = "Helvetica",color = "black"),
    legend.text = element_text(size = 10, family = "Helvetica",color = "black"))+
  theme(legend.position = "top")

wid <- length(unique(number$Genes))

png(paste0(output_dir,"Figures/Monogenic_details.png"),width = wid*4,height = wid*2,res = 600,units = "cm")

p1/ p2 + plot_layout(ncol = 2, widths = c( 6, 0.5))

dev.off()

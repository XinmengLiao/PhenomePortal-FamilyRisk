library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(fmsb)
library(viridis)
library(fmsb)
library(VennDiagram)
library(xlsx)
library(png)
library(grid)

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
pgx_list <- list.files(paste0(output_dir, "PGx"), full.names = TRUE, pattern = "PGx")


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

# load comparison files 
compare_positive <- read.xlsx(genelist_paths[['compare_file']],sheetIndex = 1) %>% filter(!is.na(Reported.Disease)) %>% select(1:8)
compare_guardian_positive <- read.xlsx(genelist_paths[['compare_file']],sheetIndex = 2) %>% filter(!is.na(Reported.Disease))
compare_babyseq_carrier <- read.xlsx(genelist_paths[['compare_file']],sheetIndex = 3,startRow = 2) %>% filter(!is.na(Reported.Disease))
compare_earlycheck_carrier <- read.xlsx(genelist_paths[['compare_file']],sheetIndex = 4) %>% filter(!is.na(Reported.Disease))


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
  carrier <- data %>% filter(Inheritance %in% c("AR","XLR")) %>% 
    filter(grepl("Pathogenic|Likely_pathogneic",ClinVar_CLNSIG) &
             grepl("Pathogenic|Likely_pathogneic",acmg_classification))
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


## ---- Compare results with other NBSeq projects (Table) ----
#### ---- positive monogenic ----
compare_positive_result1 <- plp.positive %>%
  left_join(., compare_positive, by = c("Genes", "HGVSc" = "Transcript")) %>% distinct() %>% 
  filter(!is.na(Project))

if(nrow(compare_positive_result1) > 0 ){
  write.table(compare_positive_result1, 
              paste0(output_dir, "/compare_positive_results.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
}

compare_positive_result2 <- plp.positive %>% 
  left_join(., compare_guardian_positive, by = c("Genes", "HGVSc" = "Transcript")) %>% distinct() %>% 
  filter(!is.na(Project))

if(nrow(compare_positive_result2) > 0 ){
  write.table(compare_positive_result2, 
              paste0(output_dir, "/compare_guardian_positive_results.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
}


#### ---- carrier status ----
compare_carrier_result1 <- carrier.status %>%
  left_join(., compare_babyseq_carrier, by = c("Genes", "HGVSc" = "Transcript")) %>% distinct() %>% 
  filter(!is.na(Project))

if(nrow(compare_carrier_result1) > 0 ){
  write.table(compare_carrier_result1, 
              paste0(output_dir, "/compare_babyseq_carrier_results.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
}

compare_carrier_result2 <- carrier.status %>% 
  left_join(., compare_earlycheck_carrier, by = c("Genes", "HGVSc" = "Transcript")) %>% distinct() %>% 
  filter(!is.na(Project))

if(nrow(compare_carrier_result2) > 0 ){
  write.table(compare_carrier_result2, 
              paste0(output_dir, "/compare_earlycheck_carrier_results.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
}


#### ---- Compare results with other NBSeq projects (Venn for overlapped genes) ----
compare_positive_list <- list(
  BabyDetect = compare_positive %>% filter(Project == "BabyDetect") %>% select(Genes) %>% pull(),
  BabySeq = compare_positive %>% filter(Project == "BabySeq") %>% select(Genes) %>% pull(),
  EarlyCheck = compare_positive %>% filter(Project == "EarlyCheck") %>% select(Genes) %>% pull(),
  Guardian = compare_guardian_positive %>% select(Genes) %>% pull(),
  Cohort = plp.positive %>% select(Genes) %>% pull()
)

compare_carrier_list <- list(
  BabySeq = compare_babyseq_carrier %>% select(Genes) %>% pull(),
  EarlyCheck = compare_earlycheck_carrier %>% select(Genes) %>% pull(),
  Cohort = carrier.status %>% select(Genes) %>% pull()
)


#### ---- Venn Plots ----
# for positive monogenic
venn.diagram(
  x = compare_positive_list,
  
  category.names = rep("", 5), 
  filename =paste0(output_dir, "Figures/Positive_comparison_venn.png"),
  
  main = "Monogenenic positives shared genes",
  main.cex = 1,
  main.fontface = "bold", 
  main.fontfamily = "Helvetica",
  
  scaled = T,
  
  # Output features  
  imagetype="png",
  units = "cm",
  height = 12,
  width = 16,  # 增加宽度为图例留空间
  resolution = 600,
  compression = "lzw",
  margin = 0.15,  # 增加边距
  
  # Circles
  lwd = 2,
  lty = 'blank', 
  fill = c("#f38181", "#c1b2f2","#fce38a", "#95e1d3","#8abdf0"),
  
  # Numbers
  cex = 1.6,
  fontfamily = "Helvetica",
  
  # 隐藏类别标签
  cat.cex = 0
)

#library(png)
library(grid)

# 读取生成的图片
img_path <- paste0(output_dir, "Figures/Positive_comparison_venn.png")
img <- readPNG(img_path)

png(paste0(output_dir, "Figures/Positive_comparison_venn_with_legend.png"),
    units = "cm", height = 12, width = 16, res = 600)

# set the margins
par(mar = c(0, 0, 0, 0))
plot(0:1, 0:1, type = "n", axes = FALSE, xlab = "", ylab = "")

# set the image margins
rasterImage(img, 0, 0, 0.75, 1)

par(new = TRUE, fig = c(0.75, 1, 0, 1), mar = c(0, 0, 0, 0))
plot(0:1, 0:1, type = "n", axes = FALSE, xlab = "", ylab = "")

# 图例内容
colors <- c("#f38181", "#c1b2f2","#fce38a", "#95e1d3","#8abdf0")
labels <- c("BabyDetect", "BabySeq", "EarlyCheck", "Guardian", "Cohort")

# 手动绘制图例
for(i in 1:5) {
  y_pos <- 0.6 - (i-1) * 0.06  # 缩短间距从0.12到0.08，并调整起始位置
  
  # 绘制小圆形（无边框）
  symbols(0.15, y_pos, circles = 0.05, 
          fg = colors[i], bg = colors[i], 
          inches = FALSE, add = TRUE)
  
  # 添加标签
  text(0.25, y_pos, labels[i], adj = 0, cex = 0.7)
}

dev.off()


# for carrier status 
venn.diagram(
  x = compare_carrier_list,
  
  category.names = rep("", 3), 
  filename =paste0(output_dir, "Figures/Carrier_comparison_venn.png"),
  
  main = "Carrier status shared genes",
  main.cex = 1,
  main.fontface = "bold", 
  main.fontfamily = "Helvetica",
  
  scaled = T,
  
  # Output features  
  imagetype="png",
  units = "cm",
  height = 12,
  width = 16,  # 增加宽度为图例留空间
  resolution = 600,
  compression = "lzw",
  margin = 0.15,  # 增加边距
  
  # Circles
  lwd = 2,
  lty = 'blank', 
  fill = c("#f38181","#c1b2f2", "#fce38a"),
  
  # Numbers
  cex = 1.6,
  fontfamily = "Helvetica",
  
  # 隐藏类别标签
  cat.cex = 0
)

# 读取生成的图片
img_path <- paste0(output_dir, "Figures/Carrier_comparison_venn.png")
img <- readPNG(img_path)

png(paste0(output_dir, "Figures/Carrier_comparison_venn_with_legend.png"),
    units = "cm", height = 12, width = 16, res = 600)

# set the margins
par(mar = c(0, 0, 0, 0))
plot(0:1, 0:1, type = "n", axes = FALSE, xlab = "", ylab = "")

# set the image margins
rasterImage(img, 0, 0, 0.75, 1)

par(new = TRUE, fig = c(0.75, 1, 0, 1), mar = c(0, 0, 0, 0))
plot(0:1, 0:1, type = "n", axes = FALSE, xlab = "", ylab = "")

# 图例内容
colors <- c("#f38181","#c1b2f2", "#fce38a")
labels <- c("BabyDetect", "BabySeq", "Cohort")

# 手动绘制图例
for(i in 1:5) {
  y_pos <- 0.6 - (i-1) * 0.06  # 缩短间距从0.12到0.08，并调整起始位置
  
  # 绘制小圆形（无边框）
  symbols(0.15, y_pos, circles = 0.05, 
          fg = colors[i], bg = colors[i], 
          inches = FALSE, add = TRUE)
  
  # 添加标签
  text(0.25, y_pos, labels[i], adj = 0, cex = 0.7)
}

dev.off()


## ---- Bar plots ----
# bar plots for carrier and monogenic disease count per sample
df.carrier <- carrier.plp(result_carrier)$carrier_raw %>% 
  select(variant_info, het_sample) %>% unique() %>% 
  group_by(het_sample) %>%
  summarise(count = n()) %>% ungroup() %>% 
  group_by(count) %>% 
  summarise(count1 = n()) %>% ungroup() %>% 
  rename(Var_num = count, Sample_num = count1) %>% mutate(Type = "Carrier status")
print(df.carrier)

df.monogenic =  positive.monogenic.plp(result)$raw_data %>% 
  pivot_longer(cols = c(het_sample, homo_sample), values_to = "sample_id") %>% 
  filter(sample_id != "" & !is.na(sample_id)) %>% 
  select(sample_id, variant_info) %>% unique() %>% 
  group_by(sample_id) %>% 
  summarise(Var_num = n_distinct(variant_info)) %>% ungroup() %>% 
  group_by(Var_num) %>% 
  summarise(Sample_num = n()) %>% ungroup() %>% mutate(Type = "Monogenic Diseases")
print(df.monogenic)

df_all <- rbind(df.carrier,df.monogenic)
df_all$Type = factor(df_all$Type, levels = unique(df_all$Type))
print(df_all)

get_plot_dimensions <- function(plot_obj) {
  plot_data <- plot_obj$data
  
  x_unique <- length(unique(plot_data$Var_num))
  y_max <- max(plot_data$Sample_num)
  
  width <- 3 + x_unique * 0.2
  height <- 3 + sqrt(y_max/100) * 1.5
  
  return(list(width = width, height = height))
}

# combined figure
max_var_num <- max(df_all$Var_num)
max_sample_num <- max(df_all$Sample_num)
p.all <- ggplot(df_all, aes(x = Var_num, y = Sample_num, fill = Type)) +
  geom_bar(stat = 'identity',position = 'dodge') +
  theme_classic() +
  scale_fill_manual(values = c("Carrier status" = "#B3DDF2", "Monogenic Diseases" = "#F2B3C6"))+
  labs(x = "No. of variant", y = "No. of Newborn") +
  scale_x_continuous(breaks = seq(0, max_var_num, 1)) +
  scale_y_continuous(breaks = seq(0, max_sample_num, round(max_sample_num/5))) +
  theme(
    axis.title.x = element_text(size = 20, family = "Helvetica", color = "black"),
    axis.title.y = element_text(size = 20, family = "Helvetica", color = "black"),
    axis.text.x = element_text(size = 18, family = "Helvetica", color = "black", margin = margin(b = 8)),
    axis.text.y = element_text(size = 18, family = "Helvetica", color = "black", margin = margin(l = 8)),
    strip.background = element_rect(fill = "#f0f0f0", color = "black"),
    strip.text = element_text(size = 15, color = "black", family = "Helvetica"),
    panel.grid.major.x = element_blank(),legend.title = element_blank(),
    legend.text = element_text(size = 10, family = "Helvetica", color = "black"),
    legend.position = "top"
  )

all_auto_dims <- get_plot_dimensions(p.all)

write.table(df_all, paste0(output_dir, "/monogenic_carrier_count.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
ggsave(p.all, filename = paste0(output_dir,"Figures/monogenic_carrier_count.pdf"), 
       width = all_auto_dims$width, height = all_auto_dims$height, dpi = 600,device = pdf)


## ---- Ontology (Radar plot) ----
# Pie chart 
#### ---- For Childhood onset diseases ----
on.child <- plp.positive %>% 
  mutate(variant_info = paste(Chr, Pos, Ref, Alt, sep = "_")) %>% 
  separate_rows(sep = "; ",Disease) %>% 
  left_join(., result %>% select(variant_info, Disorder_Group, Genes, Disease), by = c("Genes","Disease","variant_info")) %>% 
  separate_rows(sep = ";", Disorder_Group) %>% 
  select(Disease, Disorder_Group) %>% unique() %>% 
  group_by(Disorder_Group) %>% 
  summarise(count = n()) %>%
  arrange(desc(count)) %>% 
  mutate(ratio = count / sum(count)) %>% 
  mutate(
    percent = paste0(round(ratio * 100, 1), "%"), 
    cum_ratio = cumsum(ratio) - ratio / 2
  ) %>% as.data.frame() %>% mutate(ratio = as.numeric(ratio)) %>% na.omit(Disorder_Group)

write.table(on.child %>% select(Disorder_Group, count),
            file = paste0(output_dir,"/Positive.ontology.count.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

#### ---- For Carrier status ----
on.carrier<- carrier.status %>% 
  mutate(variant_info = paste(Chr, Pos, Ref, Alt, sep = "_")) %>% 
  separate_rows(sep = "; ",Disease) %>% 
  left_join(., result %>% select(variant_info, Disorder_Group, Genes, Disease), by = c("Genes","Disease","variant_info")) %>%
  separate_rows(sep = ";", Disorder_Group) %>% 
  select(Disease, Disorder_Group) %>% unique() %>% 
  mutate(Disorder_Group = if_else(grepl(";",Disorder_Group),"Complex genetic syndromes",Disorder_Group)) %>% unique() %>% 
  group_by(Disorder_Group) %>% 
  summarise(count = n()) %>%
  arrange(desc(count)) %>% 
  mutate(ratio = count / sum(count)) %>% 
  mutate(
    percent = paste0(round(ratio * 100, 1), "%"),  # 计算百分比标签
    cum_ratio = cumsum(ratio) - ratio / 2          # 扇形中间位置
  ) %>% as.data.frame() %>% mutate(ratio = as.numeric(ratio)) %>% na.omit(Disorder_Group) 

write.table(on.carrier %>% select(Disorder_Group, count),
            file = paste0(output_dir,"/Carrrier.ontology.count.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

# extract the number of ontology (Ontology count)
# 为什么这里有空的Disorder Group？
wide_data <- rbind(on.child %>% mutate(group = "child"), on.carrier %>% mutate(group = "carrier")) %>% 
  select(Disorder_Group, group, count) %>%
  # 修复：过滤掉空值和NA值，并为空值提供默认名称
  filter(!is.na(Disorder_Group) & Disorder_Group != "" & !is.null(Disorder_Group)) %>%
  mutate(Disorder_Group = ifelse(is.na(Disorder_Group) | Disorder_Group == "", "Unknown", Disorder_Group)) %>%
  pivot_wider(names_from = Disorder_Group, values_from = count, values_fill = list(count = 0)) %>% 
  unique()

child_data <- wide_data %>% filter(group == "child") %>% select(-group)
carrier_data <- wide_data %>% filter(group == "carrier") %>% select(-group)
colnames(child_data) <-  str_wrap(colnames(child_data), width = 20) 
colnames(carrier_data) <-  str_wrap(colnames(carrier_data), width = 20) 

max_count <- max(c(max(carrier_data), max(child_data)), na.rm = TRUE)
plot_data <- rbind(
  rep(max_count, ncol(carrier_data)), 
  rep(0, ncol(carrier_data)), 
  carrier_data,
  child_data
) 

na.col <- which(colnames(plot_data) == "NA")
if(length(na.col)> 0 ){
  plot_data <- plot_data[,-na.col]
}

# Red for positive and grey for carrier status 
colors <- c(rgb(0.541, 0.549, 0.749, 0.5), rgb(0.88, 0.44, 0.4,0.7)) # carrier, child
lines <- c( rgb(0.541, 0.549, 0.749, 1), rgb(0.88, 0.44, 0.4,1))

# positive and carrier onotologies 
png(filename = paste0(output_dir, "Figures/Positive_Carrier_Ontology.png"),width = 28,height = 28,units = "cm",res = 600)
radarchart(
  plot_data,
  axistype = 1,
  pcol = lines,          # Line color 
  pfcol = colors,         # Fill color 
  plwd = 2,               # Line width
  plty = 1,               # Line type 
  cglcol = "darkgrey",        # Grid color 
  cglty = 1,              # Grid type 
  caxislabels=seq(0,max_count, round(max_count/4)), # Customize the labels
  calcex = 0.8, # Axis value label size
  cglwd = 0.5,            # Grid width 
  axislabcol = "black",    # Axis color
  vlcex = 1,           # Text size 
  title = "Disease ontology distribution of positive Monogenic Diseases and Carrier Status"
)

legend("topright",
       legend = c("Carrier Status", "Monogenic Positives"),  
       col =  lines,
       pch = 20,
       pt.bg = colors,
       bty = "n",
       cex = 0.8
)

dev.off()

# Carrier status only 
plot_data_carrier <- rbind(
  rep(max(carrier_data), ncol(carrier_data)), # max value
  rep(0, ncol(carrier_data)),   # min value 
  carrier_data
) 

png(filename = paste0(output_dir, "Figures/Carrier_Ontology.png"),width = 28,height = 28,units = "cm",res = 600)
radarchart(
  plot_data_carrier,
  axistype = 1,
  pcol = rgb(0.541, 0.549, 0.749, 0.5),          # Line color 
  pfcol = rgb(0.541, 0.549, 0.749, 0.5),         # Fill color 
  plwd = 2,               # Line width
  plty = 1,               # Line type 
  cglcol = "darkgrey",        # Grid color 
  cglty = 1,              # Grid type 
  caxislabels=seq(0, max(carrier_data), round(max(carrier_data)/4)), # Customize the labels
  calcex = 0.8, # Axis value label size
  cglwd = 0.5,            # Grid width 
  axislabcol = "black",    # Axis color
  vlcex = 1,           # Text size 
  title = "Disease ontology distribution of Diseases Carrier Status"
)

legend("topright",
       legend = "Carrier Status",  
       col =  rgb(0.541, 0.549, 0.749, 0.5),
       pch = 20,
       pt.bg = rgb(0.541, 0.549, 0.749, 0.5),
       bty = "n",
       cex = 0.8
)

dev.off()


# Monogenic disease only 
plot_data_monogenic <- rbind(
  rep(max(child_data), ncol(child_data)),  # max value
  rep(0, ncol(child_data)),    # min value 
  child_data
) 

png(filename = paste0(output_dir, "Figures/Positive_Ontology.png"),width = 28,height = 28,units = "cm",res = 600)
radarchart(
  plot_data_monogenic,
  axistype = 1,
  pcol = rgb(0.88, 0.44, 0.4,0.7),          # Line color 
  pfcol = rgb(0.88, 0.44, 0.4,0.7),         # Fill color 
  plwd = 2,               # Line width
  plty = 1,               # Line type 
  cglcol = "darkgrey",        # Grid color 
  cglty = 1,              # Grid type 
  caxislabels=seq(0,max(child_data), max(child_data)/4), # Customize the labels
  calcex = 0.8, # Axis value label size
  cglwd = 0.5,            # Grid width 
  axislabcol = "black",    # Axis color
  vlcex = 2,           # Text size 
  title = "Disease ontology distribution of Positive Monogenic Disease"
)

legend("topright",
       legend = "Monogenic Positives",  
       col =  rgb(0.88, 0.44, 0.4,0.7),
       pch = 20,
       pt.bg = rgb(0.88, 0.44, 0.4,0.7),
       bty = "n",
       cex = 0.8
)

dev.off()

## ---- PGx ----
# For PGx, bar plot and compare heatmap compare with other populations (only for the significant levels of drugs and pediatric) drugs 
pgx.all <- data.frame()

for(file in pgx_list){
  print(file)
  filename = basename(file)
  tpm <- read.csv(file, header = T, sep = "\t") %>% 
    filter(Specialty.Population == "Pediatric") %>% 
    filter(Level.of.Evidence %in% c("1A","2A","1B","2B"),
           Testing.Level %in% c("Testing Recommended","Testing Required","Actionable PGx")) %>% 
    unique %>% select(-Zygous) %>% 
    mutate(Sample = filename)
  pgx.all <- rbind(pgx.all, tpm)
}

# bar plot 
individual.pgx <- pgx.all %>% 
  group_by(Sample) %>% 
  summarise(`No. of variant` = n_distinct(Variant.Haplotypes),.groups = "drop") %>% 
  group_by(`No. of variant`) %>% 
  summarise(`No. of sample` = n_distinct(Sample),.groups = "drop") 

write.table(individual.pgx, 
            file = paste0(output_dir, "/PGx_samples_variants_number.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

p1 <- ggplot(individual.pgx, aes(x = `No. of variant`, y = `No. of sample`)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "No. of PGx Haplotypes", y = "No. of Samples") +
  scale_x_continuous(breaks = seq(min(individual.pgx$`No. of variant`), max(individual.pgx$`No. of variant`), by = 1)) +
  theme_classic()

p1.w <- length(unique(individual.pgx$`No. of sample`))
p1.h <- length(unique(individual.pgx$`No. of variant`))

ggsave(p1, filename = paste0(output_dir, "Figures/PGx_samples_variants_number.png"),
       height = p1.w, width = p1.h,dpi= 600, units = "cm")

pgx.stat <- pgx.all %>% 
  group_by(Variant.Haplotypes) %>% 
  summarise(count = n_distinct(Sample)) %>% 
  arrange(desc(count)) %>% rename(`No. of samples` = count) 
write.table(pgx.stat, 
            file = paste0(output_dir, "/PGx_statistics.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

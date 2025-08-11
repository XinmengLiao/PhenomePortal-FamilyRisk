library(kinship2)
library(tidyr)
library(dplyr)
library(ggplot2)
for (i in c("select","filter", "mutate","rename", "left_join", "slice")){
  conflicted::conflict_prefer(i, "dplyr")
}
rm(i)
conflicted::conflicts_prefer(stats::sd)
conflicted::conflicts_prefer(httr::content)
conflicted::conflicts_prefer(plotly::layout)

## ---- 1. Load files ----
#user input
args <- commandArgs(trailingOnly = TRUE)
sampleid <- args[1]
input_file <- args[2]
ped_file <- args[3]
output_dir = args[4]
print(paste0("Now pocessing ", sampleid))
print(paste0("Input file: ", input_file))
print(paste0("Ped file: ", ped_file))
print(paste0("Output directory: ", output_dir))

# load files
result <- read.csv(input_file,header = T,sep = "\t") %>% unique()
colnames(result) <- gsub("am_","AlphaMissense_", colnames(result))
ped <- read.csv(ped_file, header = F,sep = "\t") %>% unique()
colnames(ped) <- c("FID","Individual","Father","Mother","Sex","Genotype")
familyid <- ped$FID %>% unique()

kid_genotype_col <- colnames(result)[grep("Kid.*Genotype",colnames(result))]
kid_pattern_col <- colnames(result)[grep("Kid.*Pattern",colnames(result))]
kidid <- gsub("Kid.*Genotype_", "", kid_genotype_col)
kidcount = length(kidid)
kid_cols <- c(colnames(result)[grep("Kid.*",colnames(result))])

# create detailed variant table from result
detailed_result <- result %>% select(-variant_info,-INFO)
write.table(detailed_result, paste0(output_dir,"/", sampleid, "_detailed_table.txt"),
            sep = "\t", quote = F, row.names = F )

# create family information table from ped 
family_info <- ped %>% 
  mutate(Relationship = case_when((Father == 0 & Mother == 0 & Sex == 2) ~ "Mother",
                                  (Father == 0 & Mother == 0 & Sex == 1) ~ "Father", 
                                  (Father != 0 & Mother != 0) ~ paste0("Kid", cumsum(Father != 0 & Mother != 0)), T ~ NA),
         Sex = if_else(Sex == 2, "Female", "Male"))
write.table(family_info, paste0(output_dir,"/", sampleid, "_family_information.txt"), 
            sep = "\t", quote = F, row.name = F )

df <- result %>% mutate(FamilyID = familyid) %>% mutate(FatherGender = "Male", MotherGender = "Female")
df1 <- df %>% 
  pivot_longer(
    cols = all_of(kid_cols),
    names_to = "Text",
    values_to = "Value",
    values_drop_na = FALSE
  )

if (unique(c("LoF","LoF_filter","LoF_flags","LoF_info") %in% names(df1))){
  del_cols <- c(names(df1)[grep("P001_*",names(df1))],"LoF","LoF_filter","LoF_flags","LoF_info")  
  df1 <- df1 %>% select(-all_of(del_cols))
}else{
  del_cols <- c(names(df1)[grep("P001_*",names(df1))])  
  df1 <- df1 %>% select(-all_of(del_cols))
}


# for ClinVar P/LP and acmg P/LP
clinvar.df <- df1 %>% filter(grepl("Pathogenic|Likely_pathogenic",ClinVar_CLNSIG) & 
                        grepl("Pathogenic|Likely_pathogenic",acmg_classification)) 

write.table(clinvar.df, paste0(output_dir,"/", sampleid, "_plp.txt"),
            sep = "\t", quote = F, row.names = F )

brief <- function(df){
	index <- grep("Kid", colnames(df))
	kid_pattern_col <- index[grep("Pattern", colnames(df)[index])]
	kid_genotype_col <- index[grep("Genotype", colnames(df)[index])]
	kid_cols_to_keep <- colnames(df)[c(kid_genotype_col, kid_pattern_col)]
    
	brief <- df %>%
		filter(Inheritance != "Unknown") %>%
      		filter(rowSums(!is.na(select(., all_of(kid_pattern_col)))) > 0) %>%
      		mutate(
        		ref = sapply(strsplit(variant_info, split = "_"), `[`, 3),
        		alt = sapply(strsplit(variant_info, split = "_"), `[`, 4),
        		FatherGenotype = mapply(function(fg, r, a) gsub("1", a, gsub("0", r, fg)), FatherGenotype, ref, alt),
        		MotherGenotype = mapply(function(mg, r, a) gsub("1", a, gsub("0", r, mg)), MotherGenotype, ref, alt),
        	across(all_of(kid_genotype_col),
               	~ mapply(function(kid, r, a) gsub("1", a, gsub("0", r, kid)), ., ref, alt),
               .names = "{.col}"),
        ClinVar_CLNSIG = gsub("_", " ", ClinVar_CLNSIG),
        acmg_classification = gsub("_", " ", acmg_classification),
        ClinVar_CLNSIG = sapply(strsplit(split = "&", ClinVar_CLNSIG), `[`, 1),
        ClinVar_CLNSIG = if_else(ClinVar_CLNSIG == "Pathogenic/Likely pathogenic/Pathogenic",
                                 "Pathogenic/Likely pathogenic",
                                 as.character(ClinVar_CLNSIG)),
        Existing_variation = sapply(strsplit(split = "&", Existing_variation), `[`, 1)
      ) %>%
      		select(Disease, Genes, Inheritance, variant_info, Existing_variation, ClinVar_CLNSIG, acmg_classification, FatherGenotype, MotherGenotype, all_of(kid_cols_to_keep)) %>%
      		rename(`Existing variant` = Existing_variation, ClinVar = ClinVar_CLNSIG, ACMG = acmg_classification) %>% distinct()
	return(brief)
}

clinvar.brief <- brief(clinvar.df)
write.table(clinvar.brief, paste0(output_dir,"/", sampleid, "_plp_brief.txt"),
            sep = "\t", quote = F, row.name = F )

# for Predicted high variants
high.df <- df1 %>% filter(IMPACT =="HIGH" & ClinVar_CLNSIG == "") 
high.brief <- brief(high.df)
write.table(clinvar.brief, paste0(output_dir,"/", sampleid, "_high_brief.txt"),
            sep = "\t", quote = F, row.name = F )

# for other Predicted variants
predict.df <- df1 %>% filter(ClinVar_CLNSIG == "" ) %>% unique()
predict.brief <- brief(predict.df)
write.table(predict.brief, paste0(output_dir,"/", sampleid, "_predict_brief.txt"),
            sep = "\t", quote = F, row.name = F )

## ---- 2. ClinVar Variants by family ----
# recessive 
clinvar.re <- clinvar.df %>% 
  filter(Value == "recessive") %>% select(variant_info, FamilyID,Text, Value) %>% unique() %>% mutate(Type = "P/LP") 

# de novo
clinvar.denovo <- clinvar.df %>% 
  filter(Value == "de novo") %>% select(variant_info, FamilyID,Text, Value) %>% unique() %>% mutate(Type = "P/LP") 

# dominant
clinvar.do <- clinvar.df %>% 
  filter(grepl("dominant", Value)) %>% select(variant_info, FamilyID,Text, Value) %>% unique() %>% mutate(Type = "P/LP") 


## ---- 3. Predicted HIGH Variants by family ----
# recessive 
high.re <- high.df %>% 
  filter(Value == "recessive") %>% select(variant_info, FamilyID,Text, Value) %>% unique() %>% mutate(Type = "pLoF") 

# de novo
high.denovo <- high.df %>% 
  filter(Value == "de novo") %>% select(variant_info, FamilyID,Text, Value) %>% unique() %>% mutate(Type = "pLoF") 

# dominant
high.do <- high.df %>% 
  filter(grepl("dominant", Value)) %>% select(variant_info, FamilyID,Text, Value) %>% unique() %>% mutate(Type = "pLoF") 


## ---- 4. Predicted missense variants by family ----
# recessive 
missense.re <- predict.df %>% filter(grepl("^missense_variant", Consequence)) %>% 
  filter(Value == "recessive") %>% select(variant_info, FamilyID,Text, Value) %>% unique() %>% mutate(Type = "pMissense") 

# de novo
missense.denovo <- predict.df %>% filter(grepl("^missense_variant", Consequence)) %>% 
  filter(Value == "de novo") %>% select(variant_info, FamilyID,Text, Value) %>% unique() %>% mutate(Type = "pMissense") 

# dominant
missense.do <- predict.df %>% filter(grepl("^missense_variant", Consequence)) %>% 
  filter(grepl("dominant", Value)) %>% select(variant_info, FamilyID,Text, Value) %>% unique() %>% mutate(Type = "pMissense") 


## ---- 5. Other predicted variants by family ----
# recessive 
other.re <- predict.df %>% filter(!grepl("^missense_variant", Consequence) & IMPACT != "HIGH") %>% 
  filter(Value == "recessive") %>% select(variant_info, FamilyID,Text, Value) %>% unique() %>% mutate(Type = "Other") 

# de novo
other.denovo <- predict.df %>% filter(!grepl("^missense_variant", Consequence) & IMPACT != "HIGH") %>% 
  filter(Value == "de novo") %>% select(variant_info, FamilyID,Text, Value) %>% unique() %>% mutate(Type = "Other") 

# dominant
other.do <- predict.df %>% filter(!grepl("^missense_variant", Consequence) & IMPACT != "HIGH") %>% 
  filter(grepl("dominant", Value)) %>% select(variant_info, FamilyID,Text, Value) %>% unique() %>% mutate(Type = "Other") 


## ---- 6. variant bar plot ----
recessive.data <- rbind(clinvar.re,high.re, missense.re,other.re) %>% mutate(Group = "Recessive")
denovo.data <- rbind(clinvar.denovo,high.denovo, missense.denovo, other.denovo) %>% mutate(Group = "de novo")
dominant.data <- rbind(clinvar.do, high.do, missense.do,other.do) %>% mutate(Group = "Dominant")
data <- rbind(recessive.data, denovo.data,dominant.data)

for(i in unique(data$Text)){
  print(paste0("Now ploting variant statistics for ", i))
  tpm <- data %>%
    filter(Text == i) %>% unique() %>% 
    mutate(Type = if_else(is.na(Type), "Other", Type)) %>% 
    group_by(FamilyID, Type, Group) %>%
    summarise(count = n(), .groups = 'drop') %>% unique()
  tpm$Type = factor(tpm$Type, levels = c("P/LP", "pLoF", "pMissense", "Other"))
  tpm$Group = factor(tpm$Group, levels = c("Dominant","Recessive", "de novo"))
  
  p1 <- ggplot(tpm, aes(x = Group, y = count, fill = Type))+
    geom_bar(position = 'stack',stat = 'identity')+
    theme_classic() +
    scale_fill_manual(values = c(
      "pMissense" = "#7fcdbb",  # 薄荷青绿色，柔和明亮
      "pLoF"      = "#fdbb84",  # 柑橘浅橙色，温暖又有活力
      "P/LP"      = "#c2a5cf",  # 柔紫薰衣草，科学感柔和
      "Other"     = "#a6bddb"   # 浅蓝灰，保持中性平衡
    ))+
    labs(x = NULL, y = "No. of Variant per Child", fill = "Variant type")+
    theme(axis.text.x = element_text(angle = 0))
  
  height <- sum(tpm$count)
  
  ggsave(p1, filename=paste0(output_dir,"/", i, "_Variant_statistics.png"), 
         width = height*1.5, height = height, units = "cm",dpi = 600)
  
}


## ---- 7. Pedigree plot ----
clinvar_plp_var <- clinvar.df %>% pull(variant_info) %>% unique()

for ( i in clinvar_plp_var){
  variant = i
  print(paste0("Now ploting pedigree plot for ", i))
  
  row <- result %>% filter(variant_info == variant)
  ped_data <- ped
  transcript <- unlist(strsplit(row$HGVSc,split = ":"))[2]
  protein <- unlist(strsplit(row$HGVSp,split = ":"))[2]
  gene <- row$SYMBOL %>% unique()
  kid_gender_col <- colnames(result)[grep("Kid.*Gender",colnames(result))]
  
  for (k in 1:length(kid_genotype_col)) {
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
  
  png(paste0(output_dir,"/", familyid, "_", gene,"_",transcript, "_pedigree_plot.png"), 
      width=800, height=900,res = 300)
  plot(ped1, col=ifelse(ped_data$avail, 2, 1), cex=0.6)
  title(main="Pedigree analysis",cex.main = 0.7)
  mtext(paste(gene, transcript, paste0("(", protein, ")")), 
        side=3, line=0.5, cex=0.5)
  dev.off()
}



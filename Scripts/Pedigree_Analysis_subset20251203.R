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
result <- args[2]
ped_file <- args[3]
output_dir = args[4]
variant_info_file <- args[2]

sampleid = "rwgs_F1"
result = "/Users/xinmengliao/Documents/Project/20250710_NewbornRisk/rwgs/rwgs_F1.txt"
ped_file = "/Users/xinmengliao/Documents/Project/20250710_NewbornRisk/rwgs/rwgs_F1.ped"
output_dir = "/Users/xinmengliao/Documents/Project/20250710_NewbornRisk/rwgs"
variant_info_file = "/Users/xinmengliao/Documents/Project/20250710_NewbornRisk/rwgs/variant_info.txt"

print(paste0("Now pocessing ", sampleid))
print(paste0("Result file: ", result))
print(paste0("Ped file: ", ped_file))
print(paste0("Output directory: ", output_dir))
print(paste0("variant_info_file: ", variant_info_file))


# read files
result <- read.csv(result,header = T,sep = "\t") %>% unique()
variant_info_df <- read.csv(variant_info_file,header = F,sep = "\t") %>% unique()
colnames(variant_info_df) <- "variant_info"
ped <- read.csv(ped_file, header = F,sep = "\t") %>% unique()
colnames(ped) <- c("FID","Individual","Father","Mother","Sex","Genotype")
familyid <- ped$FID %>% unique()


## ---- 2. Pedigree plot Functions ----

plot_pedigree_subset <- function(df, variant_info, ped){
  
  for ( i in variant_info){
    variant = i
    print(paste0("Now ploting pedigree plot for ", i))
    
    row <- df %>% filter(variant_info == variant)
    ped_data <- ped
    transcript <- if_else(is.na(unlist(strsplit(row$HGVSc,split = ":"))[2]), "", unlist(strsplit(row$HGVSc,split = ":"))[2])
    protein <- if_else(is.na(unlist(strsplit(row$HGVSp,split = ":"))[2]), "", paste0("(",unlist(strsplit(row$HGVSp,split = ":"))[2],")"))
    gene <- row$SYMBOL %>% unique()
    kid_gender_col <- colnames(df)[grep("Kid.*Gender",colnames(df))]
    kid_genotype_col <- colnames(df)[grep("Kid.*Genotype",colnames(df))]
    
    for (k in 1:length(kid_gender_col)) {
      kidid <- unlist(strsplit(kid_gender_col[k], split = "Gender_"))[2]
      genotype = row[k,kid_genotype_col]
      ped_data$Genotype[which(ped_data$Individual == kidid)] = genotype
    }
    
    # for father and mother genotype 
    f_genotype = row[1,"FatherGenotype"]
    m_genotype = row[1,"MotherGenotype"]
    ped_data$Genotype[which(ped_data$Father ==0 & ped_data$Mother == 0 & ped_data$Sex == 1)] = f_genotype
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
    
    png(paste0(output_dir,"/", familyid, "_", gene,"_", transcript, "_pedigree_plot.png"), 
        width=800,  height=900,res = 300)
    plot(ped1, col=ifelse(ped_data$avail, 2, 1), cex=0.6)
    #title(main="Pedigree analysis",cex.main = 0.7)
    mtext(paste(gene, transcript, protein), side=3, line=0.5, cex=0.5)
    dev.off()
  }
}


## ---- 4. Output ClinVar variants with different patterns to show as default ----

plot_pedigree_subset(df = result, variant_info = variant_info_df$variant_info, ped = ped)


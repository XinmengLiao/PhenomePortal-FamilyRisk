## ---- Load libraries ----
library(dplyr)
library(ggplot2)

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
pop_file <- args[3]
output_dir <- args[4]
cohort <- args[5] # family or population

# for family, the family PGS will be compared with the entire reference populations;
# for population, the population PGS will be compared with each of the super-population. 

result_file <- "/Users/xinmengliao/Documents/Project/20250710_NewbornRisk/Batch/newborn103-asthma-PGS002248-score/newborn103_pgs.txt"
pop_file <- '/Users/xinmengliao/Documents/Project/20250710_NewbornRisk/Batch/newborn103-asthma-PGS002248-score/newborn103_popsimilarity.txt'
output_dir <- '/Users/xinmengliao/Documents/Project/20250710_NewbornRisk/Batch/Results20251203'
cohort <- "population"

## ---- Load files ----
result <- read.csv(result_file,header = T,sep = "\t") 
pop <- read.csv(pop_file,header = T,sep = "\t") %>% select(sampleset, FID, IID, SuperPop, Group) %>% unique()
result <- result %>% left_join(., pop)

## ---- Family ----
if(cohort == "family"){
  
  # output the cohort PGS score
  df <- result %>% 
    filter(sampleset != "reference") %>% unique()
  write.table(df, paste0(output_dir, "/cohort_PGS.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  
  # draw the density plot 
  pop_colors <- c(
    "Cohort" = "#F4A7B9",
    "1KG+HGDP" = "#6A8DBF"
  )
  
  plot.df <- result %>% 
    mutate(Group = ifelse(sampleset != "reference", "Cohort", "1KG+HGDP"))
  
  p1 <- ggplot(plot.df, aes(x = Z_norm1, color = Group,fill = Group)) +
    geom_density(alpha = 0.15) +
    scale_color_manual(values = pop_colors) +
    scale_fill_manual(values = pop_colors) +
    theme_bw() +
    theme(plot.title = element_text(size = 20, face = "bold",hjust = 0.5),
          legend.position = "right", axis.title = element_text(size = 18), axis.text = element_text(size = 15),
          legend.text = element_text(size = 13), legend.title = element_text(size = 15),
          strip.text = element_text(size = 16, face = "bold")
    ) +
    labs(x = NULL, y = "Density", title = "Z_norm1")
  
  p2 <- ggplot(plot.df, aes(x = Z_norm2, color = Group,fill = Group)) +
    geom_density(alpha = 0.15) +
    scale_color_manual(values = pop_colors) +
    scale_fill_manual(values = pop_colors) +
    theme_bw() +
    theme(plot.title = element_text(size = 20, face = "bold",hjust = 0.5),
          legend.position = "right", axis.title = element_text(size = 18), axis.text = element_text(size = 15),
          legend.text = element_text(size = 13), legend.title = element_text(size = 15),
          strip.text = element_text(size = 16, face = "bold")
    ) +
    labs(x = NULL, y = "Density", title = "Z_norm2")
  
  ggsave(p1, filename = paste0(output_dir, "/PGS_Znorm1_DensityPlot.png"), width = 8, height = 6, dpi = 300)
  ggsave(p2, filename = paste0(output_dir, "/PGS_Znorm2_DensityPlot.png"), width = 8, height = 6, dpi = 300)
  
}

## ---- Population ----
if(cohort == "population"){
  
  # output the cohort PGS score
  df <- result %>% 
    filter(sampleset != "reference") %>% unique()
  write.table(df, paste0(output_dir, "/cohort_PGS.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  
  # draw the density plot 
  pop_colors <- c(
    "AFR" = "#e8c6a7",
    "AMR" = "#d2c8f0",
    "EAS" = "#c5e7da",
    "EUR" = "#c9e3ec",
    "MID" = "#f3c3d2",
    "SAS" = "#f0d8b4",
    "Cohort"  = "#2c3e50"  # 深蓝色
  )
  
  plot.df <- result %>% 
    mutate(SuperPop = ifelse(sampleset != "reference", "Cohort", SuperPop))
  
  p1 <- ggplot(plot.df, aes(x = Z_norm1, color = SuperPop,fill = SuperPop)) +
    geom_density(alpha = 0.15) +
    scale_color_manual(values = pop_colors) +
    scale_fill_manual(values = pop_colors) +
    theme_bw() +
    theme(plot.title = element_text(size = 20, face = "bold",hjust = 0.5),
      legend.position = "right", axis.title = element_text(size = 18), axis.text = element_text(size = 15),
      legend.text = element_text(size = 13), legend.title = element_text(size = 15),
      strip.text = element_text(size = 16, face = "bold")
    ) +
    labs(x = NULL, y = "Density", title = "Z_norm1")
  
  p2 <- ggplot(plot.df, aes(x = Z_norm2, color = SuperPop,fill = SuperPop)) +
    geom_density(alpha = 0.15) +
    scale_color_manual(values = pop_colors) +
    scale_fill_manual(values = pop_colors) +
    theme_bw() +
    theme(plot.title = element_text(size = 20, face = "bold",hjust = 0.5),
          legend.position = "right", axis.title = element_text(size = 18), axis.text = element_text(size = 15),
          legend.text = element_text(size = 13), legend.title = element_text(size = 15),
          strip.text = element_text(size = 16, face = "bold")
    ) +
    labs(x = NULL, y = "Density", title = "Z_norm2")
  
  ggsave(p1, filename = paste0(output_dir, "/PGS_Znorm1_DensityPlot.png"), width = 8, height = 6, dpi = 300)
  ggsave(p2, filename = paste0(output_dir, "/PGS_Znorm2_DensityPlot.png"), width = 8, height = 6, dpi = 300)
  
}

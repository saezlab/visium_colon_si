# Copyright (c) [2020] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#'Makes dotplots of differential expression analysis
#'

library(tidyverse)
library(Seurat)
library(clustree)
library(progeny)
library(dorothea)
library(cowplot)

source("visiumtools/reading2clustering.R")
source("visiumtools/funcomics.R")
source("visiumtools/differential_test.R")
source("visiumtools/misty_pipelines.R")
source("visiumtools/misty_utils.R")

sample_ids = set_names(c("V19S23-097_A1","V19S23-097_B1"))


# Fix dea results

for(slide in sample_ids){
  
  print(slide)
  
  dea_file = sprintf("results/single_slide/%s/%s_diff_features_ldvgann.rds",
                     slide,slide)
  
  dea_res = readRDS(dea_file)
  
  dea_res_fixed = lapply(dea_res,function(x){
    
    new_res = x %>% dplyr::mutate(cluster = as.numeric(as.character(cluster))) %>%
      arrange(cluster,p_val_adj)
    
    return(new_res)
      
  })
  
  saveRDS(dea_res_fixed, file = dea_file)
  
}


# Plot general features
for(slide in sample_ids){
  
  print(slide)
  
  slide_file = sprintf("./results/single_slide/%s/%s.rds",
                       slide,slide)
  
  visium_slide = readRDS(slide_file)
  
  plot_file = sprintf("./results/single_slide/%s/%s_funcomics_feat.pdf",
                      slide,slide)
  
  func_assays = c("progeny","ctscores")
  
  pdf(height = 12, width = 15,
      file = plot_file)
  
  for(a in Assays(visium_slide)[Assays(visium_slide) %in% func_assays]){
    
    DefaultAssay(visium_slide) = a
    
    plot(SpatialPlot(visium_slide,image.alpha = 0, 
                     features = rownames(visium_slide),
                     pt.size.factor = 6))
    
  }
  
  dev.off()
  
}

# First generate excel file with summaries:

for(slide in sample_ids){
  
  print(slide)
  
  dea_file = sprintf("results/single_slide/%s/%s_diff_features_ldvgann.rds",
                     slide,slide)
  
  dea_excel = sprintf("results/single_slide/%s/%s_diff_features_ldvgann.xlsx",
                     slide,slide)
  
  dea_res = readRDS(dea_file)
  
  diff_expr_xlsx(differential_features = dea_res,
                 only_pos = T,
                 excel_out = dea_excel)
}

# Generate dotplots

for(slide in sample_ids){
  
  print(slide)
  
  slide_file = sprintf("./results/single_slide/%s/%s.rds",
                       slide,slide)
  
  visium_slide = readRDS(slide_file)
  
  dea_file = sprintf("results/single_slide/%s/%s_diff_features_ldvgann.rds",
                     slide,slide)
  
  dea_res = readRDS(dea_file)
  
  dea_plots = sprintf("results/single_slide/%s/%s_diff_features_dplots_ldvgann.pdf",
                      slide,slide)
  
  dot_plots = get_dot_list(visium_slide = visium_slide,
               dea_res = dea_res,
               top_genes = 5,
               p_val_thrsh = 0.05)
  
  pdf(file = dea_plots, height = 15, width = 10)
  
  lapply(dot_plots,plot)
  
  dev.off()
  
}




# Copyright (c) [2020] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#'Generates visium objects + qc + funcomics + differential expression analysis

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

sample_ids = set_names(c("V19S23-097_A1","V19S23-097_B1",
               "V19S23-097_C1", "V19S23-097_D1"))

colon_slides = c("V19S23-097_A1","V19S23-097_B1")

marker_df_e = read_delim(file = "markers/Epithelial_Marker_genes_colon.csv", 
                       delim = ",") %>%
  arrange(cluster, p_val_adj) %>%
  select(gene, avg_logFC, cluster) %>%
  group_by(cluster) %>%
  slice(1:200) %>%
  ungroup() %>%
  mutate(cluster = paste0("Epithelial_c",cluster))

marker_df_s = read_delim(file = "markers/Stromal_Marker_genes_colon.csv", 
                         delim = ",") %>%
  arrange(cluster, p_val_adj) %>%
  select(gene, avg_logFC, cluster) %>%
  group_by(cluster) %>%
  slice(1:200) %>%
  ungroup() %>%
  mutate(cluster = paste0("Stromal_c",cluster))

marker_df = bind_rows(marker_df_e, marker_df_s)

for(slide in sample_ids){
  
  print(slide)
  
  slide_dir = sprintf("./results/single_slide/%s",slide)
  #system(paste0("mkdir results/single_slide/",slide))
  
  slide_out = sprintf("results/single_slide/%s/%s.rds",
                      slide,slide)
  
  slide_dea = sprintf("results/single_slide/%s/%s_diff_features.rds",
                      slide,slide)
  
  visium_slide = process_visium(dir_path = slide,
                                var_features = "seurat", 
                                p_value_thrsh = NULL,
                                out_dir = slide_dir,
                                resolution = 1,
                                verbose = FALSE)
  
  if(slide %in% colon_slides){

  visium_slide = add_funcomics(visium_slide = visium_slide,
                               species = "mouse",
                               confidence_lbls = c("A","B","C"),
                               top = 1000,
                               marker_df = marker_df,
                               verbose = FALSE)
  
  differential_features = find_allfeat(visium_slide = visium_slide,
                                       assays_collection = c("SCT",
                                                             "dorothea",
                                                             "progeny",
                                                             "ctscores"))
  
  }else{
    
    visium_slide = add_funcomics(visium_slide = visium_slide,
                                 species = "mouse",
                                 confidence_lbls = c("A","B","C"),
                                 top = 1000,
                                 marker_df = NULL,
                                 verbose = FALSE)
    
    differential_features = find_allfeat(visium_slide = visium_slide,
                                         assays_collection = c("SCT",
                                                               "dorothea",
                                                               "progeny"))
    
  }
  
  # Perform differential expression analysis of features in all assays
  # Warning: By default I am running Wilcox tests, because they seem
  # to work good enough for characterization purposes
  saveRDS(visium_slide,
          file = slide_out)
  
  saveRDS(differential_features, file = slide_dea)
  
}


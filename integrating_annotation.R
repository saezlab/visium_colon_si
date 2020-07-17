# Copyright (c) [2020] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Maps to integrated identities

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

colon_slides = tibble(slide_id = c("V19S23-097_A1","V19S23-097_B1"),
                      day = c("d0","d14"))

cell_clst_ann = read.csv(file = "./meta_data_integrated_analysis.tsv",header = T,
                row.names = 1, sep = "\t", stringsAsFactors = F) %>%
  rownames_to_column("spot_id") %>%
  dplyr::left_join(colon_slides) %>%
  mutate(spot_id = unlist(lapply(strsplit(spot_id,"_"), 
                                 function(x) x[1])))

for(slide in colon_slides$slide_id){
  
  print(slide)
  
  slide_file = sprintf("./results/single_slide/%s/%s.rds",slide,slide)
  
  clust_out = sprintf("./results/single_slide/%s/%s_clustering_int.pdf",slide,slide)
  
  slide_dea = sprintf("results/single_slide/%s/%s_diff_features_ldvgann.rds",
                    slide,slide)
  
  visium_slide = readRDS(slide_file)
  
  slide_ann = cell_clst_ann %>% 
    dplyr::filter(slide_id == slide)
  
  spot_ann = set_names((slide_ann$SCT_snn_res.0.8),
                      slide_ann$spot_id)
  
  spot_ann_2 = factor(spot_ann,levels = sort(unique(spot_ann),decreasing = F))
    
  visium_slide = AddMetaData(visium_slide,
                             spot_ann_2[colnames(visium_slide)],
                             "intgrtd_ident")
  
  Idents(visium_slide) = "intgrtd_ident"
  
  ident_cols = sample(cols,length(levels(Idents(visium_slide))))
  
  pdf(file = clust_out, height = 6,width = 10)
  
  p1 = DimPlot(visium_slide, reduction = "umap", label = FALSE) +
    scale_color_manual(values = ident_cols)
  p2 = Seurat::SpatialPlot(visium_slide,
                           image.alpha = 0,
                           pt.size.factor = 6,label = F) + 
    scale_fill_manual(values = ident_cols)
  
  plot(plot_grid(p1,p2,nrow = 1,align = "hv"))
  
  dev.off()
  
  differential_features = find_allfeat(visium_slide = visium_slide,
                                       assays_collection = c("SCT",
                                                             "dorothea",
                                                             "progeny",
                                                             "ctscores"))
  saveRDS(visium_slide,file = slide_file)
  
  saveRDS(differential_features,file = slide_dea)
  
}

































































# Copyright (c) [2020] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Factor analysis correlation with pathways

library(tidyverse)
library(Seurat)
library(cowplot)
library(corrplot)

source("visiumtools/reading2clustering.R")
source("visiumtools/funcomics.R")
source("visiumtools/differential_test.R")
source("visiumtools/misty_pipelines.R")
source("visiumtools/misty_utils.R")

get_order = function(cor_mat_feat_factor){
  cor_factors = cor(cor_mat_feat_factor)
  clust_factors = hclust(as.dist(1-cor_factors))
  
  cor_feats = cor(t(cor_mat_feat_factor))
  clust_feats = hclust(as.dist(1-cor_feats))
  
  return(list("factors" = clust_factors$labels[clust_factors$order],
              "features" = clust_feats$labels[clust_feats$order]))
}

# Slide 1 is A, slide 2 is B

all_samples = read.table(file = "./markers/NMF_Colon_combined.tsv",
           sep = "\t",header = T,row.names = 1)

A1_combined_factors = all_samples[grepl("1_1",rownames(all_samples)),]
A1_factors = read.table(file = "./markers/NMF_Colon_d0.tsv",
                        sep = "\t",header = T,row.names = 1)

B1_combined_factors = all_samples[grepl("1_2",rownames(all_samples)),]
B1_factors = read.table(file = "./markers/NMF_Colon_d14.tsv",
                        sep = "\t",header = T,row.names = 1)

slides = set_names(c("V19S23-097_A1","V19S23-097_B1"))

for(slide in slides){
  
  print(slide)
  
  slide_file = sprintf("./results/single_slide/%s/%s.rds",
                       slide,slide)
  
  path_file = sprintf("./results/single_slide/%s/pathway_factor_cor.pdf",
                      slide,slide)
  
  tf_file = sprintf("./results/single_slide/%s/tf_factor_cor.pdf",
                      slide,slide)
  
  visium_slide = readRDS(slide_file)
  DefaultAssay(visium_slide) = "SCT"
  
  PROGENy_mat = t(as.matrix(visium_slide@assays$progeny@data))
  TF_mat = t(as.matrix(visium_slide@assays$dorothea@data))
  
  #Correlations of factors vs funcomics
  
  if(slide == "V19S23-097_A1"){
   
    rownames(PROGENy_mat) = paste0(rownames(PROGENy_mat),"_1")
    rownames(TF_mat) = paste0(rownames(TF_mat),"_1")
    
    combined_cor_path = cor(PROGENy_mat[rownames(A1_combined_factors),],A1_combined_factors)
    combined_cor_TFs = cor(TF_mat[rownames(A1_combined_factors),],A1_combined_factors)
    
    single_cor_path = cor(PROGENy_mat[rownames(A1_factors),],A1_factors)
    single_cor_TFs = cor(TF_mat[rownames(A1_factors),],A1_factors)
    
    
  }else{
    
    rownames(PROGENy_mat) = paste0(rownames(PROGENy_mat),"_2")
    rownames(TF_mat) = paste0(rownames(TF_mat),"_2")
    
    combined_cor_path = cor(PROGENy_mat[rownames(B1_combined_factors),],B1_combined_factors)
    combined_cor_TFs = cor(TF_mat[rownames(B1_combined_factors),],B1_combined_factors)
    
    single_cor_path = cor(PROGENy_mat[rownames(B1_factors),],B1_factors)
    single_cor_TFs = cor(TF_mat[rownames(B1_factors),],B1_factors)
    
  }
  
  # PLOT correlation heatmaps
  
  # Pathways
  
  pdf(file = path_file, width = 10,height = 9)
  
  corrplot(combined_cor_path[get_order(combined_cor_path)$features,
                             get_order(combined_cor_path)$factors],
                method = "color",is.corr = F,
                tl.col = "black",title = "combined_factors",
           col=colorRampPalette(c("darkblue","white","darkred"))(100))
  
  corrplot(single_cor_path[get_order(single_cor_path)$features,
                                get_order(single_cor_path)$factors],
                method = "color",is.corr = F,
                tl.col = "black",title = "individual_factors",
           col=colorRampPalette(c("darkblue","white","darkred"))(100))
  
  dev.off()
  
  # TFs
  
  top_combined = names(sort(rowSums(abs(combined_cor_TFs)),decreasing = T)[1:25])
  top_individual = names(sort(rowSums(abs(single_cor_TFs)),decreasing = T)[1:25])
  
  pdf(file = tf_file, width = 10,height = 14)
  
  order_combined = get_order(combined_cor_TFs[top_combined,])
  order_individual = get_order(single_cor_TFs[top_individual,])
  
  corrplot(combined_cor_TFs[order_combined$features,
                            order_combined$factors],
                method = "color",is.corr = F,
                tl.col = "black",
                title = "combined_factors",
           col=colorRampPalette(c("darkblue","white","darkred"))(100))
  
  corrplot(single_cor_TFs[order_individual$features,
                          order_individual$factors],
           method = "color",is.corr = F,
           tl.col = "black",
           title = "individual_factors",
           col=colorRampPalette(c("darkblue","white","darkred"))(100))
  
  dev.off()
  
    
}
  
  





































# Copyright (c) [2020] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Run MISTy with ligands
#' 

library(OmnipathR)
library(tidyverse)
library(Seurat)

source("visium_exploratory/MISTy_ct_runs.R")
source("visiumtools/reading2clustering.R")
source("visiumtools/funcomics.R")
source("visiumtools/differential_test.R")
source("visiumtools/misty_pipelines.R")
source("visiumtools/misty_utils.R")

# Basic function to convert human to mouse gene names
convertHumanGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), 
                   filters = "hgnc_symbol", 
                   values = x , mart = human, 
                   attributesL = c("mgi_symbol"),
                   martL = mouse, uniqueRows=T)
  
  humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  return(genesV2)
}

lig_rec = import_intercell_network(
  interactions_param = list(datasets = c('ligrecextra', 'omnipath', 'pathwayextra'),
                            transmitter_param = list(parent = 'ligand'),
                            receiver_param = list(parent = 'receptor')
  ))

lig_rec_red = lig_rec %>% dplyr::ungroup() %>%
  dplyr::select(target_genesymbol,source_genesymbol,
                is_directed,category_intercell_target) %>%
  dplyr::filter(category_intercell_target == "receptor")


ligand_mgenes = convertHumanGeneList(unique(lig_rec_red$source_genesymbol))
colnames(ligand_mgenes) = c("source_genesymbol","source_mgenesymbol")
lig_rec_red = left_join(lig_rec_red, ligand_mgenes)
lig_rec_red = lig_rec_red %>% 
  dplyr::filter(!is.na(source_mgenesymbol))
ligs_m = unique(lig_rec_red$source_mgenesymbol)

#Cytokines

gene_sets = readRDS(file = "Genesets_Dec19.rds")
cytokines = gene_sets$MSIGDB_KEGG$KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION
cytokines_m = convertHumanGeneList(unique(cytokines))
cytokines_m = cytokines_m$MGI.symbol

# Main pipeline

slides = set_names(c("V19S23-097_A1","V19S23-097_B1"))

for(slide in slides){
  print(slide)
  
  misty_cytok_out = sprintf("./results/single_slide/%s/misty/%s_ctk",
                            slide,slide)
  
  misty_opath_out = sprintf("./results/single_slide/%s/misty/%s_opath",
                           slide,slide)
  
  #dea_file = sprintf("./results/single_slide/%s/%s_diff_features_ldvgann.rds",
  #                   slide,slide)
  
  #dea_res = readRDS(dea_file)
  
  slide_file = sprintf("./results/single_slide/%s/%s.rds",
                       slide,slide)
  
  visium_slide = readRDS(slide_file)
  
  gene_filter = VariableFeatures(visium_slide)
  
  slide_ctk = cytokines_m[cytokines_m %in% gene_filter]
  slide_ligs = ligs_m[ligs_m %in% gene_filter]
  
  # Run MISTy pipelines
  
  ls = c(2,5,10,20)
  
  # This is PROGENy high level analysis
  
  test_opath = lapply(ls,para_ppln_seurat,
                     visium_slide = visium_slide,
                     intra_assay = "progeny",
                     intra_features = NULL,
                     para_assay = "SCT",
                     para_features = slide_ligs,
                     spot_ids = NULL,
                     out_alias = misty_opath_out)
  
  test_cytok = lapply(ls,para_ppln_seurat,
                      visium_slide = visium_slide,
                      intra_assay = "progeny",
                      intra_features = NULL,
                      para_assay = "SCT",
                      para_features = slide_ctk,
                      spot_ids = NULL,
                      out_alias = misty_cytok_out)
  
}

# Lets optimize and summarize
# Watch out for the bug in get optimal, that's why
# I run it separately

for(slide in slides){
  
  print(slide)
  
  ls = c(2,5)
  
  misty_cytok_out = sprintf("./results/single_slide/%s/misty/%s_ctk",
                            slide,slide)
  
  misty_opath_out = sprintf("./results/single_slide/%s/misty/%s_opath",
                            slide,slide)
  print("path")
  get_optimal(out_dir_name = misty_opath_out,
              ls = ls)
  
  print("cytok")
  get_optimal(out_dir_name = misty_cytok_out,
              ls = ls)
  
}

# Let's plot the summaries of the best predicted features
# In this case, since I am predicting PROGENy, I can inclue everything.
# This is relevant for gene vs gene MISTy
for(slide in slides){
  
  print(slide)
  
  misty_cytok_out = sprintf("./results/single_slide/%s/misty/%s_ctk",
                            slide,slide)
  
  misty_opath_out = sprintf("./results/single_slide/%s/misty/%s_opath",
                            slide,slide)
  
  folders = paste0(c(misty_cytok_out,misty_opath_out),"_optim")
  
  #Here I plot all predicted genes, but you can play with the thrsholds
  #The p-value measures the difference of using a single view vs a multiview
  #When predicting many features, I use p_value < 0.15, R2_trsh = 0.05 and
  #importance_cut = 0.5
  plot_misty_bic(MISTy_out_folders = folders,
                 p_value_trsh = 1,
                 importance_cut = 2,
                 R2_trsh = 0.05)
  
}


# Some auxiliary plots to get a paraview matrix to plot

para_test = get_para_matrix(visium_slide,para_assay = "SCT",
                para_features = c("Il11ra1","Il10ra","Ccl21a"),
                l = 5)

visium_slide[['para']] = CreateAssayObject(data = para_test)

DefaultAssay(visium_slide) = "para"


plot(SpatialPlot(visium_slide,image.alpha = 0, 
                 features = c("Il11ra1","Il10ra","Ccl21a"),
                 pt.size.factor = 6))














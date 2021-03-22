# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Analyze cycling stem cells

library(tidyverse)
library(Seurat)
library(ggrepel)
source("visiumtools/differential_test.R")
source("visiumtools/misty_pipelines.R")
source("visiumtools/misty_utils.R")

#' Generate a module score matrix
#' 
#' @param seurat_object: seurat object
#' @param ms_regulon: named list of gene sets
#' @return matrix with module scores per cell
get_ms_matrix <- function(seurat_object, 
                          ms_regulon, 
                          assay = NULL) {
  # Dealing with weird set names ------------------
  names_vect <- gsub("[.]", "_", names(ms_regulon))
  names_vect <- gsub("-", "_", names_vect)
  
  tf_act_mat <- AddModuleScore(seurat_object,
                               features = ms_regulon,
                               name = paste0(names_vect,"__"),
                               assay = assay)
  
  tf_act_mat <- tf_act_mat@meta.data
  
  cell_ids <- rownames(tf_act_mat)
  calculated_regulons <- colnames(tf_act_mat)[grepl("__", colnames(tf_act_mat))]
  
  tf_act_mat <- tf_act_mat[,calculated_regulons, drop=FALSE]
  
  colnames(tf_act_mat) <- unlist(map(strsplit(colnames(tf_act_mat), split = "__"), 
                                     function(x) x[1]))
  
  rownames(tf_act_mat) <- cell_ids
  
  tf_act_mat <- t(as.matrix(tf_act_mat))
  
  return(tf_act_mat)
}

# Main -------------------------------------------------

stem_cells_mrkrs <- c("Hmgb2", "Ube2c", "Pclaf", "Stmn1",
                      "Top2a", "Tubb5", "Birc5", "Mki67",
                      "Cenpf", "Tuba1b", "Cenpa", "Ccdc34",
                      "Tmpo", "Cdca3", "Ccna2", "Cdk1",
                      "Nucks1", "Smc4", "Spc24", "Cdca8",
                      "Nusap1", "Racgap1", "Pbk", "Kif15",
                      "Mad2l1")

slide <- "V19S23-097_B1"

slide_file <- sprintf("./results/single_slide/%s/%s.rds",
                     slide,slide)

visium_slide <- readRDS(slide_file)

DefaultAssay(visium_slide) <- "SCT"

cell_type_ms <- get_ms_matrix(seurat_object = visium_slide,
                              ms_regulon = list("proliferative_stems" = stem_cells_mrkrs),
                              assay = "SCT")

visium_slide$stem_score <- scale(cell_type_ms[1,])
stem <- visium_slide$stem_score

progeny_mod <- rbind(as.matrix(visium_slide@assays$progeny@data),
                     stem)

visium_slide[['progeny_mod']] = CreateAssayObject(data = progeny_mod)

# Correlation with pathway activities ----------------------------------------
stem_score <- tibble("cell_type" = colnames(visium_slide), 
                     "stem_score" = cell_type_ms[1,])

progeny_scores <- as.data.frame(visium_slide@assays$progeny@data) %>%
  rownames_to_column("pathway") %>%
  pivot_longer(-pathway, 
               names_to = "cell_type",
               values_to = "path_score") %>%
  left_join(stem_score)
  
progeny_scores_cors <- progeny_scores %>%
  group_by(pathway) %>%
  nest() %>%
  mutate(cor_res = map(data, function(x) { 
    stest <- cor.test(x$path_score, x$stem_score)
    return(broom::tidy(stest))
    })) %>% 
  dplyr::select(-data) %>%
  unnest()

progeny_scores_cors <- progeny_scores_cors %>%
  ungroup() %>%
  mutate(p.value = ifelse(p.value == 0, min(p.value), p.value)) %>%
  mutate(corr_p.value = p.adjust(p.value)) %>%
  mutate(logpvalue = -log10(corr_p.value))

cor_plt <- ggplot(progeny_scores_cors, aes(x = estimate, 
                                y = logpvalue,
                                label = pathway)) +
  geom_point() +
  ggrepel::geom_text_repel() + 
  theme_classic() +
  xlab("Pearson Correlation") +
  ylab("-log10(corrected p-value)") +
  ylim(c(0,200)) +
  xlim(-1,1)

pdf(file = paste0("./results/single_slide/V19S23-097_B1/","_correlation_stemscore_progeny.pdf"), width = 5, height = 5)

plot(cor_plt)

dev.off()

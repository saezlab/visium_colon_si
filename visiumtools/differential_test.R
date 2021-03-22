# Copyright (c) [2020] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#'Performs wilcoxon test to identify differential features
#'of identities (by default from clustering) of selected slides

find_allfeat = function(visium_slide, assays_collection = c("SCT","dorothea","progeny",
                                                            "ctscores")){
  require(Seurat)
  require(purrr)
  require(tidyr)
  require(tibble)
  require(dplyr)
  
  possible_assays = set_names(assays_collection)
  
  diff_features = lapply(possible_assays, function(x){
    
    DefaultAssay(visium_slide) = x
    
    vis_wilcox_markers = FindAllMarkers(visium_slide, 
                                        logfc.threshold = 0.05,
                                        test.use = "wilcox")
  })
  
  return(diff_features)
  
} 


## Generates excel output with differential expression results

#' Generates excel output with differential expression results
#' @param differential_features: result list from find_allfeat
#' @param p_val_thrsh: pvalue filter
#' @param abs_lfc_thrsh: abs log fold change filter
#' @param only_pos: TRUE if only up-regulations are kept
#' @param excel_out: output file

diff_expr_xlsx = function(differential_features,
                          p_val_thrsh = 0.05,
                          abs_lfc_thrsh = 0,
                          only_pos = F,
                          excel_out = "default.xls",
                          top = 100){
  require(xlsx)
  
  # All assays with differential features
  
  a_names = names(differential_features)
  
  for(assay in a_names){
    
    res_df = differential_features[[assay]] %>%
      dplyr::filter(abs(avg_logFC) >= abs_lfc_thrsh,
             p_val_adj <= p_val_thrsh) %>% 
      arrange(cluster, -abs(avg_logFC))
    
    if(only_pos){
      res_df = res_df %>%
        dplyr::filter(avg_logFC > 0)
    }
    
    res_df = res_df %>%
      dplyr::group_by(cluster) %>%
      dplyr::slice(1:top) %>%
      ungroup() %>%
      as.data.frame()
    
    write.xlsx(res_df, file = excel_out, 
               row.names = F,
                 sheetName = assay, append = TRUE)
    
  }
  
  return(NULL)
}

#' Returns a list of dot plots to visualize differential expression analysis
#' @param visium_slide: seurat visium object with identities
#' @param dea_res: result list from find_allfeat
#' @param top_genes: how many top genes to include
#' @param possible_assays: assays to generate summary (should be in dea_res)
#' @param p_val_thrsh: pvalue filter od dea_res

get_dot_list = function(visium_slide, 
                        dea_res, 
                        top_genes = 5,
                        possible_assays = c("SCT","progeny","dorothea",
                                            "ECM","ctscores",
                                            "ctscores_match"),
                        p_val_thrsh = 0.05
){
  
  possible_assays = set_names(possible_assays)
  
  # Dot plot of top 5 per cluster
  
  possible_assays = possible_assays[possible_assays %in%
                                      Assays(visium_slide)]
  
  if(!is.infinite(top_genes)){
    plot_summary_data = lapply(dea_res[possible_assays], function(x){
      x %>% 
        dplyr::filter(avg_logFC > 0,
                      p_val_adj<=p_val_thrsh) %>%
        arrange(cluster,-avg_logFC)%>%
        group_by(cluster) %>%
        slice(1:top_genes)
    })
  } else{
    
    plot_summary_data = lapply(dea_res[possible_assays], function(x){
      x %>% 
        dplyr::filter(avg_logFC > 0,
                      p_val_adj<=p_val_thrsh) %>%
        dplyr::arrange(cluster,-avg_logFC)
    })
    
  }
  
  plot_list = lapply(names(plot_summary_data), function(x){
    
    if(nrow(plot_summary_data[[x]])>0){
      
      feat_list =  unique(plot_summary_data[[x]]$gene)
      
      DotPlot(visium_slide,features = feat_list,
              assay = x) +
        coord_flip() + theme(axis.text.y = element_text(size=7),
                             axis.text.x = element_text(size=8,angle = 90,
                                                        vjust = 0.5),
                             axis.title = element_blank(),
                             legend.text = element_text(size=8),
                             legend.title = element_text(size=8)) +
        ggtitle(x)
    }else{
      NULL
    }
    
  })
  
  return(plot_list)
  
}











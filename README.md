# The spatial transcriptomic landscape of the healing intestine following damage

Supporting code of section "Predictive algorithms revealed coordinated signaling pathways depending on location." from Parigi, Larsson, Das, et. al.

**Abstract**
The intestinal barrier is composed of a complex cell network defining highly compartmentalized and specialized structures. However, how transcriptional activity is spatially organized and how it adapts to fundamental biological processes, such as tissue regeneration, remains largely unexplored. Here, we devise spatial transcriptomics (ST) to characterize the transcriptomic landscape of the murine colon in steady state conditions and during mucosal healing following injury. We characterized a previously unappreciated molecular regionalization of the murine colon at steady state conditions, which dramatically changes during mucosal healing in the distal, but not proximal, colon. We identified eight spatially-organized transcriptional programs defining compartmentalized mucosal healing, and we identified regions with dominant wired pathways. Furthermore, we showed that decreased p53 activation defined areas with increased presence of proliferating epithelial stem cells. Finally, we used our resource to map transcriptomics modules associated with progression of experimental colitis, human intestinal cell transcriptomic profiles, clinical subclasses of ulcerative colitis (UC) patients and inflammatory bowel disease risk genes demonstrating that ST can be used to inform clinical practice. Overall, we provide a publicly available resource defining principles of transcriptomic regionalization of the colon during mucosal healing and a framework to develop and progress further hypotheses. 

**Methods**
Estimation of signalling pathway activities in Spatial Transcriptomics
For each slide, pathway activities of each spot were estimated with PROGENy (Schubert et al.) using the top 1,000 genes of each transcriptional footprint. Input expression matrices of each slide were normalized with sctransform implemented in Seurat 3.1.4.9 (Hafemeister and Satija). To further annotate the factors identified in the NNMF of the combined slides, we calculated the Pearson correlation of the spot-level scores of each individual factor with their corresponding spot-level pathway activity in each slide separately.

Integration of scRNAseq data with ST data set
To map the gene expression signature of proliferating stem cells from our independent dataset of intestinal epithelial cells on the d14 colonic tissue, we calculated module scores as implemented in Seurat’s function AddModuleScore (Stuart et al.; Tirosh et al.). Next, we calculated the Pearson correlation between each PROGENy pathway score and the stem cell score. P-values were corrected using a Benjamini-Hochberg procedure. The genes used to build the signature were the following: Hmgb2, Ube2c, Pclaf, Stmn1, Top2a, Tubb5, Birc5, Mki67, Cenpf, Tuba1b, Cenpa, Ccdc34, Tmpo, Cdca3, Ccna2, Cdk1, Nucks1, Smc4, Spc24, Cdca8, Nusap1, Racgap1, Pbk, Kif15, and Mad2l1.





# LISTA (LIver Spatio-Temporal Atlas)
## Codes used in LISTA project.

**Step 1:** Run 00.Stereo_seq_Matrix2SeuratObject-pipeline.R to create Seurat objects. You can run SAW analysis pipeline (https://github.com/STOmics/SAW) to create matrix or download the processed matrix from our database: https://db.cngb.org/stomics/lista/download/

**Step 2:** Run 00.cut_zonation_layer_and_pathway_module_score.R to split spots into 9 zonation layers. The pathway module score can be add simutaneously. 

**Step 3:** Run 00.Ligand_receptor_interaction_zonation_analysis.R to calculate interaction score of ligand & receptor pairs. The ligand receptor pairs used in our study was provided in the mouse_lr_pair.txt file. You can investigate all of them or a part.

**Step 4:** Run 00.Find_zonation_pathway_phyper.test.R to detect pathways enriched with zonation genes.

**Step 5:** Run 00.run_RCTD.R to calculate cell type projection score of scRNAseq on Stereo-seq data. You can find a detailed tutorial from RCTD offtial website: https://github.com/dmcable/spacexr

**Step 6:** Run 00.run_scenic.py to calculate gene regulatory network. You can find a detailed tutorial from SCENIC offtial website: https://pyscenic.readthedocs.io/en/latest/

**Step 7:** Run 00.run_hotspot.py to calculate gene coexpression modules base on their expression pattern. You can find a detailed tutorial from Hotspot offtial website: https://yoseflab.github.io/Hotspot/

**Step 8:** Run 01.TBL1XR1.bulkRNAseq.preprocessing.sh to get raw bulk RNAseq matrix of TBL1XR1 purturbation data.

**Step 9:** Run 01.TBL1XR1.bulkRNAseq.DEseq.r to get differential expressed genes of TBL1XR1 purturbation data.

**Step 10:** Run 02.ATAC_chipseq_preprocessing.sh to get chromatin modification regions of ATAC or Chip-seq data.

**Step 11:** Run 02.ATAC_chipseq_Motif_scan.r to get putative binding motifs of specific chromatin modification regions.

**Step 12:** Run 03.scRNAseq_clustering_scanpy.py to cluster scRNAseq datasets.

![image](https://github.com/haoshijie13/LISTA/assets/59014440/92db2bcd-39fd-4bbb-906c-ed2e4b0f0e5c)

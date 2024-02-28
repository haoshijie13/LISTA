# LISTA (LIver Spatio-Temporal Atlas)
## Codes used in LISTA project.

**Step 1:** Run Stereo_seq_Matrix2SeuratObject-pipeline.R to create Seurat objects. You can run SAW analysis pipeline (https://github.com/STOmics/SAW) to create matrix or download the processed matrix from our database: https://db.cngb.org/stomics/lista/download/

**Step 2:** Run cut_zonation_layer_and_pathway_module_score.R to split spots into 9 zonation layers. The pathway module score can be add simutaneously. 

**Step 3:** Run Ligand_receptor_interaction_zonation_analysis.R to calculate interaction score of ligand & receptor pairs. The ligand receptor pairs used in our study was provided in the mouse_lr_pair.txt file. You can investigate all of them or a part.

**Step 4:** Run Find_zonation_pathway_phyper.test.R to detect pathways enriched with zonation genes.

**Step 5:** Run run_RCTD.R to calculate cell type projection score of scRNAseq on Stereo-seq data. You can find a detailed tutorial from RCTD offtial website: https://github.com/dmcable/spacexr

**Step 6:** Run run_scenic.py to calculate gene regulatory network. You can find a detailed tutorial from SCENIC offtial website: https://pyscenic.readthedocs.io/en/latest/

**Step 7:** Run run_hotspot.py to calculate gene coexpression modules base on their expression pattern. You can find a detailed tutorial from Hotspot offtial website: https://yoseflab.github.io/Hotspot/

![image](https://github.com/haoshijie13/LISTA/assets/59014440/92db2bcd-39fd-4bbb-906c-ed2e4b0f0e5c)

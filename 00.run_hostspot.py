import os
import glob
import sys
import numpy as np
import pandas as pd
import hotspot
import matplotlib.pyplot as plt
import matplotlib.colors
import seaborn as sns
import pickle
# Load the counts and positions
counts_file = sys.argv[2]      # input the cell×gene matrix
pos_file = sys.argv[3]         # input the cell×spatial position matrix
OUTDIR = "."
NAME = sys.argv[1]             # output prefix

HS_RESULTS = ''.join([OUTDIR,"/",NAME,"_hs_results.p"])
LCZ = ''.join([OUTDIR, "/", NAME, "_lcz.p"])
MODULES = ''.join([OUTDIR, "/", NAME, "_modules.p"])
HOTSPOT = ''.join([OUTDIR, "/", NAME, "_hotspot.p"])
pos = pd.read_csv(pos_file, index_col=0)
counts = pd.read_csv(counts_file, index_col=0) # Takes a while, ~10min
# Align the indices
counts = counts.loc[:, pos.index]
barcodes = pos.index.values
# Swap position axes
# We swap x'=y and y'=-x to match the slides in the paper
pos = pd.DataFrame(
    {
        'X': pos.X,
        'Y': pos.Y,
    }, index=pos.index
)
num_umi = counts.sum(axis=0)
# Filter genes
#gene_counts = (counts > 0).sum(axis=1)
#valid_genes = gene_counts >= 50
#counts = counts.loc[valid_genes]


# Create the Hotspot object and the neighborhood graph
hs = hotspot.Hotspot(counts, model='normal', latent=pos)

hs.create_knn_graph(
    weighted_graph=False, n_neighbors=5,
)


hs_results = hs.compute_autocorrelations(jobs=20)

with open(HS_RESULTS, "wb") as f:
    pickle.dump(hs_results,f)


#select the genes with significant spatial autocorrelation
hs_genes = hs_results.index[hs_results.FDR < 0.05]

# Compute pair-wise local correlations between these genes
lcz = hs.compute_local_correlations(hs_genes, jobs=20)

with open(LCZ, "wb") as f:
    pickle.dump(lcz,f)


modules = hs.create_modules(
    min_gene_threshold=5, core_only=False, fdr_threshold=0.05
)

with open(MODULES, "wb") as f:
    pickle.dump(modules, f)


#with open(HOTSPOT, "wb") as f:
#    pickle.dump(hs,f)

results = hs.results.join(hs.modules)
results.to_csv("".join([sys.argv[1],"-Regulon2Gene.csv"]))

module_scores = hs.calculate_module_scores()
module_scores.to_csv("".join([sys.argv[1],"-module_score.csv"]))

plt.rcParams['figure.figsize'] = (15.0, 12.0)
hs.plot_local_correlations()
plt.savefig("".join([sys.argv[1],"-regulon_module_number.pdf"]), dpi = 600)

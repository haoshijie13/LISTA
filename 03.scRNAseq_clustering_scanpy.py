#!/usr/bin/env python
# coding: utf-8


import scanpy as sc
import os
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scanpy.external as sce


adata = sc.read('./01.qc/PHx.combine.h5ad')

adata.raw.var.index = pd.Index(adata.var['features'])
adata.var.index = pd.Index(adata.var['features'])
adata.obs['Library'] = adata.obs['split'].tolist()
adata.var['mt'] = adata.var_names.str.startswith('mt-') 
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, use_raw=True, log1p=False, inplace=True)

sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts','pct_counts_mt'],
             jitter=0.4, multi_panel=True)

bdata = adata[(adata.obs.n_genes_by_counts < 4000) &( adata.obs.n_genes_by_counts > 500 )& (adata.obs.total_counts < 20000) & (adata.obs.total_counts > 1000) & (adata.obs.pct_counts_mt < 10), :]
sc.pl.violin(bdata, ['n_genes_by_counts', 'total_counts','pct_counts_mt'],
             jitter=0.4, multi_panel=True)

sc.pp.highly_variable_genes(bdata, min_mean=0.0125, max_mean=3, min_disp=0.5)
bdata = bdata[:, bdata.var.highly_variable]
sc.pp.regress_out(bdata, ['total_counts', 'pct_counts_mt'], n_jobs = 50)
sc.pp.scale(bdata, max_value=10)
sc.tl.pca(bdata, svd_solver='arpack', use_highly_variable = True)
sce.pp.harmony_integrate(bdata, key = 'Library', basis='X_pca', adjusted_basis='X_pca_harmony')
sc.pp.neighbors(bdata, use_rep = 'X_pca_harmony', n_neighbors=10, n_pcs=40)
sc.tl.umap(bdata)
bdata.write_h5ad('PHx.combine_filter.h5ad')





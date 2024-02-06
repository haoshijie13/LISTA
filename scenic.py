import os
import glob
import pickle
import pandas as pd
import numpy as np
import sys

args = sys.argv
print(args[1])

from dask.diagnostics import ProgressBar
if __name__ == '__main__':
    ProgressBar = ProgressBar()
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2
from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell
import seaborn as sns
import dask
dask.config.set(num_workers=80)
if __name__ == '__main__':
        print(os.getcwd())
        ex_matrix = pd.read_table(args[1], sep='\t', header=0, index_col=0).T
        print(ex_matrix)
        DATA_FOLDER="./"
        RESOURCES_FOLDER="/hwfssz1/ST_MCHRI/P18Z10200N0350_tmp/PROJECT/MiguelLab/liver_zonation/SCENIC"
        DATABASE_FOLDER="/hwfssz1/ST_MCHRI/P18Z10200N0350_tmp/PROJECT/MiguelLab/liver_zonation/SCENIC/motif"
#       SCHEDULER="123.122.8.24:8786"
        DATABASES_GLOB = os.path.join(DATABASE_FOLDER, "mm10_*.mc9nr.feather")
        MOTIF_ANNOTATIONS_FNAME = os.path.join(RESOURCES_FOLDER, "motifs-v9-nr.mgi-m0.001-o0.0.tbl")
        MM_TFS_FNAME = os.path.join(RESOURCES_FOLDER, 'mm_mgi_tfs.txt')
        REGULONS_FNAME = ''.join([ DATA_FOLDER, args[1],"_regulons.p"])
        MOTIFS_FNAME = ''.join([ DATA_FOLDER, args[1], "_motifs.txt"])
        ADJACENCIES_FNAME = ''.join([DATA_FOLDER, args[1], "_adjacencies.txt"])
        MODULES_FNAME = ''.join([DATA_FOLDER, args[1], "_modules.p"])
        AUC_FNAME = ''.join([DATA_FOLDER, args[1], "_AUC.txt"])

        print(REGULONS_FNAME)
        print(MOTIFS_FNAME)

        tf_names = load_tf_names(MM_TFS_FNAME)
        db_fnames = glob.glob(DATABASES_GLOB)
        def name(fname):
            return os.path.splitext(os.path.basename(fname))[0]
        dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]

        print(dbs)
        adjancencies = grnboost2(expression_data=ex_matrix, tf_names=tf_names, verbose=True)
        adjancencies.to_csv(ADJACENCIES_FNAME, index=False, sep='\t')
        modules = list(modules_from_adjacencies(adjancencies, ex_matrix))
        with open(MODULES_FNAME, 'wb') as f:
            pickle.dump(modules, f)
        df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME)
        df.head()
        df.to_csv(MOTIFS_FNAME)
        regulons = df2regulons(df)
        with open(REGULONS_FNAME, 'wb') as f:
            pickle.dump(regulons, f)
        auc_mtx = aucell(ex_matrix, regulons, num_workers=5)
        auc_mtx.to_csv(AUC_FNAME, sep = "\t")

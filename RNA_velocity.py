# -*- coding: utf-8 -*-

# RNA velocity analysis has been done by my colleague: Houyu Zhang

#Version: 1.1
#Copyright (c) 2024 __CarlosLab@PKU__. All rights reserved.

import anndata
import numpy as np
import pandas as pd
pd.core.common.is_list_like = pd.api.types.is_list_like
from scipy import io
import scanpy as sc
import scvelo as scv
scv.set_figure_params()
scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
NUMBA_DISABLE_INTEL_SVML=1

#============================
# 1. Interbrain
#============================
fileprefix = 'Interbrain_scvelo_'
X = io.mmread(fileprefix + "counts.mtx")
adata = anndata.AnnData(X=X.transpose().tocsr())
cell_meta = pd.read_csv(fileprefix + "metadata.csv")
with open(fileprefix + "GeneNames.csv", 'r') as f:
    gene_names = f.read().splitlines()
adata.obs = cell_meta
adata.obs.index = adata.obs['barcode']
adata.var.index = gene_names
pca = pd.read_csv(fileprefix + "pca.csv")
pca.index = adata.obs.index
adata.obsm['X_pca'] = pca.to_numpy()
adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T
adata.write(fileprefix + '.h5ad')
adata = sc.read_h5ad(fileprefix + '.h5ad')
Mu1 = scv.read('./loom_file/Mu93BS_1.loom', cache=True)
Mu2 = scv.read('./loom_file/Mu93BS_2.loom', cache=True)
Mu3 = scv.read('./loom_file/Mu93BS_3.loom', cache=True)
LI1 = scv.read('./loom_file/MuLIBS_1.loom', cache=True)
LI2 = scv.read('./loom_file/MuLIBS_2.loom', cache=True)
LI3 = scv.read('./loom_file/MuLIBS_3.loom', cache=True)
barcodes = [bc.split(':')[1] for bc in Mu1.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_10' for bc in barcodes]
Mu1.obs.index = barcodes
barcodes = [bc.split(':')[1] for bc in Mu2.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_11' for bc in barcodes]
Mu2.obs.index = barcodes
barcodes = [bc.split(':')[1] for bc in Mu3.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_9' for bc in barcodes]
Mu3.obs.index = barcodes
barcodes = [bc.split(':')[1] for bc in LI1.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_8' for bc in barcodes]
LI1.obs.index = barcodes
barcodes = [bc.split(':')[1] for bc in LI2.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_7' for bc in barcodes]
LI2.obs.index = barcodes
barcodes = [bc.split(':')[1] for bc in LI3.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_6' for bc in barcodes]
LI3.obs.index = barcodes

Mu1.var_names_make_unique()
Mu2.var_names_make_unique()
Mu3.var_names_make_unique()
LI1.var_names_make_unique()
LI2.var_names_make_unique()
LI3.var_names_make_unique()
ldata = Mu1.concatenate([Mu2, Mu3,LI1,LI2,LI3])
scv.utils.clean_obs_names(adata)
scv.utils.clean_obs_names(ldata)

Interbrain_adata = scv.utils.merge(adata, ldata)
scv.pl.proportions(Interbrain_adata, groupby='orig.ident', save=fileprefix + 'proportions.pdf')
scv.pp.filter_and_normalize(Interbrain_adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(Interbrain_adata, n_pcs=30, n_neighbors=30)
scv.tl.velocity(Interbrain_adata, mode='stochastic')
scv.tl.velocity_graph(Interbrain_adata, n_jobs=20)
scv.pl.velocity_embedding_stream(Interbrain_adata, basis='umap', color=['seurat_clusters'], title='',dpi=500, density=5, size=100, linewidth=1, figsize=(8,8),
                                 palette = ["#F8766D","#00BA38"], save = fileprefix + 'embedding_stream.pdf')

scv.tl.velocity_confidence(Interbrain_adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(Interbrain_adata, c=keys, cmap='coolwarm', perc=[5, 95],size=80, figsize=(8,8),dpi=500, save = fileprefix + 'speed_scatter.pdf')
scv.tl.velocity_pseudotime(Interbrain_adata)
scv.pl.scatter(Interbrain_adata, color='velocity_pseudotime', cmap='gnuplot', size=80, figsize=(8,8),dpi=500, save = fileprefix + 'pseudotime.pdf')
scv.tl.recover_dynamics(Interbrain_adata)
scv.tl.velocity(Interbrain_adata, mode='dynamical')
scv.tl.velocity_graph(Interbrain_adata)
scv.pl.velocity_embedding_stream(Interbrain_adata, basis='umap')
scv.pl.velocity_embedding_stream(Interbrain_adata, basis='umap', color=['seurat_clusters'], title='', dpi=500, density=5, size=100, linewidth=1, figsize=(8,8),
                                 palette = ["#F8766D","#00BA38"], save = fileprefix + 'embedding_stream_dynamics.pdf')
scv.tl.latent_time(Interbrain_adata)
scv.pl.scatter(Interbrain_adata, color='latent_time', color_map='gnuplot', size=80, figsize=(8,8),dpi=500, save = fileprefix + 'latent_time.pdf')

#============================
# 2. Brainstem
#============================
fileprefix = 'Brainstem_scvelo_'
X = io.mmread(fileprefix + "counts.mtx")
adata = anndata.AnnData(X=X.transpose().tocsr())
cell_meta = pd.read_csv(fileprefix + "metadata.csv")
with open(fileprefix + "GeneNames.csv", 'r') as f:
    gene_names = f.read().splitlines()
adata.obs = cell_meta
adata.obs.index = adata.obs['barcode']
adata.var.index = gene_names
pca = pd.read_csv(fileprefix + "pca.csv")
pca.index = adata.obs.index
adata.obsm['X_pca'] = pca.to_numpy()
adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T
adata.write(fileprefix + '.h5ad')
adata = sc.read_h5ad(fileprefix + '.h5ad')
Mu1 = scv.read('./loom_file/Mu93BS_1.loom', cache=True)
Mu2 = scv.read('./loom_file/Mu93BS_2.loom', cache=True)
Mu3 = scv.read('./loom_file/Mu93BS_3.loom', cache=True)
LI1 = scv.read('./loom_file/MuLIBS_1.loom', cache=True)
LI2 = scv.read('./loom_file/MuLIBS_2.loom', cache=True)
LI3 = scv.read('./loom_file/MuLIBS_3.loom', cache=True)
barcodes = [bc.split(':')[1] for bc in Mu1.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_10' for bc in barcodes]
Mu1.obs.index = barcodes
barcodes = [bc.split(':')[1] for bc in Mu2.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_11' for bc in barcodes]
Mu2.obs.index = barcodes
barcodes = [bc.split(':')[1] for bc in Mu3.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_9' for bc in barcodes]
Mu3.obs.index = barcodes
barcodes = [bc.split(':')[1] for bc in LI1.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_8' for bc in barcodes]
LI1.obs.index = barcodes
barcodes = [bc.split(':')[1] for bc in LI2.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_7' for bc in barcodes]
LI2.obs.index = barcodes
barcodes = [bc.split(':')[1] for bc in LI3.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_6' for bc in barcodes]
LI3.obs.index = barcodes

Mu1.var_names_make_unique()
Mu2.var_names_make_unique()
Mu3.var_names_make_unique()
LI1.var_names_make_unique()
LI2.var_names_make_unique()
LI3.var_names_make_unique()
ldata = Mu1.concatenate([Mu2, Mu3,LI1,LI2,LI3])
scv.utils.clean_obs_names(adata)
scv.utils.clean_obs_names(ldata)

Brainstem_adata = scv.utils.merge(adata, ldata)

scv.pl.proportions(Brainstem_adata, groupby='orig.ident', save=fileprefix + 'proportions.pdf')
scv.pp.filter_and_normalize(Brainstem_adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(Brainstem_adata, n_pcs=30, n_neighbors=30)
scv.tl.velocity(Brainstem_adata, mode='stochastic')
scv.tl.velocity_graph(Brainstem_adata, n_jobs=20)
scv.pl.velocity_embedding_stream(Brainstem_adata, basis='umap', color=['seurat_clusters'], title='',dpi=500, density=5, size=100, linewidth=1, figsize=(8,8),
                                 palette = ["#F8766D","#00BA38"], save = fileprefix + 'embedding_stream.pdf')

scv.tl.velocity_confidence(Brainstem_adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(Brainstem_adata, c=keys, cmap='coolwarm', perc=[5, 95],size=80, figsize=(8,8),dpi=500, save = fileprefix + 'speed_scatter.pdf')
scv.tl.velocity_pseudotime(Brainstem_adata)
scv.pl.scatter(Brainstem_adata, color='velocity_pseudotime', cmap='gnuplot', size=80, figsize=(8,8),dpi=500, save = fileprefix + 'pseudotime.pdf')
scv.tl.recover_dynamics(Brainstem_adata)
scv.tl.velocity(Brainstem_adata, mode='dynamical')
scv.tl.velocity_graph(Brainstem_adata)
scv.pl.velocity_embedding_stream(Brainstem_adata, basis='umap')
scv.pl.velocity_embedding_stream(Brainstem_adata, basis='umap', color=['seurat_clusters'], title='', dpi=500, density=5, size=100, linewidth=1, figsize=(8,8),
                                 palette = ["#F8766D","#00BA38"], save = fileprefix + 'embedding_stream_dynamics.pdf')
scv.tl.latent_time(Brainstem_adata)
scv.pl.scatter(Brainstem_adata, color='latent_time', color_map='gnuplot', size=80, figsize=(8,8),dpi=500, save = fileprefix + 'latent_time.pdf')

#============================
# 3. Granule
#============================
fileprefix = 'Granule_scvelo_'
X = io.mmread(fileprefix + "counts.mtx")
adata = anndata.AnnData(X=X.transpose().tocsr())
cell_meta = pd.read_csv(fileprefix + "metadata.csv")
with open(fileprefix + "GeneNames.csv", 'r') as f:
    gene_names = f.read().splitlines()
adata.obs = cell_meta
adata.obs.index = adata.obs['barcode']
adata.var.index = gene_names
pca = pd.read_csv(fileprefix + "pca.csv")
pca.index = adata.obs.index
adata.obsm['X_pca'] = pca.to_numpy()
adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T
adata.write(fileprefix + '.h5ad') # save dataset as anndata format

adata = sc.read_h5ad(fileprefix + '.h5ad')
# load loom files for spliced/unspliced matrices for each sample:
Mu1 = scv.read('./loom_file/Mu93BS_1.loom', cache=True)
Mu2 = scv.read('./loom_file/Mu93BS_2.loom', cache=True)
Mu3 = scv.read('./loom_file/Mu93BS_3.loom', cache=True)
LI1 = scv.read('./loom_file/MuLIBS_1.loom', cache=True)
LI2 = scv.read('./loom_file/MuLIBS_2.loom', cache=True)
LI3 = scv.read('./loom_file/MuLIBS_3.loom', cache=True)

# rename barcodes in order to merge:
barcodes = [bc.split(':')[1] for bc in Mu1.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_10' for bc in barcodes]
Mu1.obs.index = barcodes
barcodes = [bc.split(':')[1] for bc in Mu2.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_11' for bc in barcodes]
Mu2.obs.index = barcodes
barcodes = [bc.split(':')[1] for bc in Mu3.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_9' for bc in barcodes]
Mu3.obs.index = barcodes
barcodes = [bc.split(':')[1] for bc in LI1.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_8' for bc in barcodes]
LI1.obs.index = barcodes
barcodes = [bc.split(':')[1] for bc in LI2.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_7' for bc in barcodes]
LI2.obs.index = barcodes
barcodes = [bc.split(':')[1] for bc in LI3.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_6' for bc in barcodes]
LI3.obs.index = barcodes

# make variable names unique
Mu1.var_names_make_unique()
Mu2.var_names_make_unique()
Mu3.var_names_make_unique()
LI1.var_names_make_unique()
LI2.var_names_make_unique()
LI3.var_names_make_unique()
ldata = Mu1.concatenate([Mu2, Mu3,LI1,LI2,LI3]) # merge matrices into the original adata object
scv.utils.clean_obs_names(adata)
scv.utils.clean_obs_names(ldata)
Granule_adata = scv.utils.merge(adata, ldata)

# plot umap to check
#sc.pl.umap(adata, color='Myannotation', frameon=False, legend_loc='on data', title='', save = fileprefix + 'celltypes.pdf')
scv.pl.proportions(Granule_adata, groupby='orig.ident', save=fileprefix + 'proportions.pdf')
scv.pp.filter_and_normalize(Granule_adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(Granule_adata, n_pcs=30, n_neighbors=30)
scv.tl.velocity(Granule_adata, mode='stochastic')
scv.tl.velocity_graph(Granule_adata, n_jobs=20)
scv.pl.velocity_embedding_stream(Granule_adata, basis='umap', color=['seurat_clusters'], title='',dpi=500, density=5, size=100, linewidth=1, figsize=(8,8),
                                 palette = ["#F8766D","#ABA300"],save = fileprefix + 'embedding_stream.pdf')

scv.tl.velocity_confidence(Granule_adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(Granule_adata, c=keys, cmap='coolwarm', perc=[5, 95],size=80, figsize=(8,8),dpi=500, save = fileprefix + 'speed_scatter.pdf')
scv.tl.velocity_pseudotime(Granule_adata)
scv.pl.scatter(Granule_adata, color='velocity_pseudotime', cmap='gnuplot', size=80, figsize=(8,8),dpi=500, save = fileprefix + 'pseudotime.pdf')
scv.tl.recover_dynamics(Granule_adata)
scv.tl.velocity(Granule_adata, mode='dynamical')
scv.tl.velocity_graph(Granule_adata)
scv.pl.velocity_embedding_stream(Granule_adata, basis='umap')
scv.pl.velocity_embedding_stream(Granule_adata, basis='umap', color=['seurat_clusters'], title='', dpi=500, density=5, size=100, linewidth=1, figsize=(8,8),
                                 palette = ["#F8766D","#ABA300"], save = fileprefix + 'embedding_stream_dynamics.pdf')
scv.tl.latent_time(Granule_adata)
scv.pl.scatter(Granule_adata, color='latent_time', color_map='gnuplot', size=80, figsize=(8,8),dpi=500, save = fileprefix + 'latent_time.pdf')


 















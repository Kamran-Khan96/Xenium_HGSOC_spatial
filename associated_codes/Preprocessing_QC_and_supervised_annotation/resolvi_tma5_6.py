import os
import tempfile
import pandas as pd
import anndata as ad
import numpy as np
import scanpy as sc
import scvi
import torch

scvi.settings.seed = 0

adata = sc.read(
    "/group/ms014/kkhan/Analysis/OVA_TMA/tma5_6_together_scVI/tma5_6_merged_tumour.h5ad"
)

if "spatial" in adata.obsm:
    adata.obsm["X_spatial"] = adata.obsm["spatial"].copy()

scvi.external.RESOLVI.setup_anndata(
    adata,
    labels_key="decoupler_fine",
    layer="raw_counts",
    batch_key="tma"
)

# Train model
supervised_resolvi = scvi.external.RESOLVI(adata, semisupervised=True)
supervised_resolvi.train(max_epochs=200,
			 early_stopping=True)

# Predictions + latent space
adata.obsm["resolvi_celltypes"] = supervised_resolvi.predict(
    adata, num_samples=3, soft=True
)

adata.obs["resolvi_predicted"] = (
    adata.obsm["resolvi_celltypes"].idxmax(axis=1)
)

adata.obsm["X_resolVI"] = supervised_resolvi.get_latent_representation(adata)

sc.pp.neighbors(adata, use_rep="X_resolVI")
sc.tl.umap(adata)

adata.write_h5ad(
    "/group/ms014/kkhan/Analysis/OVA_TMA/tma5_6_together_scVI/tma5_6_merged_tum_resolvi.h5ad"
)




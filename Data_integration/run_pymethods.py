# scarches env
import scanpy as sc
import scvi
import pyreadr
import anndata
import numpy as np
import pandas as pd
from scvi.model.utils import mde
from sklearn.metrics import classification_report
from scarches.models.scpoli import scPoli
import warnings
warnings.filterwarnings('ignore')

data_dir = "Jurkat_293t"
condition_key = 'Batch'
cell_type_key = 'CellType'
npcs = 50

early_stopping_kwargs = {
    "early_stopping_metric": "val_prototype_loss",
    "mode": "min",
    "threshold": 0,
    "patience": 20,
    "reduce_lr": True,
    "lr_patience": 13,
    "lr_factor": 0.1,
}

########################################## Run on HVGs
# prepare adata object
dataset = pyreadr.read_r("/home/server/zy/group_scripts/datasets_preparation/" + data_dir + "/count_dataset.rds")
meta = pyreadr.read_r("/home/server/zy/group_scripts/datasets_preparation/" + data_dir + "/meta.rds")
dataset = dataset[None]
meta = meta[None]
adata = anndata.AnnData(X=dataset.T, obs=meta)
adata.X = adata.X.astype(np.float32)

#######################
# scvi model
#######################
scvi.model.SCVI.setup_anndata(adata, batch_key=condition_key)
vae = scvi.model.SCVI(adata, n_layers=3, n_latent=npcs, gene_likelihood="nb")
vae.train()#set batch size for this dataset batch_size = 100
res_scvi = pd.DataFrame(vae.get_latent_representation().T)
res_scvi.to_csv("/home/server/zy/group_scripts/data_results/" + data_dir + "/res_scvi_raw.csv")

#######################
# scanvi model
#######################
adata = anndata.AnnData(X=dataset.T, obs=meta)
adata.X = adata.X.astype(np.float32)
scvi.model.SCVI.setup_anndata(adata, batch_key=condition_key)
vae = scvi.model.SCVI(adata, n_layers=3, n_latent=npcs, gene_likelihood="nb")
lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    adata=adata,
    labels_key=cell_type_key,
    unlabeled_category="Unknown",
)
lvae.train(max_epochs=20, n_samples_per_label=100)#set batch size for this dataset , batch_size = 100
res_scanvi = pd.DataFrame(lvae.get_latent_representation(adata).T)
res_scanvi.to_csv("/home/server/zy/group_scripts/data_results/" + data_dir + "/res_scanvi_raw.csv")

#######################
# scpoli model
#######################
scpoli_model = scPoli(
    adata=adata,
    condition_keys=condition_key,
    cell_type_keys=cell_type_key,
    latent_dim=npcs,
    embedding_dims=5,
    recon_loss='nb',
)
scpoli_model.train(
    n_epochs=50,
    pretraining_epochs=40,
    early_stopping_kwargs=early_stopping_kwargs,
    eta=5,
)
res_scpoli = pd.DataFrame(scpoli_model.get_latent(adata).T)
res_scpoli.to_csv("/home/server/zy/group_scripts/data_results/" + data_dir + "/res_scpoli_raw.csv")

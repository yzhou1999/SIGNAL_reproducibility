import scanpy as sc
import scvi
import numpy as np
import pandas as pd
import pyreadr
# from scipy.sparse import csr_matrix
import anndata
from sklearn.metrics import classification_report
from scarches.models.scpoli import scPoli
import warnings
warnings.filterwarnings('ignore')

data_dir = "Human_MTG"
condition_key = 'Batch'
cell_type_key = 'SubClass'
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

# prepare ref and query adata object
ref_dataset = pyreadr.read_r("/home/server/zy/group_scripts/datasets_preparation/" + data_dir + "/count_ref_dataset.rds")
ref_meta = pyreadr.read_r("/home/server/zy/group_scripts/datasets_preparation/" + data_dir + "/ref_meta.rds")
ref_dataset = ref_dataset[None]
ref_meta = ref_meta[None]
ref_adata = anndata.AnnData(X=ref_dataset.T, obs=ref_meta)
ref_adata.X = ref_adata.X.astype(np.float32)
query_dataset = pyreadr.read_r("/home/server/zy/group_scripts/datasets_preparation/" + data_dir + "/count_query_dataset.rds")
query_meta = pyreadr.read_r("/home/server/zy/group_scripts/datasets_preparation/" + data_dir + "/query_meta.rds")
query_dataset = query_dataset[None]
query_meta = query_meta[None]
query_adata = anndata.AnnData(X=query_dataset.T, obs=query_meta)
query_adata.X = query_adata.X.astype(np.float32)

# scvi model
scvi.model.SCVI.setup_anndata(ref_adata, batch_key=condition_key)
vae = scvi.model.SCVI(ref_adata, n_layers=3, n_latent=npcs, gene_likelihood="nb")
vae.train()
ref_adata.obs["labels_scanvi"] = ref_adata.obs[cell_type_key].values
# scanvi model
lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    adata=ref_adata,
    labels_key="labels_scanvi",
    unlabeled_category="Unknown",
)
lvae.train(max_epochs=20, n_samples_per_label=100)
dir_path_scan = "/model_scanvi/"
lvae.save("/home/server/zy/group_scripts/datasets_preparation/" + data_dir + dir_path_scan, overwrite=True)
# scvi.model.SCANVI.prepare_query_anndata(query_adata, "/home/server/zy/group_scripts/datasets_preparation/" + data_dir + dir_path_scan)
# reference mapping
vae_q = scvi.model.SCANVI.load_query_data(
    query_adata,
    "/home/server/zy/group_scripts/datasets_preparation/" + data_dir + dir_path_scan,
)
vae_q.train(
    max_epochs=100,
    plan_kwargs=dict(weight_decay=0.0),
    check_val_every_n_epoch=10,
)
# save mapped result
res_scanvi = pd.DataFrame(np.concatenate((lvae.get_latent_representation(ref_adata).T, vae_q.get_latent_representation(query_adata).T), axis = 1))
res_scanvi.to_csv("/home/server/zy/group_scripts/data_results/" + data_dir + "/mapped_scanvi.csv")
# save transfered labels
pred_scanvi = pd.DataFrame(vae_q.predict())
pred_scanvi.to_csv("/home/server/zy/group_scripts/data_results/" + data_dir + "/pred_labels_scanvi.csv")

# scpoli model
scpoli_model = scPoli(
    adata=ref_adata,
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
# reference mapping
scpoli_query = scPoli.load_query_data(
    adata=query_adata,
    reference_model=scpoli_model,
    labeled_indices=[],
)
scpoli_query.train(
    n_epochs=50,
    pretraining_epochs=40,
    eta=10
)
results_dict = scpoli_query.classify(query_adata, scale_uncertainties=True)
# save mapped result
res_scpoli = pd.DataFrame(np.concatenate((scpoli_query.get_latent(ref_adata, mean = True).T, 
                                          scpoli_query.get_latent(query_adata, mean = True).T), axis = 1))
res_scpoli.to_csv("/home/server/zy/group_scripts/data_results/" + data_dir + "/mapped_scpoli.csv")
# save transfered labels
pred_scpoli = pd.DataFrame(results_dict[cell_type_key]["preds"])
pred_scpoli.to_csv("/home/server/zy/group_scripts/data_results/" + data_dir + "/pred_labels_scpoli.csv")
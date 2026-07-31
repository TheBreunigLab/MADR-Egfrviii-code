#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  1 16:43:22 2026

@author: JJB
"""

# %%
import scanpy as sc
from pathlib import Path
import re

# Base directory
# replace "" contents and '*' with path to matrix files
base = Path("/mnt/****/***/*[parse plate]")

# Sample folders (excluding "all-sample")
samples = [
    "EvIII-PC_3wk_Female_D20",
    "EvIII-PC_3wk_Male_D20",
    "EvIII-PC_5wk_Female_D34",
    "EvIII-PC_5wk_Male_D34",
    "EvIII-PC_7_wk_135057",
    "EvIII-PC_7_wk_Male_135059",
]

# %%



def parse_metadata(sample_name):
    week_match = re.search(r"(\d+)\s*_?\s*wk", sample_name)
    week = f"{week_match.group(1)}wk" if week_match else "NA"

    if sample_name == "EvIII-PC_7_wk_135057":
        sex = "Female"
    elif "Female" in sample_name:
        sex = "Female"
    elif "Male" in sample_name:
        sex = "Male"
    else:
        sex = "NA"

    return week, sex

adatas = []

for s in samples:
    d = base / s / "DGE_filtered"

    print(f"Loading {s}")

    ad = sc.read_mtx(d / "count_matrix.mtx")

    genes = pd.read_csv(d / "all_genes.csv")
    cells = pd.read_csv(d / "cell_metadata.csv")

    gene_col = "gene_name" if "gene_name" in genes.columns else genes.columns[0]
    cell_col = "barcode" if "barcode" in cells.columns else cells.columns[0]

    ad.var_names = genes[gene_col].astype(str).values
    ad.obs_names = cells[cell_col].astype(str).values
    ad.var_names_make_unique()

    # Add cell metadata first, avoiding duplicated columns
    cells2 = cells.set_index(cell_col)
    cells2 = cells2.loc[:, ~cells2.columns.isin(ad.obs.columns)]
    ad.obs = ad.obs.join(cells2, how="left")

    week, sex = parse_metadata(s)

    # Overwrite/add clean metadata columns
    ad.obs["sample"] = s
    ad.obs["Week"] = week
    ad.obs["Sex"] = sex

    adatas.append(ad)

adata = sc.concat(
    adatas,
    label="batch",
    keys=samples,
    index_unique="-"
)

print(adata)
print(adata.obs[["sample", "Week", "Sex"]].drop_duplicates())

# %%
# Quick sanity check
print(adata)
print(adata.obs[["sample", "Week", "Sex"]].drop_duplicates())


scrubbed = []

for s in adata.obs["sample"].unique():
    print(f"Running Scrublet on {s}")

    ad_s = adata[adata.obs["sample"] == s].copy()

    sc.pp.scrublet(
        ad_s,
        batch_key=None,
        expected_doublet_rate=0.06
    )

    scrubbed.append(ad_s.obs[["doublet_score", "predicted_doublet"]])

# combine results
scrublet_obs = pd.concat(scrubbed)

# add back to full adata
adata.obs["doublet_score"] = scrublet_obs.loc[adata.obs_names, "doublet_score"]
adata.obs["predicted_doublet"] = scrublet_obs.loc[adata.obs_names, "predicted_doublet"]

# check by sample
pd.crosstab(adata.obs["sample"], adata.obs["predicted_doublet"])

# %%
# %%
sc.pl.violin(
    adata,
    "doublet_score",
    groupby="sample",
    rotation=90
)

sc.pl.umap(
    adata,
    color=["sample", "doublet_score", "predicted_doublet"]
)


# Mito genes (mouse)
adata.var["mt"] = adata.var_names.str.startswith(("mt-", "MT-"))

sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    percent_top=None,
    log1p=False,
    inplace=True
)

# Quick look
sc.pl.violin(
    adata,
    ["total_counts", "n_genes_by_counts", "pct_counts_mt"],
    groupby="sample",
    rotation=90
)




# %%
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(10, 4))

for ax, sex in zip(axes, ['Male', 'Female']):
    sc.pl.umap(
        adata[adata.obs['Sex'] == sex],
        color='Week',
        ax=ax,
        show=False,
        title=sex
    )

plt.tight_layout()
plt.show()

#%%
fig, axes = plt.subplots(1, 3, figsize=(15, 4))

for ax, week in zip(axes, ['3wk', '5wk', '7wk']):
    sc.pl.umap(
        adata[adata.obs['Week'] == week],
        color='Sex',
        ax=ax,
        show=False,
        title=week
    )

plt.tight_layout()
plt.show()

#%%
# %%
import matplotlib.pyplot as plt
import scanpy as sc

weeks = ["3wk", "5wk", "7wk"]

fig, axes = plt.subplots(
    1,
    len(weeks),
    figsize=(4 * len(weeks), 4),
    constrained_layout=True
)

for ax, week in zip(axes, weeks):
    sc.pl.umap(
        adata[adata.obs["Week"] == week],
        color="Sex",
        ax=ax,
        show=False,
        title=week,
        frameon=False,
        legend_loc="right margin"
    )

plt.show()

# %%
sexes = ["Male", "Female"]

fig, axes = plt.subplots(
    1,
    len(sexes),
    figsize=(4 * len(sexes), 4),
    constrained_layout=True
)

for ax, sex in zip(axes, sexes):
    sc.pl.umap(
        adata[adata.obs["Sex"] == sex],
        color="Week",
        ax=ax,
        show=False,
        title=sex,
        frameon=False,
        legend_loc="right margin"
    )

plt.show()

# %%
samples = sorted(adata.obs["sample"].unique())

fig, axes = plt.subplots(
    2,
    3,
    figsize=(12, 8),
    constrained_layout=True
)

axes = axes.flatten()

for ax, sample in zip(axes, samples):
    sc.pl.umap(
        adata[adata.obs["sample"] == sample],
        color="Week",
        ax=ax,
        show=False,
        title=sample,
        frameon=False,
        legend_loc=None
    )

plt.show()

# %%
import numpy as np
import matplotlib.pyplot as plt

samples = sorted(adata.obs["sample"].unique())

n_cols = 2
n_rows = int(np.ceil(len(samples) / n_cols))

fig, axes = plt.subplots(
    n_rows,
    n_cols,
    figsize=(6 * n_cols, 4 * n_rows),
    constrained_layout=True
)

axes = axes.flatten()

for ax, sample in zip(axes, samples):
    sc.pl.umap(
        adata[adata.obs["sample"] == sample],
        color="Week",
        ax=ax,
        show=False,
        title=sample,
        frameon=False,
        legend_loc=None
    )

# Turn off any unused panels
for ax in axes[len(samples):]:
    ax.axis("off")

plt.show()


# %%
# Rename Leiden clusters into cell-type labels

cluster_to_celltype = {
    "0": "DAM MG",
    "4": "DAM MG",

    "1": "Endothelial",
    "2": "Endothelial",

    "3": "Tumor OPCs",

    "5": "Ependymal",

    "6": "Pericytes",
    "9": "Pericytes",
    "13": "Pericytes",

    "7": "Tumor cycling",

    "8": "Mature oligos",

    "10": "Tumor astro-like",

    "11": "Trans. Glial Progenitor",

    "12": "Neurons",

    "14": "Hom. MG",

    "15": "T cells",
}

# make sure leiden is string
adata.obs["leiden"] = adata.obs["leiden"].astype(str)

# create new annotation column
adata.obs["celltype"] = adata.obs["leiden"].map(cluster_to_celltype)

# check for any unmapped clusters
print(adata.obs[["leiden", "celltype"]].drop_duplicates().sort_values("leiden"))

# plot
sc.pl.umap(
    adata,
    color="celltype",
    legend_loc="right margin",
    legend_fontsize=8,
    frameon=False
)
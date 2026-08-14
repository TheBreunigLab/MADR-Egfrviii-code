import numpy as np
from anndata import AnnData
import pandas as pd
import scanpy as sc
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=600, dpi_save=1200, color_map='OrRd')
sc.settings.vector_friendly = False
sc.settings.n_jobs = 56
import scanpy.external as sce
import seaborn as sns
import loompy
import os, sys
import sklearn
import scipy
import matplotlib.colors
import matplotlib.pyplot as plt
plt.rc('font', size=14)
import plotly.io as pio
pio.renderers.default='browser' #This set to plot the graph on the browser
#pio.renderers.default='svg' #This set to plot the graph in Spyder
import plotly.express as px
import plotly.graph_objects as go
os.chdir('/media/david/4TBNvMe/scRNAseq/Analysis/Adult Glioma')


Human_final= sc.read("/media/david/4TBNvMe/scRNAseq/Analyzed datasets/2026_05_06 - Human EgfrvIII final object.h5ad")
Mouse_final= sc.read('/media/david/4TBNvMe/scRNAseq/Analyzed datasets/2026_05_06 - Final EgfrvIII mouse tumor 10X and parse integration datasets.h5ad')


#Create a function to generate a dataframe with the summarized information
def get_cluster_proportions(adata,
                            cluster_key="cluster_final",
                            sample_key="replicate",
                            drop_values=None):
    """
    Input
    =====
    adata : AnnData object
    cluster_key : key of `adata.obs` storing cluster info
    sample_key : key of `adata.obs` storing sample/replicate info
    drop_values : list/iterable of possible values of `sample_key` that you don't want
    
    Returns
    =======
    pd.DataFrame with samples as the index and clusters as the columns and 0-100 floats
    as values
    """
    
    adata_tmp = adata.copy()
    sizes = adata_tmp.obs.groupby([cluster_key, sample_key]).size()
    sizes = sizes.dropna()
    props = sizes.groupby(level=sizes.index.names[1], observed=False).apply(lambda x: 100 * x / x.sum())
    props.index = props.index.droplevel(0)
    props_df = props.reset_index(name='proportion')
    props_df = props_df.pivot(columns=sample_key, index=cluster_key, values='proportion').T
    props_df.fillna(0, inplace=True)
    
    if drop_values is not None:
        for drop_value in drop_values:
            props_df.drop(drop_value, axis=0, inplace=True)
    return props_df

datasets = Mouse_final.concatenate(Human_final,join='outer')
Mouse_proportions = get_cluster_proportions(datasets, cluster_key='Clustering', sample_key='species')

#Merge in a numpy array all columns to simplify the generation of the figure
Mouse_proportions = Mouse_proportions.iloc[:,[11,12,13,14,15,5,6,1,4,10,0,7,8,2,3,9]]
Mouse_proportions = Mouse_proportions.iloc[[1,0],:]

values = np.array([Mouse_proportions.iloc[:,0], Mouse_proportions.iloc[:,1], Mouse_proportions.iloc[:,2],
                   Mouse_proportions.iloc[:,3], Mouse_proportions.iloc[:,4], Mouse_proportions.iloc[:,5],
                   Mouse_proportions.iloc[:,6], Mouse_proportions.iloc[:,7], Mouse_proportions.iloc[:,8],
                   Mouse_proportions.iloc[:,9], Mouse_proportions.iloc[:,10], Mouse_proportions.iloc[:,11],
                   Mouse_proportions.iloc[:,12], Mouse_proportions.iloc[:,13], Mouse_proportions.iloc[:,14],
                   Mouse_proportions.iloc[:,15]])

plt.rcParams["figure.figsize"] = [2.5, 7] #Determine the proportions of the figure

         
colors = ["#FF7043",  # Tumor Astro-like (orange-red)
          "#FF5252",  # Tumor Cycling (bright red-pink)
          "#e54918",  # Tumor Mesenchymal-like
          "#C62828",  # Tumor Neuron-like (deep crimson)
          "#FFCDD2",  # Tumor OPC-like (light coral)]
          "#984EA3",  # Microglia (deep purple)
          "#b29cdb",  # Microglia Cycling (light green)
          "#F781BF",  # Dendritic Cells (pink)
          "#CE93D8",  # Peripheral Immune (orchid)
          "#E1BEE7",  # Peripheral Macrophages (pale lavender)
          "#4F81BD",  # Astrocytes (blue)
          "#FFD54F",  # Neurons (yellow)
          "#A6CEE3",  # Oligodendrocytes (light blue)
          "#43A047",  # Endothelial Cells (medium green)
          "#81C784",  # Ependymal Cells (soft pastel green)
          "#66BB6A"] # Pericytes (leaf green)



# Run all remaining line together from here
fig, ax = plt.subplots(dpi=600)
# Stacked bar chart with loop
for i in range(Mouse_proportions.shape[1]):
  ax.bar(Mouse_proportions.index, Mouse_proportions.iloc[:,i], 
         width=0.8,edgecolor = "black", linewidth = 0.25,
         bottom = np.sum(values[:i], axis = 0), 
         label=Mouse_proportions.columns.values.tolist()[i],
         color=colors[i])
ax.grid(False)
ax.set_xlabel(Mouse_proportions.index.name.capitalize())
ax.set_ylabel("Proportion")
sns.despine(fig, ax)
ax.legend(bbox_to_anchor=(1.01, 0.97), frameon=False, title="Cell Type")
plt.show()

# Only Immune populations

immune_cols = ['Microglia',
               'Microglia Cycling',
               'Dendritic Cells',
               'Lymphocytes',
               'Peripheral Macrophages']

Mouse_immune = Mouse_proportions[immune_cols].copy()
Mouse_immune = Mouse_immune.div(Mouse_immune.sum(axis=1), axis=0) * 100
# Prepare values for stacking
values_immune = np.array([Mouse_immune.iloc[:,i] for i in range(Mouse_immune.shape[1])])


immune_colors = ['#984EA3', #'Microglia',
                 '#b29cdb', #'Microglia Cycling',
                 '#F781BF', #'Dendritic Cells',
                 '#CE93D8', #'Lymphocytes',
                 '#E1BEE7'] #'Peripheral Macrophages',

fig, ax = plt.subplots(dpi=600)

for i in range(Mouse_immune.shape[1]):
    ax.bar(Mouse_immune.index,
           Mouse_immune.iloc[:,i],
           width=0.8,
           edgecolor="black",
           linewidth=0.25,
           bottom=np.sum(values_immune[:i], axis=0),
           label=Mouse_immune.columns[i],
           color=immune_colors[i])

ax.grid(False)
ax.set_xlabel(Mouse_immune.index.name.capitalize())
ax.set_ylabel("Proportion (immune only)")
sns.despine(fig, ax)
ax.legend(bbox_to_anchor=(1.01, 0.97), frameon=False, title="Immune Cell Type")
plt.show()
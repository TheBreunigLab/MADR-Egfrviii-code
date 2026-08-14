#Human EGFr analysis for Katie

import numpy as np
import pandas as pd
import scanpy as sc
sc.settings.verbosity = 3    # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=600, dpi_save=600, color_map='OrRd')
sc.settings.n_jobs = 56
import scanpy.external as sce
import scvelo as scv
scv.logging.print_version()
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.settings.set_figure_params('scvelo')  # for beautified visualization
import seaborn as sns
import loompy
import os, sys
import graphtools as gt
import phate
import scprep
import sklearn
import scipy
import anndata
import matplotlib.pyplot as plt
plt.rc('font', size=14)
import plotly.io as pio
pio.renderers.default='browser' #This set to plot the graph on the browser
#pio.renderers.default='svg' #This set to plot the graph in Spyder
import plotly.express as px
#Velocity cones 3D graph
import plotly.graph_objects as go


os.chdir('/mnt/4TBNvMe/scRNAseq/Analysis/Adult Glioma')  ###E18 10X

IDH_1 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/01.-IDH/GSM5319533_SF9715/GSM5319533_SF9715_matrix.mtx.gz')
IDH_1 = IDH_1.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/01.-IDH/GSM5319533_SF9715/GSM5319533_SF9715_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/01.-IDH/GSM5319533_SF9715/GSM5319533_SF9715_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
IDH_1.obs.index = obs1['cells'].astype(str)
IDH_1.var.index = var1['var_names'].astype(str)

IDH_2 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/01.-IDH/GSM6199405_SF9715v2/GSM6199405_SF9715v2_matrix.mtx.gz')
IDH_2 = IDH_2.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/01.-IDH/GSM6199405_SF9715v2/GSM6199405_SF9715v2_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/01.-IDH/GSM6199405_SF9715v2/GSM6199405_SF9715v2_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
IDH_2.obs.index = obs1['cells'].astype(str)
IDH_2.var.index = var1['var_names'].astype(str)

IDH_EGFRsnv = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/02.-IDH & EGFR-SNV/GSM5319526_SF5581/GSM5319526_SF5581_matrix.mtx.gz')
IDH_EGFRsnv = IDH_EGFRsnv.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/02.-IDH & EGFR-SNV/GSM5319526_SF5581/GSM5319526_SF5581_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/02.-IDH & EGFR-SNV/GSM5319526_SF5581/GSM5319526_SF5581_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
IDH_EGFRsnv.obs.index = obs1['cells'].astype(str)
IDH_EGFRsnv.var.index = var1['var_names'].astype(str)

IDH_CDKN2A_1 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/03.-IDH & CDKN2A/GSM5319529_SF8963/GSM5319529_SF8963_matrix.mtx')
IDH_CDKN2A_1 = IDH_CDKN2A_1.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/03.-IDH & CDKN2A/GSM5319529_SF8963/GSM5319529_SF8963_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/03.-IDH & CDKN2A/GSM5319529_SF8963/GSM5319529_SF8963_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
IDH_CDKN2A_1.obs.index = obs1['cells'].astype(str)
IDH_CDKN2A_1.var.index = var1['var_names'].astype(str)

IDH_CDKN2A_2 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/03.-IDH & CDKN2A/GSM5319554_SF6621/GSM5319554_SF6621_matrix.mtx.gz')
IDH_CDKN2A_2 = IDH_CDKN2A_2.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/03.-IDH & CDKN2A/GSM5319554_SF6621/GSM5319554_SF6621_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/03.-IDH & CDKN2A/GSM5319554_SF6621/GSM5319554_SF6621_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
IDH_CDKN2A_2.obs.index = obs1['cells'].astype(str)
IDH_CDKN2A_2.var.index = var1['var_names'].astype(str)

IDH_CDKN2A_3 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/03.-IDH & CDKN2A/GSM5319567_SF12165/GSM5319567_SF12165_matrix.mtx.gz')
IDH_CDKN2A_3 = IDH_CDKN2A_3.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/03.-IDH & CDKN2A/GSM5319567_SF12165/GSM5319567_SF12165_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/03.-IDH & CDKN2A/GSM5319567_SF12165/GSM5319567_SF12165_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
IDH_CDKN2A_3.obs.index = obs1['cells'].astype(str)
IDH_CDKN2A_3.var.index = var1['var_names'].astype(str)

IDH_CDKN2A_3b = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/03.-IDH & CDKN2A/GSM5319567_SF12165_batch2/GSM5319567_SF12165_batch2_matrix.mtx.gz')
IDH_CDKN2A_3b = IDH_CDKN2A_3b.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/03.-IDH & CDKN2A/GSM5319567_SF12165_batch2/GSM5319567_SF12165_batch2_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/03.-IDH & CDKN2A/GSM5319567_SF12165_batch2/GSM5319567_SF12165_batch2_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
IDH_CDKN2A_3b.obs.index = obs1['cells'].astype(str)
IDH_CDKN2A_3b.var.index = var1['var_names'].astype(str)

NF1_1 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/04.-NF1/GSM5319517_SF2501/GSM5319517_SF2501_matrix.mtx.gz')
NF1_1 = NF1_1.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/04.-NF1/GSM5319517_SF2501/GSM5319517_SF2501_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/04.-NF1/GSM5319517_SF2501/GSM5319517_SF2501_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
NF1_1.obs.index = obs1['cells'].astype(str)
NF1_1.var.index = var1['var_names'].astype(str)

NF1_2 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/04.-NF1/GSM5319520_SF3076/GSM5319520_SF3076_matrix.mtx.gz')
NF1_2 = NF1_2.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/04.-NF1/GSM5319520_SF3076/GSM5319520_SF3076_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/04.-NF1/GSM5319520_SF3076/GSM5319520_SF3076_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
NF1_2.obs.index = obs1['cells'].astype(str)
NF1_2.var.index = var1['var_names'].astype(str)

NF1_3 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/04.-NF1/GSM5319548_SF2979/GSM5319548_SF2979_matrix.mtx.gz')
NF1_3 = NF1_3.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/04.-NF1/GSM5319548_SF2979/GSM5319548_SF2979_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/04.-NF1/GSM5319548_SF2979/GSM5319548_SF2979_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
NF1_3.obs.index = obs1['cells'].astype(str)
NF1_3.var.index = var1['var_names'].astype(str)

NF1_4 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/04.-NF1/GSM5319552_SF6098/GSM5319552_SF6098_matrix.mtx.gz')
NF1_4 = NF1_4.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/04.-NF1/GSM5319552_SF6098/GSM5319552_SF6098_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/04.-NF1/GSM5319552_SF6098/GSM5319552_SF6098_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
NF1_4.obs.index = obs1['cells'].astype(str)
NF1_4.var.index = var1['var_names'].astype(str)

NF1_5 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/04.-NF1/GSM6199404_SF11873/GSM6199404_SF11873_matrix.mtx.gz')
NF1_5 = NF1_5.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/04.-NF1/GSM6199404_SF11873/GSM6199404_SF11873_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/04.-NF1/GSM6199404_SF11873/GSM6199404_SF11873_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
NF1_5.obs.index = obs1['cells'].astype(str)
NF1_5.var.index = var1['var_names'].astype(str)

NF1_EGFRsnv = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/05.-NF1 & EGFR-SNV/GSM5319538_SF10514/GSM5319538_SF10514_matrix.mtx.gz')
NF1_EGFRsnv = NF1_EGFRsnv.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/05.-NF1 & EGFR-SNV/GSM5319538_SF10514/GSM5319538_SF10514_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/05.-NF1 & EGFR-SNV/GSM5319538_SF10514/GSM5319538_SF10514_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
NF1_EGFRsnv.obs.index = obs1['cells'].astype(str)
NF1_EGFRsnv.var.index = var1['var_names'].astype(str)

NF1_PTEN_1 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/06.-NF1 & PTEN/GSM5319514_SF12115/GSM5319514_SF12115_matrix.mtx.gz')
NF1_PTEN_1 = NF1_PTEN_1.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/06.-NF1 & PTEN/GSM5319514_SF12115/GSM5319514_SF12115_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/06.-NF1 & PTEN/GSM5319514_SF12115/GSM5319514_SF12115_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
NF1_PTEN_1.obs.index = obs1['cells'].astype(str)
NF1_PTEN_1.var.index = var1['var_names'].astype(str)

NF1_PTEN_2 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/06.-NF1 & PTEN/GSM5319569_SF12751/GSM5319569_SF12751_matrix.mtx.gz')
NF1_PTEN_2 = NF1_PTEN_2.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/06.-NF1 & PTEN/GSM5319569_SF12751/GSM5319569_SF12751_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/06.-NF1 & PTEN/GSM5319569_SF12751/GSM5319569_SF12751_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
NF1_PTEN_2.obs.index = obs1['cells'].astype(str)
NF1_PTEN_2.var.index = var1['var_names'].astype(str)

NF1_CDKN2A = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/07.-NF1 & CDKN2A/GSM5319556_SF7062/GSM5319556_SF7062_matrix.mtx.gz')
NF1_CDKN2A = NF1_CDKN2A.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/07.-NF1 & CDKN2A/GSM5319556_SF7062/GSM5319556_SF7062_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/07.-NF1 & CDKN2A/GSM5319556_SF7062/GSM5319556_SF7062_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
NF1_CDKN2A.obs.index = obs1['cells'].astype(str)
NF1_CDKN2A.var.index = var1['var_names'].astype(str)

NF1_PTEN_CDKN2A_1 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/09.-NF1, PTEN & CDKN2A/GSM5319524_SF4810/GSM5319524_SF4810_matrix.mtx.gz')
NF1_PTEN_CDKN2A_1 = NF1_PTEN_CDKN2A_1.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/09.-NF1, PTEN & CDKN2A/GSM5319524_SF4810/GSM5319524_SF4810_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/09.-NF1, PTEN & CDKN2A/GSM5319524_SF4810/GSM5319524_SF4810_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
NF1_PTEN_CDKN2A_1.obs.index = obs1['cells'].astype(str)
NF1_PTEN_CDKN2A_1.var.index = var1['var_names'].astype(str)

NF1_PTEN_CDKN2A_1b = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/09.-NF1, PTEN & CDKN2A/GSM5319524_SF4810_batch2/GSM5319524_SF4810_batch2_matrix.mtx.gz')
NF1_PTEN_CDKN2A_1b = NF1_PTEN_CDKN2A_1b.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/09.-NF1, PTEN & CDKN2A/GSM5319524_SF4810_batch2/GSM5319524_SF4810_batch2_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/09.-NF1, PTEN & CDKN2A/GSM5319524_SF4810_batch2/GSM5319524_SF4810_batch2_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
NF1_PTEN_CDKN2A_1b.obs.index = obs1['cells'].astype(str)
NF1_PTEN_CDKN2A_1b.var.index = var1['var_names'].astype(str)

NF1_PTEN_CDKN2A_2 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/09.-NF1, PTEN & CDKN2A/GSM5319534_SF10108/GSM5319534_SF10108_matrix.mtx.gz')
NF1_PTEN_CDKN2A_2 = NF1_PTEN_CDKN2A_2.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/09.-NF1, PTEN & CDKN2A/GSM5319534_SF10108/GSM5319534_SF10108_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/09.-NF1, PTEN & CDKN2A/GSM5319534_SF10108/GSM5319534_SF10108_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
NF1_PTEN_CDKN2A_2.obs.index = obs1['cells'].astype(str)
NF1_PTEN_CDKN2A_2.var.index = var1['var_names'].astype(str)

NF1_PTEN_CDKN2A_3 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/09.-NF1, PTEN & CDKN2A/GSM5319575_SF9871/GSM5319575_SF9871_matrix.mtx.gz')
NF1_PTEN_CDKN2A_3 = NF1_PTEN_CDKN2A_3.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/09.-NF1, PTEN & CDKN2A/GSM5319575_SF9871/GSM5319575_SF9871_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/09.-NF1, PTEN & CDKN2A/GSM5319575_SF9871/GSM5319575_SF9871_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
NF1_PTEN_CDKN2A_3.obs.index = obs1['cells'].astype(str)
NF1_PTEN_CDKN2A_3.var.index = var1['var_names'].astype(str)

NF1_EGFRamp_PTEN_CDKN2A = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/10.-NF1, EGFR-amplification, PTEN & CDKN2A/GSM5319507_SF11344/GSM5319507_SF11344_matrix.mtx.gz')
NF1_EGFRamp_PTEN_CDKN2A = NF1_EGFRamp_PTEN_CDKN2A.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/10.-NF1, EGFR-amplification, PTEN & CDKN2A/GSM5319507_SF11344/GSM5319507_SF11344_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/10.-NF1, EGFR-amplification, PTEN & CDKN2A/GSM5319507_SF11344/GSM5319507_SF11344_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
NF1_EGFRamp_PTEN_CDKN2A.obs.index = obs1['cells'].astype(str)
NF1_EGFRamp_PTEN_CDKN2A.var.index = var1['var_names'].astype(str)

EGFRsnv_1 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/11.-EGFR-SNV/GSM5319512_SF11977/GSM5319512_SF11977_matrix.mtx.gz')
EGFRsnv_1 = EGFRsnv_1.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/11.-EGFR-SNV/GSM5319512_SF11977/GSM5319512_SF11977_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/11.-EGFR-SNV/GSM5319512_SF11977/GSM5319512_SF11977_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
EGFRsnv_1.obs.index = obs1['cells'].astype(str)
EGFRsnv_1.var.index = var1['var_names'].astype(str)

EGFRsnv_2 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/11.-EGFR-SNV/GSM5319519_SF2990/GSM5319519_SF2990_matrix.mtx.gz')
EGFRsnv_2 = EGFRsnv_2.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/11.-EGFR-SNV/GSM5319519_SF2990/GSM5319519_SF2990_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/11.-EGFR-SNV/GSM5319519_SF2990/GSM5319519_SF2990_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
EGFRsnv_2.obs.index = obs1['cells'].astype(str)
EGFRsnv_2.var.index = var1['var_names'].astype(str)

EGFRsnv_3 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/11.-EGFR-SNV/GSM5319532_SF9494/GSM5319532_SF9494_matrix.mtx.gz')
EGFRsnv_3 = EGFRsnv_3.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/11.-EGFR-SNV/GSM5319532_SF9494/GSM5319532_SF9494_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/11.-EGFR-SNV/GSM5319532_SF9494/GSM5319532_SF9494_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
EGFRsnv_3.obs.index = obs1['cells'].astype(str)
EGFRsnv_3.var.index = var1['var_names'].astype(str)

EGFRsnv_4 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/11.-EGFR-SNV/GSM5319549_SF3073/GSM5319549_SF3073_matrix.mtx.gz')
EGFRsnv_4 = EGFRsnv_4.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/11.-EGFR-SNV/GSM5319549_SF3073/GSM5319549_SF3073_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/11.-EGFR-SNV/GSM5319549_SF3073/GSM5319549_SF3073_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
EGFRsnv_4.obs.index = obs1['cells'].astype(str)
EGFRsnv_4.var.index = var1['var_names'].astype(str)

EGFRsnv_CDKN2A = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/12.-EGFR-SNV & CDKN2A/GSM5319527_SF6809/GSM5319527_SF6809_matrix.mtx.gz')
EGFRsnv_CDKN2A = EGFRsnv_CDKN2A.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/12.-EGFR-SNV & CDKN2A/GSM5319527_SF6809/GSM5319527_SF6809_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/12.-EGFR-SNV & CDKN2A/GSM5319527_SF6809/GSM5319527_SF6809_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
EGFRsnv_CDKN2A.obs.index = obs1['cells'].astype(str)
EGFRsnv_CDKN2A.var.index = var1['var_names'].astype(str)

EGFRamp_CDKN2A_1 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/13.-EGFR-Amplification & CDKN2A/GSM5319508_SF11587/GSM5319508_SF11587_matrix.mtx.gz')
EGFRamp_CDKN2A_1 = EGFRamp_CDKN2A_1.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/13.-EGFR-Amplification & CDKN2A/GSM5319508_SF11587/GSM5319508_SF11587_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/13.-EGFR-Amplification & CDKN2A/GSM5319508_SF11587/GSM5319508_SF11587_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
EGFRamp_CDKN2A_1.obs.index = obs1['cells'].astype(str)
EGFRamp_CDKN2A_1.var.index = var1['var_names'].astype(str)

EGFRamp_CDKN2A_2 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/13.-EGFR-Amplification & CDKN2A/GSM5319515_SF12616/GSM5319515_SF12616_matrix.mtx.gz')
EGFRamp_CDKN2A_2 = EGFRamp_CDKN2A_2.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/13.-EGFR-Amplification & CDKN2A/GSM5319515_SF12616/GSM5319515_SF12616_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/13.-EGFR-Amplification & CDKN2A/GSM5319515_SF12616/GSM5319515_SF12616_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
EGFRamp_CDKN2A_2.obs.index = obs1['cells'].astype(str)
EGFRamp_CDKN2A_2.var.index = var1['var_names'].astype(str)

EGFRamp_CDKN2A_3 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/13.-EGFR-Amplification & CDKN2A/GSM5319563_SF12407/GSM5319563_SF12407_matrix.mtx.gz')
EGFRamp_CDKN2A_3 = EGFRamp_CDKN2A_3.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/13.-EGFR-Amplification & CDKN2A/GSM5319563_SF12407/GSM5319563_SF12407_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/13.-EGFR-Amplification & CDKN2A/GSM5319563_SF12407/GSM5319563_SF12407_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
EGFRamp_CDKN2A_3.obs.index = obs1['cells'].astype(str)
EGFRamp_CDKN2A_3.var.index = var1['var_names'].astype(str)

EGFRamp_CDKN2A_4 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/13.-EGFR-Amplification & CDKN2A/GSM5319573_SF11981/GSM5319573_SF11981_matrix.mtx.gz')
EGFRamp_CDKN2A_4 = EGFRamp_CDKN2A_4.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/13.-EGFR-Amplification & CDKN2A/GSM5319573_SF11981/GSM5319573_SF11981_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/13.-EGFR-Amplification & CDKN2A/GSM5319573_SF11981/GSM5319573_SF11981_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
EGFRamp_CDKN2A_4.obs.index = obs1['cells'].astype(str)
EGFRamp_CDKN2A_4.var.index = var1['var_names'].astype(str)

EGFRamp_CDKN2A_5 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/13.-EGFR-Amplification & CDKN2A/GSM5319574_SF12754/GSM5319574_SF12754_matrix.mtx.gz')
EGFRamp_CDKN2A_5 = EGFRamp_CDKN2A_5.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/13.-EGFR-Amplification & CDKN2A/GSM5319574_SF12754/GSM5319574_SF12754_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/13.-EGFR-Amplification & CDKN2A/GSM5319574_SF12754/GSM5319574_SF12754_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
EGFRamp_CDKN2A_5.obs.index = obs1['cells'].astype(str)
EGFRamp_CDKN2A_5.var.index = var1['var_names'].astype(str)

EGFRsnv_PTEN = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/14.-EGFR-SNV & PTEN/GSM5319550_SF3243/GSM5319550_SF3243_matrix.mtx.gz')
EGFRsnv_PTEN = EGFRsnv_PTEN.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/14.-EGFR-SNV & PTEN/GSM5319550_SF3243/GSM5319550_SF3243_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/14.-EGFR-SNV & PTEN/GSM5319550_SF3243/GSM5319550_SF3243_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
EGFRsnv_PTEN.obs.index = obs1['cells'].astype(str)
EGFRsnv_PTEN.var.index = var1['var_names'].astype(str)

EGFRamp_EGFRsnv_CDKN2A_1 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/15.-EGFR-Amplification, EGFR-SNV & CDKN2A/GSM5319509_SF11780/matrix.mtx.gz')
EGFRamp_EGFRsnv_CDKN2A_1 = EGFRamp_EGFRsnv_CDKN2A_1.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/15.-EGFR-Amplification, EGFR-SNV & CDKN2A/GSM5319509_SF11780/barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/15.-EGFR-Amplification, EGFR-SNV & CDKN2A/GSM5319509_SF11780/features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
EGFRamp_EGFRsnv_CDKN2A_1.obs.index = obs1['cells'].astype(str)
EGFRamp_EGFRsnv_CDKN2A_1.var.index = var1['var_names'].astype(str)

EGFRamp_EGFRsnv_CDKN2A_2 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/15.-EGFR-Amplification, EGFR-SNV & CDKN2A/GSM5319510_SF11815/GSM5319510_SF11815_matrix.mtx.gz')
EGFRamp_EGFRsnv_CDKN2A_2 = EGFRamp_EGFRsnv_CDKN2A_2.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/15.-EGFR-Amplification, EGFR-SNV & CDKN2A/GSM5319510_SF11815/GSM5319510_SF11815_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/15.-EGFR-Amplification, EGFR-SNV & CDKN2A/GSM5319510_SF11815/GSM5319510_SF11815_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
EGFRamp_EGFRsnv_CDKN2A_2.obs.index = obs1['cells'].astype(str)
EGFRamp_EGFRsnv_CDKN2A_2.var.index = var1['var_names'].astype(str)

EGFRamp_EGFRsnv_CDKN2A_3 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/15.-EGFR-Amplification, EGFR-SNV & CDKN2A/GSM5319511_SF11916/GSM5319511_SF11916_matrix.mtx.gz')
EGFRamp_EGFRsnv_CDKN2A_3 = EGFRamp_EGFRsnv_CDKN2A_3.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/15.-EGFR-Amplification, EGFR-SNV & CDKN2A/GSM5319511_SF11916/GSM5319511_SF11916_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/15.-EGFR-Amplification, EGFR-SNV & CDKN2A/GSM5319511_SF11916/GSM5319511_SF11916_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
EGFRamp_EGFRsnv_CDKN2A_3.obs.index = obs1['cells'].astype(str)
EGFRamp_EGFRsnv_CDKN2A_3.var.index = var1['var_names'].astype(str)

EGFRamp_EGFRsnv_CDKN2A_4 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/15.-EGFR-Amplification, EGFR-SNV & CDKN2A/GSM5319543_SF12382/GSM5319543_SF12382_matrix.mtx.gz')
EGFRamp_EGFRsnv_CDKN2A_4 = EGFRamp_EGFRsnv_CDKN2A_4.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/15.-EGFR-Amplification, EGFR-SNV & CDKN2A/GSM5319543_SF12382/GSM5319543_SF12382_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/15.-EGFR-Amplification, EGFR-SNV & CDKN2A/GSM5319543_SF12382/GSM5319543_SF12382_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
EGFRamp_EGFRsnv_CDKN2A_4.obs.index = obs1['cells'].astype(str)
EGFRamp_EGFRsnv_CDKN2A_4.var.index = var1['var_names'].astype(str)

EGFRamp_EGFRsnv_CDKN2A_5 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/15.-EGFR-Amplification, EGFR-SNV & CDKN2A/GSM5319560_SF12408/GSM5319560_SF12408_matrix.mtx.gz')
EGFRamp_EGFRsnv_CDKN2A_5 = EGFRamp_EGFRsnv_CDKN2A_5.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/15.-EGFR-Amplification, EGFR-SNV & CDKN2A/GSM5319560_SF12408/GSM5319560_SF12408_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/15.-EGFR-Amplification, EGFR-SNV & CDKN2A/GSM5319560_SF12408/GSM5319560_SF12408_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
EGFRamp_EGFRsnv_CDKN2A_5.obs.index = obs1['cells'].astype(str)
EGFRamp_EGFRsnv_CDKN2A_5.var.index = var1['var_names'].astype(str)

EGFRamp_EGFRsnv_CDKN2A_6 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/15.-EGFR-Amplification, EGFR-SNV & CDKN2A/GSM5319565_SF12243/GSM5319565_SF12243_matrix.mtx.gz')
EGFRamp_EGFRsnv_CDKN2A_6 = EGFRamp_EGFRsnv_CDKN2A_6.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/15.-EGFR-Amplification, EGFR-SNV & CDKN2A/GSM5319565_SF12243/GSM5319565_SF12243_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/15.-EGFR-Amplification, EGFR-SNV & CDKN2A/GSM5319565_SF12243/GSM5319565_SF12243_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
EGFRamp_EGFRsnv_CDKN2A_6.obs.index = obs1['cells'].astype(str)
EGFRamp_EGFRsnv_CDKN2A_6.var.index = var1['var_names'].astype(str)

EGFRamp_EGFRsnv_CDKN2A_6b = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/15.-EGFR-Amplification, EGFR-SNV & CDKN2A/GSM5319565_SF12243_batch2/GSM5319565_SF12243_batch2_matrix.mtx.gz')
EGFRamp_EGFRsnv_CDKN2A_6b = EGFRamp_EGFRsnv_CDKN2A_6b.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/15.-EGFR-Amplification, EGFR-SNV & CDKN2A/GSM5319565_SF12243_batch2/GSM5319565_SF12243_batch2_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/15.-EGFR-Amplification, EGFR-SNV & CDKN2A/GSM5319565_SF12243_batch2/GSM5319565_SF12243_batch2_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
EGFRamp_EGFRsnv_CDKN2A_6b.obs.index = obs1['cells'].astype(str)
EGFRamp_EGFRsnv_CDKN2A_6b.var.index = var1['var_names'].astype(str)

EGFRamp_PTEN_CDKN2A_1 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/16.-EGFR-Amplification, PTEN & CDKN2A/GSM5319505_SF10857/GSM5319505_SF10857_matrix.mtx.gz')
EGFRamp_PTEN_CDKN2A_1 = EGFRamp_PTEN_CDKN2A_1.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/16.-EGFR-Amplification, PTEN & CDKN2A/GSM5319505_SF10857/GSM5319505_SF10857_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/16.-EGFR-Amplification, PTEN & CDKN2A/GSM5319505_SF10857/GSM5319505_SF10857_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
EGFRamp_PTEN_CDKN2A_1.obs.index = obs1['cells'].astype(str)
EGFRamp_PTEN_CDKN2A_1.var.index = var1['var_names'].astype(str)

EGFRamp_PTEN_CDKN2A_2 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/16.-EGFR-Amplification, PTEN & CDKN2A/GSM5319539_SF11248/GSM5319539_SF11248_matrix.mtx.gz')
EGFRamp_PTEN_CDKN2A_2 = EGFRamp_PTEN_CDKN2A_2.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/16.-EGFR-Amplification, PTEN & CDKN2A/GSM5319539_SF11248/GSM5319539_SF11248_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/16.-EGFR-Amplification, PTEN & CDKN2A/GSM5319539_SF11248/GSM5319539_SF11248_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
EGFRamp_PTEN_CDKN2A_2.obs.index = obs1['cells'].astype(str)
EGFRamp_PTEN_CDKN2A_2.var.index = var1['var_names'].astype(str)

EGFRamp_PTEN_CDKN2A_3 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/16.-EGFR-Amplification, PTEN & CDKN2A/GSM5319545_SF12460/GSM5319545_SF12460_matrix.mtx.gz')
EGFRamp_PTEN_CDKN2A_3 = EGFRamp_PTEN_CDKN2A_3.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/16.-EGFR-Amplification, PTEN & CDKN2A/GSM5319545_SF12460/GSM5319545_SF12460_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/16.-EGFR-Amplification, PTEN & CDKN2A/GSM5319545_SF12460/GSM5319545_SF12460_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
EGFRamp_PTEN_CDKN2A_3.obs.index = obs1['cells'].astype(str)
EGFRamp_PTEN_CDKN2A_3.var.index = var1['var_names'].astype(str)

EGFRamp_PTEN_CDKN2A_4 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/16.-EGFR-Amplification, PTEN & CDKN2A/GSM5319564_SF10565/GSM5319564_SF10565_matrix.mtx.gz')
EGFRamp_PTEN_CDKN2A_4 = EGFRamp_PTEN_CDKN2A_4.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/16.-EGFR-Amplification, PTEN & CDKN2A/GSM5319564_SF10565/GSM5319564_SF10565_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/16.-EGFR-Amplification, PTEN & CDKN2A/GSM5319564_SF10565/GSM5319564_SF10565_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
EGFRamp_PTEN_CDKN2A_4.obs.index = obs1['cells'].astype(str)
EGFRamp_PTEN_CDKN2A_4.var.index = var1['var_names'].astype(str)

EGFRamp_PTEN_CDKN2A_4b = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/16.-EGFR-Amplification, PTEN & CDKN2A/GSM5319564_SF10565_batch2/GSM5319564_SF10565_batch2_matrix.mtx.gz')
EGFRamp_PTEN_CDKN2A_4b = EGFRamp_PTEN_CDKN2A_4b.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/16.-EGFR-Amplification, PTEN & CDKN2A/GSM5319564_SF10565_batch2/GSM5319564_SF10565_batch2_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/16.-EGFR-Amplification, PTEN & CDKN2A/GSM5319564_SF10565_batch2/GSM5319564_SF10565_batch2_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
EGFRamp_PTEN_CDKN2A_4b.obs.index = obs1['cells'].astype(str)
EGFRamp_PTEN_CDKN2A_4b.var.index = var1['var_names'].astype(str)

EGFRamp_PTEN_CDKN2A_5 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/16.-EGFR-Amplification, PTEN & CDKN2A/GSM5319570_SF12704/GSM5319570_SF12704_matrix.mtx.gz')
EGFRamp_PTEN_CDKN2A_5 = EGFRamp_PTEN_CDKN2A_5.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/16.-EGFR-Amplification, PTEN & CDKN2A/GSM5319570_SF12704/GSM5319570_SF12704_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/16.-EGFR-Amplification, PTEN & CDKN2A/GSM5319570_SF12704/GSM5319570_SF12704_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
EGFRamp_PTEN_CDKN2A_5.obs.index = obs1['cells'].astype(str)
EGFRamp_PTEN_CDKN2A_5.var.index = var1['var_names'].astype(str)

EGFRamp_PTEN_CDKN2A_5b = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/16.-EGFR-Amplification, PTEN & CDKN2A/GSM6199410_SF12704v2/GSM6199410_SF12704v2_matrix.mtx.gz')
EGFRamp_PTEN_CDKN2A_5b = EGFRamp_PTEN_CDKN2A_5b.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/16.-EGFR-Amplification, PTEN & CDKN2A/GSM6199410_SF12704v2/GSM6199410_SF12704v2_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/16.-EGFR-Amplification, PTEN & CDKN2A/GSM6199410_SF12704v2/GSM6199410_SF12704v2_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
EGFRamp_PTEN_CDKN2A_5b.obs.index = obs1['cells'].astype(str)
EGFRamp_PTEN_CDKN2A_5b.var.index = var1['var_names'].astype(str)

EGFRamp_EGFRsnv_PTEN_CDKN2A_1 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/17.-EGFR-Amplification, EGFR-SNV, PTEN & CDKN2A/GSM5319503_SF10099/GSM5319503_SF10099_matrix.mtx.gz')
EGFRamp_EGFRsnv_PTEN_CDKN2A_1 = EGFRamp_EGFRsnv_PTEN_CDKN2A_1.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/17.-EGFR-Amplification, EGFR-SNV, PTEN & CDKN2A/GSM5319503_SF10099/GSM5319503_SF10099_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/17.-EGFR-Amplification, EGFR-SNV, PTEN & CDKN2A/GSM5319503_SF10099/GSM5319503_SF10099_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
EGFRamp_EGFRsnv_PTEN_CDKN2A_1.obs.index = obs1['cells'].astype(str)
EGFRamp_EGFRsnv_PTEN_CDKN2A_1.var.index = var1['var_names'].astype(str)

EGFRamp_EGFRsnv_PTEN_CDKN2A_1b = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/17.-EGFR-Amplification, EGFR-SNV, PTEN & CDKN2A/GSM5319503_SF10099_bach2/GSM5319503_SF10099_batch2_matrix.mtx.gz')
EGFRamp_EGFRsnv_PTEN_CDKN2A_1b = EGFRamp_EGFRsnv_PTEN_CDKN2A_1b.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/17.-EGFR-Amplification, EGFR-SNV, PTEN & CDKN2A/GSM5319503_SF10099_bach2/GSM5319503_SF10099_batch2_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/17.-EGFR-Amplification, EGFR-SNV, PTEN & CDKN2A/GSM5319503_SF10099_bach2/GSM5319503_SF10099_batch2_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
EGFRamp_EGFRsnv_PTEN_CDKN2A_1b.obs.index = obs1['cells'].astype(str)
EGFRamp_EGFRsnv_PTEN_CDKN2A_1b.var.index = var1['var_names'].astype(str)

EGFRamp_EGFRsnv_PTEN_CDKN2A_2 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/17.-EGFR-Amplification, EGFR-SNV, PTEN & CDKN2A/GSM5319504_SF10432/GSM5319504_SF10432_matrix.mtx.gz')
EGFRamp_EGFRsnv_PTEN_CDKN2A_2 = EGFRamp_EGFRsnv_PTEN_CDKN2A_2.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/17.-EGFR-Amplification, EGFR-SNV, PTEN & CDKN2A/GSM5319504_SF10432/GSM5319504_SF10432_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/17.-EGFR-Amplification, EGFR-SNV, PTEN & CDKN2A/GSM5319504_SF10432/GSM5319504_SF10432_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
EGFRamp_EGFRsnv_PTEN_CDKN2A_2.obs.index = obs1['cells'].astype(str)
EGFRamp_EGFRsnv_PTEN_CDKN2A_2.var.index = var1['var_names'].astype(str)

EGFRamp_EGFRsnv_PTEN_CDKN2A_3 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/17.-EGFR-Amplification, EGFR-SNV, PTEN & CDKN2A/GSM5319542_SF12008/GSM5319542_SF12008_matrix.mtx.gz')
EGFRamp_EGFRsnv_PTEN_CDKN2A_3 = EGFRamp_EGFRsnv_PTEN_CDKN2A_3.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/17.-EGFR-Amplification, EGFR-SNV, PTEN & CDKN2A/GSM5319542_SF12008/GSM5319542_SF12008_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/17.-EGFR-Amplification, EGFR-SNV, PTEN & CDKN2A/GSM5319542_SF12008/GSM5319542_SF12008_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
EGFRamp_EGFRsnv_PTEN_CDKN2A_3.obs.index = obs1['cells'].astype(str)
EGFRamp_EGFRsnv_PTEN_CDKN2A_3.var.index = var1['var_names'].astype(str)

EGFRamp_EGFRsnv_PTEN_CDKN2A_4 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/17.-EGFR-Amplification, EGFR-SNV, PTEN & CDKN2A/GSM5319544_SF12427/GSM5319544_SF12427_matrix.mtx.gz')
EGFRamp_EGFRsnv_PTEN_CDKN2A_4 = EGFRamp_EGFRsnv_PTEN_CDKN2A_4.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/17.-EGFR-Amplification, EGFR-SNV, PTEN & CDKN2A/GSM5319544_SF12427/GSM5319544_SF12427_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/17.-EGFR-Amplification, EGFR-SNV, PTEN & CDKN2A/GSM5319544_SF12427/GSM5319544_SF12427_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
EGFRamp_EGFRsnv_PTEN_CDKN2A_4.obs.index = obs1['cells'].astype(str)
EGFRamp_EGFRsnv_PTEN_CDKN2A_4.var.index = var1['var_names'].astype(str)

EGFRamp_EGFRsnv_PTEN_CDKN2A_5 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/17.-EGFR-Amplification, EGFR-SNV, PTEN & CDKN2A/GSM5319571_SF12333/GSM5319571_SF12333_matrix.mtx.gz')
EGFRamp_EGFRsnv_PTEN_CDKN2A_5 = EGFRamp_EGFRsnv_PTEN_CDKN2A_5.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/17.-EGFR-Amplification, EGFR-SNV, PTEN & CDKN2A/GSM5319571_SF12333/GSM5319571_SF12333_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/17.-EGFR-Amplification, EGFR-SNV, PTEN & CDKN2A/GSM5319571_SF12333/GSM5319571_SF12333_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
EGFRamp_EGFRsnv_PTEN_CDKN2A_5.obs.index = obs1['cells'].astype(str)
EGFRamp_EGFRsnv_PTEN_CDKN2A_5.var.index = var1['var_names'].astype(str)

CDKN2A = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/18.-CDKN2A/GSM5319555_SF7025/GSM5319555_SF7025_matrix.mtx.gz')
CDKN2A = CDKN2A.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/18.-CDKN2A/GSM5319555_SF7025/GSM5319555_SF7025_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/18.-CDKN2A/GSM5319555_SF7025/GSM5319555_SF7025_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
CDKN2A.obs.index = obs1['cells'].astype(str)
CDKN2A.var.index = var1['var_names'].astype(str)

PTEN_1 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/19.-PTEN/GSM5319521_SF3391/GSM5319521_SF3391_matrix.mtx.gz')
PTEN_1 = PTEN_1.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/19.-PTEN/GSM5319521_SF3391/GSM5319521_SF3391_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/19.-PTEN/GSM5319521_SF3391/GSM5319521_SF3391_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
PTEN_1.obs.index = obs1['cells'].astype(str)
PTEN_1.var.index = var1['var_names'].astype(str)

PTEN_2 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/19.-PTEN/GSM5319530_SF9358/GSM5319530_SF9358_matrix.mtx.gz')
PTEN_2 = PTEN_2.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/19.-PTEN/GSM5319530_SF9358/GSM5319530_SF9358_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/19.-PTEN/GSM5319530_SF9358/GSM5319530_SF9358_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
PTEN_2.obs.index = obs1['cells'].astype(str)
PTEN_2.var.index = var1['var_names'].astype(str)

PTEN_3 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/19.-PTEN/GSM5319531_SF9372/GSM5319531_SF9372_matrix.mtx.gz')
PTEN_3 = PTEN_3.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/19.-PTEN/GSM5319531_SF9372/GSM5319531_SF9372_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/19.-PTEN/GSM5319531_SF9372/GSM5319531_SF9372_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
PTEN_3.obs.index = obs1['cells'].astype(str)
PTEN_3.var.index = var1['var_names'].astype(str)

PTEN_4 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/19.-PTEN/GSM5319541_SF10592/GSM5319541_SF10592_matrix.mtx.gz')
PTEN_4 = PTEN_4.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/19.-PTEN/GSM5319541_SF10592/GSM5319541_SF10592_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/19.-PTEN/GSM5319541_SF10592/GSM5319541_SF10592_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
PTEN_4.obs.index = obs1['cells'].astype(str)
PTEN_4.var.index = var1['var_names'].astype(str)

PTEN_5 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/19.-PTEN/GSM5319551_SF3448/GSM5319551_SF3448_matrix.mtx.gz')
PTEN_5 = PTEN_5.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/19.-PTEN/GSM5319551_SF3448/GSM5319551_SF3448_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/19.-PTEN/GSM5319551_SF3448/GSM5319551_SF3448_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
PTEN_5.obs.index = obs1['cells'].astype(str)
PTEN_5.var.index = var1['var_names'].astype(str)

PTEN_6 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/19.-PTEN/GSM5319561_SF11857/GSM5319561_SF11857_matrix.mtx.gz')
PTEN_6 = PTEN_6.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/19.-PTEN/GSM5319561_SF11857/GSM5319561_SF11857_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/19.-PTEN/GSM5319561_SF11857/GSM5319561_SF11857_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
PTEN_6.obs.index = obs1['cells'].astype(str)
PTEN_6.var.index = var1['var_names'].astype(str)

PTEN_7 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/19.-PTEN/GSM6199407_SF12774/GSM6199407_SF12774_matrix.mtx.gz')
PTEN_7 = PTEN_7.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/19.-PTEN/GSM6199407_SF12774/GSM6199407_SF12774_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/19.-PTEN/GSM6199407_SF12774/GSM6199407_SF12774_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
PTEN_7.obs.index = obs1['cells'].astype(str)
PTEN_7.var.index = var1['var_names'].astype(str)

PTEN_CDKN2A_1 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/20.-PTEN & CDKN2A/GSM5319506_SF11082/GSM5319506_SF11082_matrix.mtx.gz')
PTEN_CDKN2A_1 = PTEN_CDKN2A_1.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/20.-PTEN & CDKN2A/GSM5319506_SF11082/GSM5319506_SF11082_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/20.-PTEN & CDKN2A/GSM5319506_SF11082/GSM5319506_SF11082_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
PTEN_CDKN2A_1.obs.index = obs1['cells'].astype(str)
PTEN_CDKN2A_1.var.index = var1['var_names'].astype(str)

PTEN_CDKN2A_2 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/20.-PTEN & CDKN2A/GSM5319516_SF12707/GSM5319516_SF12707_matrix.mtx.gz')
PTEN_CDKN2A_2 = PTEN_CDKN2A_2.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/20.-PTEN & CDKN2A/GSM5319516_SF12707/GSM5319516_SF12707_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/20.-PTEN & CDKN2A/GSM5319516_SF12707/GSM5319516_SF12707_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
PTEN_CDKN2A_2.obs.index = obs1['cells'].astype(str)
PTEN_CDKN2A_2.var.index = var1['var_names'].astype(str)

PTEN_CDKN2A_3 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/20.-PTEN & CDKN2A/GSM5319523_SF4209/GSM5319523_SF4209_matrix.mtx.gz')
PTEN_CDKN2A_3 = PTEN_CDKN2A_3.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/20.-PTEN & CDKN2A/GSM5319523_SF4209/GSM5319523_SF4209_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/20.-PTEN & CDKN2A/GSM5319523_SF4209/GSM5319523_SF4209_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
PTEN_CDKN2A_3.obs.index = obs1['cells'].astype(str)
PTEN_CDKN2A_3.var.index = var1['var_names'].astype(str)

PTEN_CDKN2A_3b = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/20.-PTEN & CDKN2A/GSM5319523_SF4209_batch2/GSM5319523_SF4209_batch2_matrix.mtx.gz')
PTEN_CDKN2A_3b = PTEN_CDKN2A_3b.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/20.-PTEN & CDKN2A/GSM5319523_SF4209_batch2/GSM5319523_SF4209_batch2_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/20.-PTEN & CDKN2A/GSM5319523_SF4209_batch2/GSM5319523_SF4209_batch2_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
PTEN_CDKN2A_3b.obs.index = obs1['cells'].astype(str)
PTEN_CDKN2A_3b.var.index = var1['var_names'].astype(str)

PTEN_CDKN2A_4 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/20.-PTEN & CDKN2A/GSM5319558_SF9510/GSM5319558_SF9510_matrix.mtx.gz')
PTEN_CDKN2A_4 = PTEN_CDKN2A_4.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/20.-PTEN & CDKN2A/GSM5319558_SF9510/GSM5319558_SF9510_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/20.-PTEN & CDKN2A/GSM5319558_SF9510/GSM5319558_SF9510_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
PTEN_CDKN2A_4.obs.index = obs1['cells'].astype(str)
PTEN_CDKN2A_4.var.index = var1['var_names'].astype(str)

PTEN_CDKN2A_5 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/20.-PTEN & CDKN2A/GSM5319562_SF11488/GSM5319562_SF11488_matrix.mtx.gz')
PTEN_CDKN2A_5 = PTEN_CDKN2A_5.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/20.-PTEN & CDKN2A/GSM5319562_SF11488/GSM5319562_SF11488_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/20.-PTEN & CDKN2A/GSM5319562_SF11488/GSM5319562_SF11488_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
PTEN_CDKN2A_5.obs.index = obs1['cells'].astype(str)
PTEN_CDKN2A_5.var.index = var1['var_names'].astype(str)

PTEN_CDKN2A_6 = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/20.-PTEN & CDKN2A/GSM5319572_SF4324/GSM5319572_SF4324_matrix.mtx.gz')
PTEN_CDKN2A_6 = PTEN_CDKN2A_6.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/20.-PTEN & CDKN2A/GSM5319572_SF4324/GSM5319572_SF4324_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/20.-PTEN & CDKN2A/GSM5319572_SF4324/GSM5319572_SF4324_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
PTEN_CDKN2A_6.obs.index = obs1['cells'].astype(str)
PTEN_CDKN2A_6.var.index = var1['var_names'].astype(str)

PDGFRAamp = sc.read_mtx('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/21.-PDGFRA-Amplification & TERTp/GSM5319537_SF10484/GSM5319537_SF10484_matrix.mtx.gz')
PDGFRAamp = PDGFRAamp.transpose()
obs1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/21.-PDGFRA-Amplification & TERTp/GSM5319537_SF10484/GSM5319537_SF10484_barcodes.tsv.gz', sep="\t", header=None)
obs1.columns = ['cells']
var1 = pd.read_csv('/mnt/4TBNvMe/scRNAseq/Samples/Human/Human samples EGFR gliomas/21.-PDGFRA-Amplification & TERTp/GSM5319537_SF10484/GSM5319537_SF10484_features.tsv.gz' , sep="\t", header=None)
var1.columns = ['var_names']
PDGFRAamp.obs.index = obs1['cells'].astype(str)
PDGFRAamp.var.index = var1['var_names'].astype(str)

del[obs1, var1]
# Add metadata information in .obs

IDH_1.obs["Type of tumor"]='Recurrent'
IDH_2.obs["Type of tumor"]='Recurrent'
IDH_EGFRsnv.obs["Type of tumor"]='Primary'
IDH_CDKN2A_2.obs["Type of tumor"]='Recurrent'
IDH_CDKN2A_3.obs["Type of tumor"]='Recurrent'
IDH_CDKN2A_3b.obs["Type of tumor"]='Recurrent'
NF1_1.obs["Type of tumor"]='Primary'
NF1_2.obs["Type of tumor"]='Primary'
NF1_3.obs["Type of tumor"]='Recurrent'
NF1_4.obs["Type of tumor"]='Recurrent'
NF1_5.obs["Type of tumor"]='Primary'
NF1_EGFRsnv.obs["Type of tumor"]='Recurrent'
NF1_PTEN_1.obs["Type of tumor"]='Recurrent'
NF1_PTEN_2.obs["Type of tumor"]='Recurrent'
NF1_CDKN2A.obs["Type of tumor"]='Recurrent'
NF1_PTEN_CDKN2A_1.obs["Type of tumor"]='Primary'
NF1_PTEN_CDKN2A_1b.obs["Type of tumor"]='Primary'
NF1_PTEN_CDKN2A_2.obs["Type of tumor"]='Recurrent'
NF1_PTEN_CDKN2A_3.obs["Type of tumor"]='Recurrent'
NF1_EGFRamp_PTEN_CDKN2A.obs["Type of tumor"]='Primary'
EGFRsnv_1.obs["Type of tumor"]='Primary'
EGFRsnv_2.obs["Type of tumor"]='Primary'
EGFRsnv_3.obs["Type of tumor"]='Recurrent'
EGFRsnv_4.obs["Type of tumor"]='Recurrent'
EGFRsnv_CDKN2A.obs["Type of tumor"]='Primary'
EGFRamp_CDKN2A_1.obs["Type of tumor"]='Primary'
EGFRamp_CDKN2A_2.obs["Type of tumor"]='Recurrent'
EGFRamp_CDKN2A_3.obs["Type of tumor"]='Recurrent'
EGFRamp_CDKN2A_4.obs["Type of tumor"]='Recurrent'
EGFRamp_CDKN2A_5.obs["Type of tumor"]='Recurrent'
EGFRsnv_PTEN.obs["Type of tumor"]='Recurrent'
EGFRamp_EGFRsnv_CDKN2A_1.obs["Type of tumor"]='Primary'
EGFRamp_EGFRsnv_CDKN2A_2.obs["Type of tumor"]='Primary'
EGFRamp_EGFRsnv_CDKN2A_3.obs["Type of tumor"]='Primary'
EGFRamp_EGFRsnv_CDKN2A_4.obs["Type of tumor"]='Recurrent'
EGFRamp_EGFRsnv_CDKN2A_5.obs["Type of tumor"]='Recurrent'
EGFRamp_EGFRsnv_CDKN2A_6.obs["Type of tumor"]='Recurrent'
EGFRamp_EGFRsnv_CDKN2A_6b.obs["Type of tumor"]='Recurrent'
EGFRamp_PTEN_CDKN2A_1.obs["Type of tumor"]='Primary'
EGFRamp_PTEN_CDKN2A_2.obs["Type of tumor"]='Recurrent'
EGFRamp_PTEN_CDKN2A_3.obs["Type of tumor"]='Recurrent'
EGFRamp_PTEN_CDKN2A_4.obs["Type of tumor"]='Recurrent'
EGFRamp_PTEN_CDKN2A_4b.obs["Type of tumor"]='Recurrent'
EGFRamp_PTEN_CDKN2A_5.obs["Type of tumor"]='Recurrent'
EGFRamp_PTEN_CDKN2A_5b.obs["Type of tumor"]='Recurrent'
EGFRamp_EGFRsnv_PTEN_CDKN2A_1.obs["Type of tumor"]='Primary'
EGFRamp_EGFRsnv_PTEN_CDKN2A_1b.obs["Type of tumor"]='Primary'
EGFRamp_EGFRsnv_PTEN_CDKN2A_2.obs["Type of tumor"]='Primary'
EGFRamp_EGFRsnv_PTEN_CDKN2A_3.obs["Type of tumor"]='Recurrent'
EGFRamp_EGFRsnv_PTEN_CDKN2A_4.obs["Type of tumor"]='Recurrent'
EGFRamp_EGFRsnv_PTEN_CDKN2A_5.obs["Type of tumor"]='Recurrent'
CDKN2A.obs["Type of tumor"]='Recurrent'
PTEN_1.obs["Type of tumor"]='Primary'
PTEN_2.obs["Type of tumor"]='Primary'
PTEN_3.obs["Type of tumor"]='Primary'
PTEN_4.obs["Type of tumor"]='Primary'
PTEN_5.obs["Type of tumor"]='Recurrent'
PTEN_6.obs["Type of tumor"]='Recurrent'
PTEN_7.obs["Type of tumor"]='Recurrent'
PTEN_CDKN2A_1.obs["Type of tumor"]='Primary'
PTEN_CDKN2A_2.obs["Type of tumor"]='Primary'
PTEN_CDKN2A_3.obs["Type of tumor"]='Primary'
PTEN_CDKN2A_3b.obs["Type of tumor"]='Primary'
PTEN_CDKN2A_4.obs["Type of tumor"]='Recurrent'
PTEN_CDKN2A_5.obs["Type of tumor"]='Recurrent'
PTEN_CDKN2A_6.obs["Type of tumor"]='Recurrent'
PDGFRAamp.obs["Type of tumor"]='Recurrent'

IDH_1.obs["Sex"]='Male'
IDH_2.obs["Sex"]='Male'
IDH_EGFRsnv.obs["Sex"]='Male'
IDH_CDKN2A_2.obs["Sex"]='Male'
IDH_CDKN2A_3.obs["Sex"]='Female'
IDH_CDKN2A_3b.obs["Sex"]='Female'
NF1_1.obs["Sex"]='Female'
NF1_2.obs["Sex"]='Male'
NF1_3.obs["Sex"]='Female'
NF1_4.obs["Sex"]='Female'
NF1_5.obs["Sex"]='Female'
NF1_EGFRsnv.obs["Sex"]='Female'
NF1_PTEN_1.obs["Sex"]='Male'
NF1_PTEN_2.obs["Sex"]='Male'
NF1_CDKN2A.obs["Sex"]='Female'
NF1_PTEN_CDKN2A_1.obs["Sex"]='Female'
NF1_PTEN_CDKN2A_1b.obs["Sex"]='Female'
NF1_PTEN_CDKN2A_2.obs["Sex"]='Male'
NF1_PTEN_CDKN2A_3.obs["Sex"]='Male'
NF1_EGFRamp_PTEN_CDKN2A.obs["Sex"]='Female'
EGFRsnv_1.obs["Sex"]='Female'
EGFRsnv_2.obs["Sex"]='Female'
EGFRsnv_3.obs["Sex"]='Male'
EGFRsnv_4.obs["Sex"]='Female'
EGFRsnv_CDKN2A.obs["Sex"]='Male'
EGFRamp_CDKN2A_1.obs["Sex"]='Female'
EGFRamp_CDKN2A_2.obs["Sex"]='Male'
EGFRamp_CDKN2A_3.obs["Sex"]='Female'
EGFRamp_CDKN2A_4.obs["Sex"]='Female'
EGFRamp_CDKN2A_5.obs["Sex"]='Female'
EGFRsnv_PTEN.obs["Sex"]='Male'
EGFRamp_EGFRsnv_CDKN2A_1.obs["Sex"]='Female'
EGFRamp_EGFRsnv_CDKN2A_2.obs["Sex"]='Female'
EGFRamp_EGFRsnv_CDKN2A_3.obs["Sex"]='Male'
EGFRamp_EGFRsnv_CDKN2A_4.obs["Sex"]='Male'
EGFRamp_EGFRsnv_CDKN2A_5.obs["Sex"]='Female'
EGFRamp_EGFRsnv_CDKN2A_6.obs["Sex"]='Female'
EGFRamp_EGFRsnv_CDKN2A_6b.obs["Sex"]='Female'
EGFRamp_PTEN_CDKN2A_1.obs["Sex"]='Male'
EGFRamp_PTEN_CDKN2A_2.obs["Sex"]='Female'
EGFRamp_PTEN_CDKN2A_3.obs["Sex"]='Female'
EGFRamp_PTEN_CDKN2A_4.obs["Sex"]='Female'
EGFRamp_PTEN_CDKN2A_4b.obs["Sex"]='Female'
EGFRamp_PTEN_CDKN2A_5.obs["Sex"]='Female'
EGFRamp_PTEN_CDKN2A_5b.obs["Sex"]='Female'
EGFRamp_EGFRsnv_PTEN_CDKN2A_1.obs["Sex"]='Male'
EGFRamp_EGFRsnv_PTEN_CDKN2A_1b.obs["Sex"]='Male'
EGFRamp_EGFRsnv_PTEN_CDKN2A_2.obs["Sex"]='Female'
EGFRamp_EGFRsnv_PTEN_CDKN2A_3.obs["Sex"]='Male'
EGFRamp_EGFRsnv_PTEN_CDKN2A_4.obs["Sex"]='Female'
EGFRamp_EGFRsnv_PTEN_CDKN2A_5.obs["Sex"]='Female'
CDKN2A.obs["Sex"]='Male'
PTEN_1.obs["Sex"]='Female'
PTEN_2.obs["Sex"]='Male'
PTEN_3.obs["Sex"]='Female'
PTEN_4.obs["Sex"]='Male'
PTEN_5.obs["Sex"]='Female'
PTEN_6.obs["Sex"]='Male'
PTEN_7.obs["Sex"]='Male'
PTEN_CDKN2A_1.obs["Sex"]='Female'
PTEN_CDKN2A_2.obs["Sex"]='Female'
PTEN_CDKN2A_3.obs["Sex"]='Female'
PTEN_CDKN2A_3b.obs["Sex"]='Female'
PTEN_CDKN2A_4.obs["Sex"]='Male'
PTEN_CDKN2A_5.obs["Sex"]='Female'
PTEN_CDKN2A_6.obs["Sex"]='Female'
PDGFRAamp.obs["Sex"]='Male'

IDH_1.obs["subtype"]='IDH'
IDH_2.obs["subtype"]='IDH'
IDH_EGFRsnv.obs["subtype"]='IDH & EGFRsnv'
IDH_CDKN2A_2.obs["subtype"]='IDH & CDKN2A'
IDH_CDKN2A_3.obs["subtype"]='IDH & CDKN2A'
IDH_CDKN2A_3b.obs["subtype"]='IDH & CDKN2A'
NF1_1.obs["subtype"]='NF1'
NF1_2.obs["subtype"]='NF1'
NF1_3.obs["subtype"]='NF1'
NF1_4.obs["subtype"]='NF1'
NF1_5.obs["subtype"]='NF1'
NF1_EGFRsnv.obs["subtype"]='NF1 & EGFRsnv'
NF1_PTEN_1.obs["subtype"]='NF1 & PTEN'
NF1_PTEN_2.obs["subtype"]='NF1 & PTEN'
NF1_CDKN2A.obs["subtype"]='NF1 & CDKN2A'
NF1_PTEN_CDKN2A_1.obs["subtype"]='NF1, PTEN & CDKN2A'
NF1_PTEN_CDKN2A_1b.obs["subtype"]='NF1, PTEN & CDKN2A'
NF1_PTEN_CDKN2A_2.obs["subtype"]='NF1, PTEN & CDKN2A'
NF1_PTEN_CDKN2A_3.obs["subtype"]='NF1, PTEN & CDKN2A'
NF1_EGFRamp_PTEN_CDKN2A.obs["subtype"]='NF1, EGFRamp, PTEN & CDKN2A'
EGFRsnv_1.obs["subtype"]='EGFRsnv'
EGFRsnv_2.obs["subtype"]='EGFRsnv'
EGFRsnv_3.obs["subtype"]='EGFRsnv'
EGFRsnv_4.obs["subtype"]='EGFRsnv'
EGFRsnv_CDKN2A.obs["subtype"]='EGFRsnv & CDKN2A'
EGFRamp_CDKN2A_1.obs["subtype"]='EGFRamp & CDKN2A'
EGFRamp_CDKN2A_2.obs["subtype"]='EGFRamp & CDKN2A'
EGFRamp_CDKN2A_3.obs["subtype"]='EGFRamp & CDKN2A'
EGFRamp_CDKN2A_4.obs["subtype"]='EGFRamp & CDKN2A'
EGFRamp_CDKN2A_5.obs["subtype"]='EGFRamp & CDKN2A'
EGFRsnv_PTEN.obs["subtype"]='EGFRsnv & PTEN'
EGFRamp_EGFRsnv_CDKN2A_1.obs["subtype"]='EGFRamp, EGFRsnv & CDKN2A'
EGFRamp_EGFRsnv_CDKN2A_2.obs["subtype"]='EGFRamp, EGFRsnv & CDKN2A'
EGFRamp_EGFRsnv_CDKN2A_3.obs["subtype"]='EGFRamp, EGFRsnv & CDKN2A'
EGFRamp_EGFRsnv_CDKN2A_4.obs["subtype"]='EGFRamp, EGFRsnv & CDKN2A'
EGFRamp_EGFRsnv_CDKN2A_5.obs["subtype"]='EGFRamp, EGFRsnv & CDKN2A'
EGFRamp_EGFRsnv_CDKN2A_6.obs["subtype"]='EGFRamp, EGFRsnv & CDKN2A'
EGFRamp_EGFRsnv_CDKN2A_6b.obs["subtype"]='EGFRamp, EGFRsnv & CDKN2A'
EGFRamp_PTEN_CDKN2A_1.obs["subtype"]='EGFRamp, PTEN & CDKN2A'
EGFRamp_PTEN_CDKN2A_2.obs["subtype"]='EGFRamp, PTEN & CDKN2A'
EGFRamp_PTEN_CDKN2A_3.obs["subtype"]='EGFRamp, PTEN & CDKN2A'
EGFRamp_PTEN_CDKN2A_4.obs["subtype"]='EGFRamp, PTEN & CDKN2A'
EGFRamp_PTEN_CDKN2A_4b.obs["subtype"]='EGFRamp, PTEN & CDKN2A'
EGFRamp_PTEN_CDKN2A_5.obs["subtype"]='EGFRamp, PTEN & CDKN2A'
EGFRamp_PTEN_CDKN2A_5b.obs["subtype"]='EGFRamp, PTEN & CDKN2A'
EGFRamp_EGFRsnv_PTEN_CDKN2A_1.obs["subtype"]='EGFRamp, EGFRsnv, PTEN & CDKN2A'
EGFRamp_EGFRsnv_PTEN_CDKN2A_1b.obs["subtype"]='EGFRamp, EGFRsnv, PTEN & CDKN2A'
EGFRamp_EGFRsnv_PTEN_CDKN2A_2.obs["subtype"]='EGFRamp, EGFRsnv, PTEN & CDKN2A'
EGFRamp_EGFRsnv_PTEN_CDKN2A_3.obs["subtype"]='EGFRamp, EGFRsnv, PTEN & CDKN2A'
EGFRamp_EGFRsnv_PTEN_CDKN2A_4.obs["subtype"]='EGFRamp, EGFRsnv, PTEN & CDKN2A'
EGFRamp_EGFRsnv_PTEN_CDKN2A_5.obs["subtype"]='EGFRamp, EGFRsnv, PTEN & CDKN2A'
CDKN2A.obs["subtype"]='CDKN2A'
PTEN_1.obs["subtype"]='PTEN'
PTEN_2.obs["subtype"]='PTEN'
PTEN_3.obs["subtype"]='PTEN'
PTEN_4.obs["subtype"]='PTEN'
PTEN_5.obs["subtype"]='PTEN'
PTEN_6.obs["subtype"]='PTEN'
PTEN_7.obs["subtype"]='PTEN'
PTEN_CDKN2A_1.obs["subtype"]='PTEN & CDKN2A'
PTEN_CDKN2A_2.obs["subtype"]='PTEN & CDKN2A'
PTEN_CDKN2A_3.obs["subtype"]='PTEN & CDKN2A'
PTEN_CDKN2A_3b.obs["subtype"]='PTEN & CDKN2A'
PTEN_CDKN2A_4.obs["subtype"]='PTEN & CDKN2A'
PTEN_CDKN2A_5.obs["subtype"]='PTEN & CDKN2A'
PTEN_CDKN2A_6.obs["subtype"]='PTEN & CDKN2A'
PDGFRAamp.obs["subtype"]='PDGFRAamp'

IDH_1.obs["Location"]="Parietal"
IDH_2.obs["Location"]="Parietal"
IDH_EGFRsnv.obs["Location"]="Parietal"
IDH_CDKN2A_2.obs["Location"]="Frontal"
IDH_CDKN2A_3.obs["Location"]="Frontal"
IDH_CDKN2A_3b.obs["Location"]="Frontal"
NF1_1.obs["Location"]="Frontal"
NF1_2.obs["Location"]="Frontal"
NF1_3.obs["Location"]="Temporal"
NF1_4.obs["Location"]="Frontal"
NF1_5.obs["Location"]="Temporal"
NF1_EGFRsnv.obs["Location"]="Temporal"
NF1_PTEN_1.obs["Location"]="Temporal"
NF1_PTEN_2.obs["Location"]="Other"
NF1_CDKN2A.obs["Location"]="Frontal"
NF1_PTEN_CDKN2A_1.obs["Location"]="Frontal"
NF1_PTEN_CDKN2A_1b.obs["Location"]="Frontal"
NF1_PTEN_CDKN2A_2.obs["Location"]="Frontal"
NF1_PTEN_CDKN2A_3.obs["Location"]="Frontal"
NF1_EGFRamp_PTEN_CDKN2A.obs["Location"]="Frontal"
EGFRsnv_1.obs["Location"]="Temporal"
EGFRsnv_2.obs["Location"]="Frontal"
EGFRsnv_3.obs["Location"]="Frontal"
EGFRsnv_4.obs["Location"]="Frontal"
EGFRsnv_CDKN2A.obs["Location"]="Frontal"
EGFRamp_CDKN2A_1.obs["Location"]="Frontal"
EGFRamp_CDKN2A_2.obs["Location"]="Other"
EGFRamp_CDKN2A_3.obs["Location"]="Temporal"
EGFRamp_CDKN2A_4.obs["Location"]="Frontal"
EGFRamp_CDKN2A_5.obs["Location"]="Other"
EGFRsnv_PTEN.obs["Location"]="Other"
EGFRamp_EGFRsnv_CDKN2A_1.obs["Location"]="Temporal"
EGFRamp_EGFRsnv_CDKN2A_2.obs["Location"]="Frontal"
EGFRamp_EGFRsnv_CDKN2A_3.obs["Location"]="Temporal"
EGFRamp_EGFRsnv_CDKN2A_4.obs["Location"]="Temporal"
EGFRamp_EGFRsnv_CDKN2A_5.obs["Location"]="Frontal"
EGFRamp_EGFRsnv_CDKN2A_6.obs["Location"]="Temporal"
EGFRamp_EGFRsnv_CDKN2A_6b.obs["Location"]="Temporal"
EGFRamp_PTEN_CDKN2A_1.obs["Location"]="Temporal"
EGFRamp_PTEN_CDKN2A_2.obs["Location"]="Occipital"
EGFRamp_PTEN_CDKN2A_3.obs["Location"]="Frontal"
EGFRamp_PTEN_CDKN2A_4.obs["Location"]="Frontal"
EGFRamp_PTEN_CDKN2A_4b.obs["Location"]="Frontal"
EGFRamp_PTEN_CDKN2A_5.obs["Location"]="Frontal"
EGFRamp_PTEN_CDKN2A_5b.obs["Location"]="Frontal"
EGFRamp_EGFRsnv_PTEN_CDKN2A_1.obs["Location"]="Temporal"
EGFRamp_EGFRsnv_PTEN_CDKN2A_1b.obs["Location"]="Temporal"
EGFRamp_EGFRsnv_PTEN_CDKN2A_2.obs["Location"]="Parietal"
EGFRamp_EGFRsnv_PTEN_CDKN2A_3.obs["Location"]="Parietal"
EGFRamp_EGFRsnv_PTEN_CDKN2A_4.obs["Location"]="Temporal"
EGFRamp_EGFRsnv_PTEN_CDKN2A_5.obs["Location"]="Parietal"
CDKN2A.obs["Location"]="Frontal"
PTEN_1.obs["Location"]="Temporal"
PTEN_2.obs["Location"]="Other"
PTEN_3.obs["Location"]="Parietal"
PTEN_4.obs["Location"]="Parietal"
PTEN_5.obs["Location"]="Frontal"
PTEN_6.obs["Location"]="Temporal"
PTEN_7.obs["Location"]="Other"
PTEN_CDKN2A_1.obs["Location"]="Temporal"
PTEN_CDKN2A_2.obs["Location"]="Frontal"
PTEN_CDKN2A_3.obs["Location"]="Parietal"
PTEN_CDKN2A_3b.obs["Location"]="Parietal"
PTEN_CDKN2A_4.obs["Location"]="Temporal"
PTEN_CDKN2A_5.obs["Location"]="Temporal"
PTEN_CDKN2A_6.obs["Location"]="Parietal"
PDGFRAamp.obs["Location"]="Temporal"

IDH_1.obs["sample_label"]="SF9715"
IDH_2.obs["sample_label"]="SF9715_batch2"
IDH_EGFRsnv.obs["sample_label"]="SF5581"
IDH_CDKN2A_2.obs["sample_label"]="SF6621"
IDH_CDKN2A_3.obs["sample_label"]="SF12165"
IDH_CDKN2A_3b.obs["sample_label"]="SF12165_batch2"
NF1_1.obs["sample_label"]="SF2501"
NF1_2.obs["sample_label"]="SF3076"
NF1_3.obs["sample_label"]="SF2979"
NF1_4.obs["sample_label"]="SF6098"
NF1_5.obs["sample_label"]="SF11873"
NF1_EGFRsnv.obs["sample_label"]="SF10514"
NF1_PTEN_1.obs["sample_label"]="SF12115"
NF1_PTEN_2.obs["sample_label"]="SF12751"
NF1_CDKN2A.obs["sample_label"]="SF7062"
NF1_PTEN_CDKN2A_1.obs["sample_label"]="SF4810"
NF1_PTEN_CDKN2A_1b.obs["sample_label"]="SF4810_batch2"
NF1_PTEN_CDKN2A_2.obs["sample_label"]="SF10108"
NF1_PTEN_CDKN2A_3.obs["sample_label"]="SF9871"
NF1_EGFRamp_PTEN_CDKN2A.obs["sample_label"]="SF11344"
EGFRsnv_1.obs["sample_label"]="SF11977"
EGFRsnv_2.obs["sample_label"]="SF2990"
EGFRsnv_3.obs["sample_label"]="SF9494"
EGFRsnv_4.obs["sample_label"]="SF3073"
EGFRsnv_CDKN2A.obs["sample_label"]="SF6809"
EGFRamp_CDKN2A_1.obs["sample_label"]="SF11587"
EGFRamp_CDKN2A_2.obs["sample_label"]="SF12616"
EGFRamp_CDKN2A_3.obs["sample_label"]="SF12407"
EGFRamp_CDKN2A_4.obs["sample_label"]="SF11981"
EGFRamp_CDKN2A_5.obs["sample_label"]="SF12754"
EGFRsnv_PTEN.obs["sample_label"]="SF3243"
EGFRamp_EGFRsnv_CDKN2A_1.obs["sample_label"]="SF11780"
EGFRamp_EGFRsnv_CDKN2A_2.obs["sample_label"]="SF11815"
EGFRamp_EGFRsnv_CDKN2A_3.obs["sample_label"]="SF11916"
EGFRamp_EGFRsnv_CDKN2A_4.obs["sample_label"]="SF12382"
EGFRamp_EGFRsnv_CDKN2A_5.obs["sample_label"]="SF12408"
EGFRamp_EGFRsnv_CDKN2A_6.obs["sample_label"]="SF12243"
EGFRamp_EGFRsnv_CDKN2A_6b.obs["sample_label"]="SF12243_batch2"
EGFRamp_PTEN_CDKN2A_1.obs["sample_label"]="SF10857"
EGFRamp_PTEN_CDKN2A_2.obs["sample_label"]="SF11248"
EGFRamp_PTEN_CDKN2A_3.obs["sample_label"]="SF12460"
EGFRamp_PTEN_CDKN2A_4.obs["sample_label"]="SF10565"
EGFRamp_PTEN_CDKN2A_4b.obs["sample_label"]="SF10565_batch2"
EGFRamp_PTEN_CDKN2A_5.obs["sample_label"]="SF12704"
EGFRamp_PTEN_CDKN2A_5b.obs["sample_label"]="SF12704_batch2"
EGFRamp_EGFRsnv_PTEN_CDKN2A_1.obs["sample_label"]="SF10099"
EGFRamp_EGFRsnv_PTEN_CDKN2A_1b.obs["sample_label"]="SF10099_batch2"
EGFRamp_EGFRsnv_PTEN_CDKN2A_2.obs["sample_label"]="SF10432"
EGFRamp_EGFRsnv_PTEN_CDKN2A_3.obs["sample_label"]="SF12008"
EGFRamp_EGFRsnv_PTEN_CDKN2A_4.obs["sample_label"]="SF12427"
EGFRamp_EGFRsnv_PTEN_CDKN2A_5.obs["sample_label"]="SF12333"
CDKN2A.obs["sample_label"]="SF7025"
PTEN_1.obs["sample_label"]="SF3391"
PTEN_2.obs["sample_label"]="SF9358"
PTEN_3.obs["sample_label"]="SF9372"
PTEN_4.obs["sample_label"]="SF10592"
PTEN_5.obs["sample_label"]="SF3448"
PTEN_6.obs["sample_label"]="SF11857"
PTEN_7.obs["sample_label"]="SF12774"
PTEN_CDKN2A_1.obs["sample_label"]="SF11082"
PTEN_CDKN2A_2.obs["sample_label"]="SF12707"
PTEN_CDKN2A_3.obs["sample_label"]="SF4209"
PTEN_CDKN2A_3b.obs["sample_label"]="SF4209_batch2"
PTEN_CDKN2A_4.obs["sample_label"]="SF9510"
PTEN_CDKN2A_5.obs["sample_label"]="SF11488"
PTEN_CDKN2A_6.obs["sample_label"]="SF4324"
PDGFRAamp.obs["sample_label"]="SF10484"

IDH_1.obs["sample_label2"]="SF9715"
IDH_2.obs["sample_label2"]="SF9715"
IDH_EGFRsnv.obs["sample_label2"]="SF5581"
IDH_CDKN2A_2.obs["sample_label2"]="SF6621"
IDH_CDKN2A_3.obs["sample_label2"]="SF12165"
IDH_CDKN2A_3b.obs["sample_label2"]="SF12165"
NF1_1.obs["sample_label2"]="SF2501"
NF1_2.obs["sample_label2"]="SF3076"
NF1_3.obs["sample_label2"]="SF2979"
NF1_4.obs["sample_label2"]="SF6098"
NF1_5.obs["sample_label2"]="SF11873"
NF1_EGFRsnv.obs["sample_label2"]="SF10514"
NF1_PTEN_1.obs["sample_label2"]="SF12115"
NF1_PTEN_2.obs["sample_label2"]="SF12751"
NF1_CDKN2A.obs["sample_label2"]="SF7062"
NF1_PTEN_CDKN2A_1.obs["sample_label2"]="SF4810"
NF1_PTEN_CDKN2A_1b.obs["sample_label2"]="SF4810"
NF1_PTEN_CDKN2A_2.obs["sample_label2"]="SF10108"
NF1_PTEN_CDKN2A_3.obs["sample_label2"]="SF9871"
NF1_EGFRamp_PTEN_CDKN2A.obs["sample_label2"]="SF11344"
EGFRsnv_1.obs["sample_label2"]="SF11977"
EGFRsnv_2.obs["sample_label2"]="SF2990"
EGFRsnv_3.obs["sample_label2"]="SF9494"
EGFRsnv_4.obs["sample_label2"]="SF3073"
EGFRsnv_CDKN2A.obs["sample_label2"]="SF6809"
EGFRamp_CDKN2A_1.obs["sample_label2"]="SF11587"
EGFRamp_CDKN2A_2.obs["sample_label2"]="SF12616"
EGFRamp_CDKN2A_3.obs["sample_label2"]="SF12407"
EGFRamp_CDKN2A_4.obs["sample_label2"]="SF11981"
EGFRamp_CDKN2A_5.obs["sample_label2"]="SF12754"
EGFRsnv_PTEN.obs["sample_label2"]="SF3243"
EGFRamp_EGFRsnv_CDKN2A_1.obs["sample_label2"]="SF11780"
EGFRamp_EGFRsnv_CDKN2A_2.obs["sample_label2"]="SF11815"
EGFRamp_EGFRsnv_CDKN2A_3.obs["sample_label2"]="SF11916"
EGFRamp_EGFRsnv_CDKN2A_4.obs["sample_label2"]="SF12382"
EGFRamp_EGFRsnv_CDKN2A_5.obs["sample_label2"]="SF12408"
EGFRamp_EGFRsnv_CDKN2A_6.obs["sample_label2"]="SF12243"
EGFRamp_EGFRsnv_CDKN2A_6b.obs["sample_label2"]="SF12243"
EGFRamp_PTEN_CDKN2A_1.obs["sample_label2"]="SF10857"
EGFRamp_PTEN_CDKN2A_2.obs["sample_label2"]="SF11248"
EGFRamp_PTEN_CDKN2A_3.obs["sample_label2"]="SF12460"
EGFRamp_PTEN_CDKN2A_4.obs["sample_label2"]="SF10565"
EGFRamp_PTEN_CDKN2A_4b.obs["sample_label2"]="SF10565"
EGFRamp_PTEN_CDKN2A_5.obs["sample_label2"]="SF12704"
EGFRamp_PTEN_CDKN2A_5b.obs["sample_label2"]="SF12704"
EGFRamp_EGFRsnv_PTEN_CDKN2A_1.obs["sample_label2"]="SF10099"
EGFRamp_EGFRsnv_PTEN_CDKN2A_1b.obs["sample_label2"]="SF10099"
EGFRamp_EGFRsnv_PTEN_CDKN2A_2.obs["sample_label2"]="SF10432"
EGFRamp_EGFRsnv_PTEN_CDKN2A_3.obs["sample_label2"]="SF12008"
EGFRamp_EGFRsnv_PTEN_CDKN2A_4.obs["sample_label2"]="SF12427"
EGFRamp_EGFRsnv_PTEN_CDKN2A_5.obs["sample_label2"]="SF12333"
CDKN2A.obs["sample_label2"]="SF7025"
PTEN_1.obs["sample_label2"]="SF3391"
PTEN_2.obs["sample_label2"]="SF9358"
PTEN_3.obs["sample_label2"]="SF9372"
PTEN_4.obs["sample_label2"]="SF10592"
PTEN_5.obs["sample_label2"]="SF3448"
PTEN_6.obs["sample_label2"]="SF11857"
PTEN_7.obs["sample_label2"]="SF12774"
PTEN_CDKN2A_1.obs["sample_label2"]="SF11082"
PTEN_CDKN2A_2.obs["sample_label2"]="SF12707"
PTEN_CDKN2A_3.obs["sample_label2"]="SF4209"
PTEN_CDKN2A_3b.obs["sample_label2"]="SF4209"
PTEN_CDKN2A_4.obs["sample_label2"]="SF9510"
PTEN_CDKN2A_5.obs["sample_label2"]="SF11488"
PTEN_CDKN2A_6.obs["sample_label2"]="SF4324"
PDGFRAamp.obs["sample_label2"]="SF10484"

IDH_1.obs["Chr7"]="None"
IDH_2.obs["Chr7"]="None"
IDH_EGFRsnv.obs["Chr7"]="Amplification"
IDH_CDKN2A_2.obs["Chr7"]="None"
IDH_CDKN2A_3.obs["Chr7"]="Amplification"
IDH_CDKN2A_3b.obs["Chr7"]="Amplification"
NF1_1.obs["Chr7"]="Amplification"
NF1_2.obs["Chr7"]="Amplification"
NF1_3.obs["Chr7"]="Amplification"
NF1_4.obs["Chr7"]="Amplification"
NF1_5.obs["Chr7"]="None"
NF1_EGFRsnv.obs["Chr7"]="Amplification"
NF1_PTEN_1.obs["Chr7"]="None"
NF1_PTEN_2.obs["Chr7"]="Amplification"
NF1_CDKN2A.obs["Chr7"]="Amplification"
NF1_PTEN_CDKN2A_1.obs["Chr7"]="Amplification"
NF1_PTEN_CDKN2A_1b.obs["Chr7"]="Amplification"
NF1_PTEN_CDKN2A_2.obs["Chr7"]="Amplification"
NF1_PTEN_CDKN2A_3.obs["Chr7"]="None"
NF1_EGFRamp_PTEN_CDKN2A.obs["Chr7"]="Amplification"
EGFRsnv_1.obs["Chr7"]="Amplification"
EGFRsnv_2.obs["Chr7"]="Amplification"
EGFRsnv_3.obs["Chr7"]="None"
EGFRsnv_4.obs["Chr7"]="None"
EGFRsnv_CDKN2A.obs["Chr7"]="Amplification"
EGFRamp_CDKN2A_1.obs["Chr7"]="Amplification"
EGFRamp_CDKN2A_2.obs["Chr7"]="None"
EGFRamp_CDKN2A_3.obs["Chr7"]="None"
EGFRamp_CDKN2A_4.obs["Chr7"]="Amplification"
EGFRamp_CDKN2A_5.obs["Chr7"]="Amplification"
EGFRsnv_PTEN.obs["Chr7"]="Amplification"
EGFRamp_EGFRsnv_CDKN2A_1.obs["Chr7"]="Amplification"
EGFRamp_EGFRsnv_CDKN2A_2.obs["Chr7"]="Amplification"
EGFRamp_EGFRsnv_CDKN2A_3.obs["Chr7"]="Amplification"
EGFRamp_EGFRsnv_CDKN2A_4.obs["Chr7"]="Amplification"
EGFRamp_EGFRsnv_CDKN2A_5.obs["Chr7"]="Amplification"
EGFRamp_EGFRsnv_CDKN2A_6.obs["Chr7"]="Amplification"
EGFRamp_EGFRsnv_CDKN2A_6b.obs["Chr7"]="Amplification"
EGFRamp_PTEN_CDKN2A_1.obs["Chr7"]="Amplification"
EGFRamp_PTEN_CDKN2A_2.obs["Chr7"]="Amplification"
EGFRamp_PTEN_CDKN2A_3.obs["Chr7"]="Amplification"
EGFRamp_PTEN_CDKN2A_4.obs["Chr7"]="Amplification"
EGFRamp_PTEN_CDKN2A_4b.obs["Chr7"]="Amplification"
EGFRamp_PTEN_CDKN2A_5.obs["Chr7"]="Amplification"
EGFRamp_PTEN_CDKN2A_5b.obs["Chr7"]="Amplification"
EGFRamp_EGFRsnv_PTEN_CDKN2A_1.obs["Chr7"]="Amplification"
EGFRamp_EGFRsnv_PTEN_CDKN2A_1b.obs["Chr7"]="Amplification"
EGFRamp_EGFRsnv_PTEN_CDKN2A_2.obs["Chr7"]="Amplification"
EGFRamp_EGFRsnv_PTEN_CDKN2A_3.obs["Chr7"]="Amplification"
EGFRamp_EGFRsnv_PTEN_CDKN2A_4.obs["Chr7"]="None"
EGFRamp_EGFRsnv_PTEN_CDKN2A_5.obs["Chr7"]="Amplification"
CDKN2A.obs["Chr7"]="Amplification"
PTEN_1.obs["Chr7"]="Amplification"
PTEN_2.obs["Chr7"]="Amplification"
PTEN_3.obs["Chr7"]="Amplification"
PTEN_4.obs["Chr7"]="Amplification"
PTEN_5.obs["Chr7"]="Amplification"
PTEN_6.obs["Chr7"]="None"
PTEN_7.obs["Chr7"]="None"
PTEN_CDKN2A_1.obs["Chr7"]="Amplification"
PTEN_CDKN2A_2.obs["Chr7"]="Amplification"
PTEN_CDKN2A_3.obs["Chr7"]="Amplification"
PTEN_CDKN2A_3b.obs["Chr7"]="Amplification"
PTEN_CDKN2A_4.obs["Chr7"]="None"
PTEN_CDKN2A_5.obs["Chr7"]="Amplification"
PTEN_CDKN2A_6.obs["Chr7"]="Amplification"
PDGFRAamp.obs["Chr7"]="Amplification"

IDH_1.obs["Chr10"]="None"
IDH_2.obs["Chr10"]="None"
IDH_EGFRsnv.obs["Chr10"]="Deletion"
IDH_CDKN2A_2.obs["Chr10"]="None"
IDH_CDKN2A_3.obs["Chr10"]="None"
IDH_CDKN2A_3b.obs["Chr10"]="None"
NF1_1.obs["Chr10"]="Deletion"
NF1_2.obs["Chr10"]="Deletion"
NF1_3.obs["Chr10"]="Deletion"
NF1_4.obs["Chr10"]="None"
NF1_5.obs["Chr10"]="None"
NF1_EGFRsnv.obs["Chr10"]="Deletion"
NF1_PTEN_1.obs["Chr10"]="None"
NF1_PTEN_2.obs["Chr10"]="Deletion"
NF1_CDKN2A.obs["Chr10"]="Deletion"
NF1_PTEN_CDKN2A_1.obs["Chr10"]="None"
NF1_PTEN_CDKN2A_1b.obs["Chr10"]="None"
NF1_PTEN_CDKN2A_2.obs["Chr10"]="None"
NF1_PTEN_CDKN2A_3.obs["Chr10"]="None"
NF1_EGFRamp_PTEN_CDKN2A.obs["Chr10"]="Deletion"
EGFRsnv_1.obs["Chr10"]="Deletion"
EGFRsnv_2.obs["Chr10"]="Deletion"
EGFRsnv_3.obs["Chr10"]="None"
EGFRsnv_4.obs["Chr10"]="Deletion"
EGFRsnv_CDKN2A.obs["Chr10"]="Deletion"
EGFRamp_CDKN2A_1.obs["Chr10"]="Deletion"
EGFRamp_CDKN2A_2.obs["Chr10"]="None"
EGFRamp_CDKN2A_3.obs["Chr10"]="None"
EGFRamp_CDKN2A_4.obs["Chr10"]="Deletion"
EGFRamp_CDKN2A_5.obs["Chr10"]="Deletion"
EGFRsnv_PTEN.obs["Chr10"]="Deletion"
EGFRamp_EGFRsnv_CDKN2A_1.obs["Chr10"]="Deletion"
EGFRamp_EGFRsnv_CDKN2A_2.obs["Chr10"]="Deletion"
EGFRamp_EGFRsnv_CDKN2A_3.obs["Chr10"]="Deletion"
EGFRamp_EGFRsnv_CDKN2A_4.obs["Chr10"]="Deletion"
EGFRamp_EGFRsnv_CDKN2A_5.obs["Chr10"]="Deletion"
EGFRamp_EGFRsnv_CDKN2A_6.obs["Chr10"]="Deletion"
EGFRamp_EGFRsnv_CDKN2A_6b.obs["Chr10"]="Deletion"
EGFRamp_PTEN_CDKN2A_1.obs["Chr10"]="Deletion"
EGFRamp_PTEN_CDKN2A_2.obs["Chr10"]="None"
EGFRamp_PTEN_CDKN2A_3.obs["Chr10"]="Deletion"
EGFRamp_PTEN_CDKN2A_4.obs["Chr10"]="Deletion"
EGFRamp_PTEN_CDKN2A_4b.obs["Chr10"]="Deletion"
EGFRamp_PTEN_CDKN2A_5.obs["Chr10"]="None"
EGFRamp_PTEN_CDKN2A_5b.obs["Chr10"]="None"
EGFRamp_EGFRsnv_PTEN_CDKN2A_1.obs["Chr10"]="Deletion"
EGFRamp_EGFRsnv_PTEN_CDKN2A_1b.obs["Chr10"]="Deletion"
EGFRamp_EGFRsnv_PTEN_CDKN2A_2.obs["Chr10"]="Deletion"
EGFRamp_EGFRsnv_PTEN_CDKN2A_3.obs["Chr10"]="Deletion"
EGFRamp_EGFRsnv_PTEN_CDKN2A_4.obs["Chr10"]="None"
EGFRamp_EGFRsnv_PTEN_CDKN2A_5.obs["Chr10"]="Deletion"
CDKN2A.obs["Chr10"]="None"
PTEN_1.obs["Chr10"]="Deletion"
PTEN_2.obs["Chr10"]="Deletion"
PTEN_3.obs["Chr10"]="Deletion"
PTEN_4.obs["Chr10"]="Deletion"
PTEN_5.obs["Chr10"]="Deletion"
PTEN_6.obs["Chr10"]="None"
PTEN_7.obs["Chr10"]="None"
PTEN_CDKN2A_1.obs["Chr10"]="Deletion"
PTEN_CDKN2A_2.obs["Chr10"]="Deletion"
PTEN_CDKN2A_3.obs["Chr10"]="Deletion"
PTEN_CDKN2A_3b.obs["Chr10"]="Deletion"
PTEN_CDKN2A_4.obs["Chr10"]="Deletion"
PTEN_CDKN2A_5.obs["Chr10"]="Deletion"
PTEN_CDKN2A_6.obs["Chr10"]="Deletion"
PDGFRAamp.obs["Chr10"]="Deletion"

IDH_1.var_names_make_unique()
IDH_2.var_names_make_unique()
IDH_EGFRsnv.var_names_make_unique()
IDH_CDKN2A_2.var_names_make_unique()
IDH_CDKN2A_3.var_names_make_unique()
IDH_CDKN2A_3b.var_names_make_unique()
NF1_1.var_names_make_unique()
NF1_2.var_names_make_unique()
NF1_3.var_names_make_unique()
NF1_4.var_names_make_unique()
NF1_5.var_names_make_unique()
NF1_EGFRsnv.var_names_make_unique()
NF1_PTEN_1.var_names_make_unique()
NF1_PTEN_2.var_names_make_unique()
NF1_CDKN2A.var_names_make_unique()
NF1_PTEN_CDKN2A_1.var_names_make_unique()
NF1_PTEN_CDKN2A_1b.var_names_make_unique()
NF1_PTEN_CDKN2A_2.var_names_make_unique()
NF1_PTEN_CDKN2A_3.var_names_make_unique()
NF1_EGFRamp_PTEN_CDKN2A.var_names_make_unique()
EGFRsnv_1.var_names_make_unique()
EGFRsnv_2.var_names_make_unique()
EGFRsnv_3.var_names_make_unique()
EGFRsnv_4.var_names_make_unique()
EGFRsnv_CDKN2A.var_names_make_unique()
EGFRamp_CDKN2A_1.var_names_make_unique()
EGFRamp_CDKN2A_2.var_names_make_unique()
EGFRamp_CDKN2A_3.var_names_make_unique()
EGFRamp_CDKN2A_4.var_names_make_unique()
EGFRamp_CDKN2A_5.var_names_make_unique()
EGFRsnv_PTEN.var_names_make_unique()
EGFRamp_EGFRsnv_CDKN2A_1.var_names_make_unique()
EGFRamp_EGFRsnv_CDKN2A_2.var_names_make_unique()
EGFRamp_EGFRsnv_CDKN2A_3.var_names_make_unique()
EGFRamp_EGFRsnv_CDKN2A_4.var_names_make_unique()
EGFRamp_EGFRsnv_CDKN2A_5.var_names_make_unique()
EGFRamp_EGFRsnv_CDKN2A_6.var_names_make_unique()
EGFRamp_EGFRsnv_CDKN2A_6b.var_names_make_unique()
EGFRamp_PTEN_CDKN2A_1.var_names_make_unique()
EGFRamp_PTEN_CDKN2A_2.var_names_make_unique()
EGFRamp_PTEN_CDKN2A_3.var_names_make_unique()
EGFRamp_PTEN_CDKN2A_4.var_names_make_unique()
EGFRamp_PTEN_CDKN2A_4b.var_names_make_unique()
EGFRamp_PTEN_CDKN2A_5.var_names_make_unique()
EGFRamp_PTEN_CDKN2A_5b.var_names_make_unique()
EGFRamp_EGFRsnv_PTEN_CDKN2A_1.var_names_make_unique()
EGFRamp_EGFRsnv_PTEN_CDKN2A_1b.var_names_make_unique()
EGFRamp_EGFRsnv_PTEN_CDKN2A_2.var_names_make_unique()
EGFRamp_EGFRsnv_PTEN_CDKN2A_3.var_names_make_unique()
EGFRamp_EGFRsnv_PTEN_CDKN2A_4.var_names_make_unique()
EGFRamp_EGFRsnv_PTEN_CDKN2A_5.var_names_make_unique()
CDKN2A.var_names_make_unique()
PTEN_1.var_names_make_unique()
PTEN_2.var_names_make_unique()
PTEN_3.var_names_make_unique()
PTEN_4.var_names_make_unique()
PTEN_5.var_names_make_unique()
PTEN_6.var_names_make_unique()
PTEN_7.var_names_make_unique()
PTEN_CDKN2A_1.var_names_make_unique()
PTEN_CDKN2A_2.var_names_make_unique()
PTEN_CDKN2A_3.var_names_make_unique()
PTEN_CDKN2A_3b.var_names_make_unique()
PTEN_CDKN2A_4.var_names_make_unique()
PTEN_CDKN2A_5.var_names_make_unique()
PTEN_CDKN2A_6.var_names_make_unique()
PDGFRAamp.var_names_make_unique()

sc.external.pp.scrublet(IDH_1)
sc.external.pp.scrublet(IDH_2)
sc.external.pp.scrublet(IDH_EGFRsnv)
sc.external.pp.scrublet(IDH_CDKN2A_2)
sc.external.pp.scrublet(IDH_CDKN2A_3)
sc.external.pp.scrublet(IDH_CDKN2A_3b)
sc.external.pp.scrublet(NF1_1)
sc.external.pp.scrublet(NF1_2)
sc.external.pp.scrublet(NF1_3)
sc.external.pp.scrublet(NF1_4)
sc.external.pp.scrublet(NF1_5)
sc.external.pp.scrublet(NF1_EGFRsnv)
sc.external.pp.scrublet(NF1_PTEN_1)
sc.external.pp.scrublet(NF1_PTEN_2)
sc.external.pp.scrublet(NF1_CDKN2A)
sc.external.pp.scrublet(NF1_PTEN_CDKN2A_1)
sc.external.pp.scrublet(NF1_PTEN_CDKN2A_1b)
sc.external.pp.scrublet(NF1_PTEN_CDKN2A_2)
sc.external.pp.scrublet(NF1_PTEN_CDKN2A_3)
sc.external.pp.scrublet(NF1_EGFRamp_PTEN_CDKN2A)
sc.external.pp.scrublet(EGFRsnv_1)
sc.external.pp.scrublet(EGFRsnv_2)
sc.external.pp.scrublet(EGFRsnv_3)
sc.external.pp.scrublet(EGFRsnv_4)
sc.external.pp.scrublet(EGFRsnv_CDKN2A)
sc.external.pp.scrublet(EGFRamp_CDKN2A_1)
sc.external.pp.scrublet(EGFRamp_CDKN2A_2)
sc.external.pp.scrublet(EGFRamp_CDKN2A_3)
sc.external.pp.scrublet(EGFRamp_CDKN2A_4)
sc.external.pp.scrublet(EGFRamp_CDKN2A_5)
sc.external.pp.scrublet(EGFRsnv_PTEN)
sc.external.pp.scrublet(EGFRamp_EGFRsnv_CDKN2A_1)
sc.external.pp.scrublet(EGFRamp_EGFRsnv_CDKN2A_2)
sc.external.pp.scrublet(EGFRamp_EGFRsnv_CDKN2A_3)
sc.external.pp.scrublet(EGFRamp_EGFRsnv_CDKN2A_4)
sc.external.pp.scrublet(EGFRamp_EGFRsnv_CDKN2A_5)
sc.external.pp.scrublet(EGFRamp_EGFRsnv_CDKN2A_6) 
sc.external.pp.scrublet(EGFRamp_EGFRsnv_CDKN2A_6b)
sc.external.pp.scrublet(EGFRamp_PTEN_CDKN2A_1)
sc.external.pp.scrublet(EGFRamp_PTEN_CDKN2A_2)
sc.external.pp.scrublet(EGFRamp_PTEN_CDKN2A_3)
sc.external.pp.scrublet(EGFRamp_PTEN_CDKN2A_4)
sc.external.pp.scrublet(EGFRamp_PTEN_CDKN2A_4b)
sc.external.pp.scrublet(EGFRamp_PTEN_CDKN2A_5)
sc.external.pp.scrublet(EGFRamp_PTEN_CDKN2A_5b)
sc.external.pp.scrublet(EGFRamp_EGFRsnv_PTEN_CDKN2A_1)
sc.external.pp.scrublet(EGFRamp_EGFRsnv_PTEN_CDKN2A_1b)
sc.external.pp.scrublet(EGFRamp_EGFRsnv_PTEN_CDKN2A_2)
sc.external.pp.scrublet(EGFRamp_EGFRsnv_PTEN_CDKN2A_3)
sc.external.pp.scrublet(EGFRamp_EGFRsnv_PTEN_CDKN2A_4)
sc.external.pp.scrublet(EGFRamp_EGFRsnv_PTEN_CDKN2A_5)
sc.external.pp.scrublet(CDKN2A)
sc.external.pp.scrublet(PTEN_1)
sc.external.pp.scrublet(PTEN_2)
sc.external.pp.scrublet(PTEN_3)
sc.external.pp.scrublet(PTEN_4)
sc.external.pp.scrublet(PTEN_5)
sc.external.pp.scrublet(PTEN_6)
sc.external.pp.scrublet(PTEN_7)
sc.external.pp.scrublet(PTEN_CDKN2A_1)
sc.external.pp.scrublet(PTEN_CDKN2A_2)
sc.external.pp.scrublet(PTEN_CDKN2A_3)
sc.external.pp.scrublet(PTEN_CDKN2A_3b)
sc.external.pp.scrublet(PTEN_CDKN2A_4)
sc.external.pp.scrublet(PTEN_CDKN2A_5)
sc.external.pp.scrublet(PTEN_CDKN2A_6)
sc.external.pp.scrublet(PDGFRAamp)



adata = IDH_1.concatenate(IDH_2,IDH_EGFRsnv,IDH_CDKN2A_2,IDH_CDKN2A_3,IDH_CDKN2A_3b,NF1_1,NF1_2,NF1_3,NF1_4,NF1_5,NF1_EGFRsnv,NF1_PTEN_1,NF1_PTEN_2,NF1_CDKN2A,
                         NF1_PTEN_CDKN2A_1,NF1_PTEN_CDKN2A_1b,NF1_PTEN_CDKN2A_2,NF1_PTEN_CDKN2A_3,NF1_EGFRamp_PTEN_CDKN2A,EGFRsnv_1,EGFRsnv_2,EGFRsnv_3,EGFRsnv_4,
                         EGFRsnv_CDKN2A,EGFRamp_CDKN2A_1,EGFRamp_CDKN2A_2,EGFRamp_CDKN2A_3,EGFRamp_CDKN2A_4,EGFRamp_CDKN2A_5,EGFRsnv_PTEN,EGFRamp_EGFRsnv_CDKN2A_1,
                         EGFRamp_EGFRsnv_CDKN2A_2,EGFRamp_EGFRsnv_CDKN2A_3,EGFRamp_EGFRsnv_CDKN2A_4,EGFRamp_EGFRsnv_CDKN2A_5,EGFRamp_EGFRsnv_CDKN2A_6,
                         EGFRamp_EGFRsnv_CDKN2A_6b,EGFRamp_PTEN_CDKN2A_1,EGFRamp_PTEN_CDKN2A_2,EGFRamp_PTEN_CDKN2A_3,EGFRamp_PTEN_CDKN2A_4,EGFRamp_PTEN_CDKN2A_4b,
                         EGFRamp_PTEN_CDKN2A_5,EGFRamp_PTEN_CDKN2A_5b,EGFRamp_EGFRsnv_PTEN_CDKN2A_1,EGFRamp_EGFRsnv_PTEN_CDKN2A_1b,EGFRamp_EGFRsnv_PTEN_CDKN2A_2,
                         EGFRamp_EGFRsnv_PTEN_CDKN2A_3,EGFRamp_EGFRsnv_PTEN_CDKN2A_4,EGFRamp_EGFRsnv_PTEN_CDKN2A_5,CDKN2A,PTEN_1,PTEN_2,PTEN_3,PTEN_4,PTEN_5,PTEN_6,
                         PTEN_7,PTEN_CDKN2A_1,PTEN_CDKN2A_2,PTEN_CDKN2A_3,PTEN_CDKN2A_3b,PTEN_CDKN2A_4,PTEN_CDKN2A_5,PTEN_CDKN2A_6,PDGFRAamp,join='outer')

adata.var_names_make_unique()
adata.var 
adata

adata.write("Human EGFR datasets complete raw.h5ad")

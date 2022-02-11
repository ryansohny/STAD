# source activate complexheatmap
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import scanpy as sc
from statsmodels.stats.multitest import multipletests
sc.settings.verbosity = 3
plt.rcParams['figure.figsize'] = (5,5)
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#104e8b", "#ffdab9", "#8b0a50"])
%matplotlib
%autoindent


# Step 1. Create Tumor Normal imputed matrix
import pandas as pd
dmr = pd.read_csv("DMR_mat_imputed_corrected.csv", index_col=0)
dmr.iloc[84:, :-1].to_csv("DMR_mat_imputed_corrected_tumor.csv")
dmr.iloc[:84, :-1].to_csv("DMR_mat_imputed_corrected_normal.csv")

# Step 2-1.
dmr_t = sc.AnnData(pd.read_csv("DMR_mat_imputed_corrected_tumor.csv", index_col=0))
dmr_t.raw = dmr_t
dmr_t.layers['Percent_met'] = dmr_t.X
clinic_info = pd.read_csv('/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/2022_WC300_clinical_information_Xadded.csv', index_col=0)
clinic_info[['T.stage', 'N.stage', 'M.stage', 'TNM.stage', 'HER2IHC']] = clinic_info[['T.stage', 'N.stage', 'M.stage', 'TNM.stage', 'HER2IHC']].astype(str)
#clinic_info = clinic_info.iloc[:,np.r_[0:17, -1]]
dmr_t.obs = clinic_info
sc.pp.scale(dmr_t)
sc.tl.pca(dmr_t, n_comps=83, zero_center=True)
sc.pl.pca(dmr_t, color='Lauren')
# sc.pl.pca_variance_ratio(dmr_t, log=True)

pca_variance_t = pd.DataFrame(dmr_t.uns['pca']['variance_ratio'], index=list(map(lambda x: 'PC' + str(x), list(range(1,84)))), columns=['Variance_ratio'])
np.sum(pca_variance_t.values.flatten()[:13])
# 0.7076739
#sc.pp.neighbors(dmr_t, n_neighbors=11, n_pcs=13) ###
#sc.pp.neighbors(dmr_t, n_neighbors=20, n_pcs=13)
sc.pp.neighbors(dmr_t, n_neighbors=11, n_pcs=13)
sc.tl.leiden(dmr_t, resolution=1.0, key_added='leiden_r1')
sc.tl.leiden(dmr_t, resolution=0.5, key_added='leiden_r05')
sc.tl.umap(dmr_t, min_dist=0.5, spread=1.0, n_components=2, alpha=1.0, gamma=1.0, init_pos='spectral', method='umap')
sc.pl.umap(dmr_t, color=['leiden_r1', 'leiden_r05', 'EpiBurden'], add_outline=False, color_map=cmap)

sc.tl.diffmap(dmr_t)
sc.pl.diffmap(dmr_t, color=['leiden_r1'], add_outline=False, color_map=cmap)

start_cell = np.isin(dmr_t.obs['leiden_r1'], '3') # leiden_r1의 3번 cluster
max_start_id = np.argmax(dmr_t.obsm['X_diffmap'][start_cell, 1])
root_id = np.arange(len(start_cell))[start_cell][max_start_id]
dmr_t.uns['iroot'] = root_id

sns.lmplot(data=dmr_t.obs, x='dpt_pseudotime', y='EpiBurden')

deg_tn_protein = pd.read_table("/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/02.RNA-seq/STAD_SNUH_Tumor_leiden_vst_T-N_DEG_protein.txt", index_col="ID")
deg_tn_protein.columns = list(map(lambda x: "X" + x, deg_tn_protein.columns))

rfh = open("DMR_tumor_pseudotime_and_DEG_TN_protein_Correlation.txt", 'w')
rfh.write("ID\tPearsonR\tPvalue\n")

for i in range(len(deg_tn_protein.index)):
    r, p = stats.pearsonr(dmr_t.obs['dpt_pseudotime'], deg_tn_protein.iloc[i])
    r_list.append(r)
    p_list.append(p)

multipletests(pvals=p_list, alpha=0.05, method='fdr_bh')[1]
pt_degtn = pd.DataFrame([r_list, multipletests(pvals=p_list, alpha=0.05, method='fdr_bh')[1]], columns=deg_tn_protein.index, index=['PearsonR', 'Padj']).T

leiden_name_change_dict1 = {'3': 'leiden_A',
                           '1': 'leiden_B',
                           '0': 'leiden_C',
                           '2': 'leiden_D'}
dmr_t.obs['DMR_leiden'] = dmr_t.obs['leiden_r05'].map(lambda x: leiden_name_change_dict1[x]).astype('category')
col_colors1 = list(dict(zip(list(dmr_t.obs['leiden_r05'].value_counts().index), dmr_t.uns['leiden_r05_colors']))[x] for x in dmr_t.obs['leiden_r05'])

dmr_t_met = pd.DataFrame(dmr_t.raw.X.T, index=dmr_t.var.index, columns=dmr_t.obs.index)
g = sns.clustermap(dmr_t_met, method='ward', metric='euclidean', z_score=None, standard_scale=None, cmap=cmap, xticklabels=True, yticklabels=False, col_colors=[col_colors1])

# DMR_leiden clusters and WHO classification
ax = pd.crosstab(dmr_t.obs['WHO'], dmr_t.obs['leiden_r05'], normalize=1).T.plot.bar(stacked=True, color=dict(zip(dmr_t.obs['WHO'].unique(), sns.color_palette("Set2", len(dmr_t.obs['WHO'].unique())))))
ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1.0))
plt.tight_layout()

# DMR_leiden clusters and Lauren's classification
ax = pd.crosstab(df['Lauren'], df['leiden_r05'], normalize=1).T.plot.bar(stacked=True, color=dict(zip(dmr_t.obs['Lauren'].unique(), sns.color_palette("Accent", len(dmr_t.obs['Lauren'].unique())))))
ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1.0))
plt.tight_layout()

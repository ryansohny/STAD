import scanpy as sc
import pandas as pd
dmr = pd.read_csv("DMR_mat_imputed.csv", index_col=0).iloc[:,:-1]
dmr = sc.AnnData(dmr)
dmr.raw = dmr
dmr.layers['Percent_met'] = dmr.X
dmr.obs['TN'] = pd.DataFrame(['Normal']*84 + ['Tumor']*84, index=dmr.obs.index)

sc.pp.scale(dmr)
sc.tl.pca(dmr, n_comps=100, zero_center=True)
sc.pl.pca(dmr, color='TN')
#sc.pl.pca_variance_ratio(dmr, log=True)

### Tests for varying degrees of n_neighbors needed ###
sc.pp.neighbors(dmr, n_neighbors=3, n_pcs=6)
sc.tl.umap(dmr, min_dist=0.5, spread=1.0, n_components=2, alpha=1.0, gamma=1.0, init_pos='spectral', method='umap')

sc.pp.neighbors(dmr, n_neighbors=15, n_pcs=6)
sc.tl.umap(dmr, min_dist=0.5, spread=1.0, n_components=2, alpha=1.0, gamma=1.0, init_pos='spectral', method='umap')


sc.pl.umap(dmr, color='TN')


sc.pp.neighbors(dmr, n_neighbors=10, n_pcs=10) ################################################################ 이게 좋은 것 같음
sc.tl.leiden(dmr, resolution=1.0, key_added='leiden_r1') ## 오우 이거 괜찮

import pandas as pd
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
import seaborn as sns

tumor_mat = pd.read_csv("DMR_mat_tumor_imputed.csv", index_col=0)
normal_mat = pd.read_csv("DMR_mat_normal_imputed.csv", index_col=0)
tn_mat = pd.read_csv("DMR_mat_tumor-normal_imputed.csv", index_col=0)

tumor_hdac5 = pd.read_table("HDAC5_onlyTumor.txt", index_col=0)
normal_hdac5 = pd.read_table("HDAC5_onlyNormal.txt", index_col=0)
diff_hdac5 = pd.read_table("HDAC5_Diff.txt", index_col=0)

rho, pvalue = list(), list()
for col in list(range(len(tumor_mat.columns))):
	res = spearmanr(tumor_hdac5.loc["HDAC5"].values, tumor_mat.iloc[:, col])
	rho.append(res[0])
	pvalue.append(res[1])

tumor_mat.loc["Spearman_R"] = rho
tumor_mat.loc["Pvalue"] = pvalue

rho, pvalue = list(), list()
for col in list(range(len(tn_mat.columns))):
    res = spearmanr(diff_hdac5.loc["HDAC5"].values, tn_mat.iloc[:, col])
    rho.append(res[0])
    pvalue.append(res[1])

tn_mat.loc["Spearman_R"] = rho
tn_mat.loc["Pvalue"] = pvalue


sns.displot(data=tumor_mat.T["Spearman_R"], kind="kde")
plt.savefig("sys.png")

plt.clf()

col_order = list(list(tumor_mat.index[:-2]).index('X' + sample) for sample in tumor_hdac5.loc["HDAC5"].sort_values(ascending=False).index)
col_order2 = list(list(tumor_mat.index[:-2]).index("X" + sample) for sample in diff_hdac5.loc["HDAC5"].sort_values(ascending=False).index)

g = sns.clustermap(tumor_mat.iloc[:-2, :].T.iloc[:, col_order], method="ward", metric="euclidean", z_score=None, standard_scale=None, row_cluster=True, col_cluster=False, xticklabels=True, yticklabels=False)
g_roworder = g.dendrogram_row.reordered_ind
plt.savefig("sys2.pdf")
plt.clf()

sys = tumor_hdac5.loc["HDAC5"].sort_values(ascending=False).plot.bar(rot=0)
sys.axes.xaxis.set_ticklabels([])
plt.savefig("sys3.png")
plt.clf()

sns.heatmap(normal_mat.T.iloc[g_roworder, col_order], xticklabels=False, yticklabels=False)
plt.savefig("sys4.png")
plt.clf()

sys = normal_hdac5.loc["HDAC5"].iloc[col_order].plot.bar(rot=0)
sys.axes.xaxis.set_ticklabels([])
plt.savefig("sys5.png")
plt.clf()

sys = diff_hdac5.loc["HDAC5"].iloc[col_order].plot.bar(rot=0)
sys.axes.xaxis.set_ticklabels([])
plt.tight_layout()
plt.savefig("sys6.png")
plt.clf()


sys = diff_hdac5.loc["HDAC5"].sort_values(ascending=False).plot.bar(rot=0)
sys.axes.xaxis.set_ticklabels([])
plt.tight_layout()
plt.savefig("sys7.png")
plt.clf()

g = sns.clustermap(tumor_mat.iloc[:-2, :].T.iloc[:, col_order2], method="ward", metric="euclidean", z_score=None, standard_scale=None, row_cluster=True, col_cluster=False, xticklabels=True, yticklabels=False)
g_roworder2 = g.dendrogram_row.reordered_ind
plt.tight_layout()
plt.savefig("sys8.png")
plt.clf()

sns.heatmap(normal_mat.T.iloc[g_roworder2, col_order2], xticklabels=False, yticklabels=False)
plt.tight_layout()
plt.savefig("sys9.png")
plt.clf()


g = sns.clustermap(tn_mat.iloc[:-2, :].T.iloc[:, col_order2], method="ward", metric="euclidean", z_score=None, standard_scale=None, row_cluster=True, col_cluster=False, xticklabels=True, yticklabels=False, vmax=30, vmin=-30, cmap="vlag")
plt.tight_layout()
plt.savefig("sys10.png")
plt.clf()



################## ################## 차이 나는 것만 가지고 ################## ################## 

g = sns.clustermap(tumor_mat.iloc[:-2, :].loc[:, abs(tumor_mat.loc['Spearman_R']) > 0.3].T.iloc[:, col_order], method="ward", metric="euclidean", z_score=None, standard_scale=None, row_cluster=True, col_cluster=False, xticklabels=True, yticklabels=False)
g_roworder3 = g.dendrogram_row.reordered_ind
plt.tight_layout()
plt.savefig("sys11.png")
plt.clf()


normal_mat.loc[:, tumor_mat.iloc[:-2, :].loc[:, abs(tumor_mat.loc["Spearman_R"]) > 0.3].columns].T.iloc[:, col_order]

sns.heatmap(normal_mat.loc[:, tumor_mat.iloc[:-2, :].loc[:, abs(tumor_mat.loc["Spearman_R"]) > 0.3].columns].T.iloc[g_roworder3, col_order], xticklabels=False, yticklabels=False)
plt.tight_layout()
plt.savefig("sys12.png")
plt.clf()


g = sns.clustermap(tumor_mat.iloc[:-2, :].loc[:, abs(tumor_mat.loc['Spearman_R']) > 0.3].T.iloc[:, col_order2], method="ward", metric="euclidean", z_score=None, standard_scale=None, row_cluster=True, col_cluster=False, xticklabels=True, yticklabels=False)
g_roworder4 = g.dendrogram_row.reordered_ind
plt.tight_layout()
plt.savefig("sys13.png")
plt.clf()

sns.heatmap(normal_mat.loc[:, tumor_mat.iloc[:-2, :].loc[:, abs(tumor_mat.loc["Spearman_R"]) > 0.3].columns].T.iloc[g_roworder4, col_order2], xticklabels=False, yticklabels=False)
plt.tight_layout()
plt.savefig("sys14.png")
plt.clf()



g = sns.clustermap(tn_mat.iloc[:-2, :].loc[:, abs(tn_mat.loc["Spearman_R"]) > 0.3].T.iloc[:, col_order2], method="ward", metric="euclidean", z_score=None, standard_scale=None, row_cluster=True, col_cluster=False, xticklabels=True, yticklabels=False, vmax=30, vmin=-30, cmap="vlag")
plt.tight_layout()
plt.savefig("sys15.png")
plt.clf()







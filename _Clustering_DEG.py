# source activate complexheatmap
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#104e8b", "#ffdab9", "#8b0a50"])

df = pd.read_csv("/mnt/mone/Project/WC300/02.RNA-seq/02.Kallisto_Analysis/DESeq2_tumor/DEG_results/2022_WC300_clinical_information_Xadded_DMRleidenadded2.csv", index_col=0)

expr = pd.read_table("STAD_SNUH_Tumor_leiden_vst_DEG_protein.txt", index_col="ID")
expr.columns = list(map(lambda x: "X" + x, expr.columns))

color_dict = {
    "leiden_A": "#9467bd",
    "leiden_B": "#ff7f0e",
    "leiden_C": "#1f77b4",
    "leiden_D": "#d62728",
    "leiden_E": '#2ca02c',
}
col_colors=list(color_dict[x] for x in df["DMR_leiden"])
col_colors=df['DMR_leiden'].map(color_dict)
legend = list(Patch(facecolor=color_dict[name]) for name in color_dict)
sns.clustermap(
    expr,
    method="ward",
    metric="euclidean",
    z_score=None,
    standard_scale=0,
    row_cluster=True,
    col_cluster=True,
    xticklabels=True,
    yticklabels=False,
    col_colors=col_colors,
    cmap=cmap,
)
plt.legend(legend, color_dict, title='DMR_leiden', bbox_to_anchor=(1,1), bbox_transform=plt.gcf().transFigure, loc='upper right')
plt.tight_layout()
plt.savefig("HC_STAD_SNUH_Tumor_leiden_vst_DEG_protein.pdf")
plt.clf()









############################################################################################################################################





# source activate complexheatmap
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#104e8b", "#ffdab9", "#8b0a50"])

df = pd.read_csv("/mnt/mone/Project/WC300/02.RNA-seq/02.Kallisto_Analysis/DESeq2_tumor/DEG_results/2022_WC300_clinical_information_Xadded_DMRleidenadded2.csv", index_col=0)

expr = pd.read_table("STAD_SNUH_Tumor_leiden_vst_DEG_protein.txt", index_col="ID")
expr.columns = list(map(lambda x: "X" + x, expr.columns))

leiden_color = {
    "leiden_A": "#9467bd",
    "leiden_B": "#ff7f0e",
    "leiden_C": "#1f77b4",
    "leiden_D": "#d62728",
    "leiden_E": '#2ca02c',
}
#col_colors1 = list(leiden_color[x] for x in df["DMR_leiden"])
col_colors1 = df['DMR_leiden'].map(leiden_color)
col_colors2 = df['LI'].map(dict(zip(df['LI'].unique(), "kw")))
col_colors3 = df['VI'].map(dict(zip(df['VI'].unique(), "kw")))
col_colors4 = df['PI'].map(dict(zip(df['PI'].unique(), "kw")))
col_colors5 = df['WHO'].map(dict(zip(df['WHO'].unique(), sns.color_palette("tab20", len(df['WHO'].unique())))))
col_colors6 = df['WHO_detail'].map(dict(zip(df['WHO_detail'].unique(), sns.color_palette("tab20", len(df['WHO_detail'].unique())))))

col_colors = pd.concat([col_colors1, col_colors2, col_colors3, col_colors4, col_colors5, col_colors6], axis=1)
legend1 = list(Patch(facecolor=leiden_color[name]) for name in leiden_color)
legend5 = list(Patch(facecolor=dict(zip(df['WHO'].unique(), sns.color_palette("tab20", len(df['WHO'].unique()))))[name]) for name in dict(zip(df['WHO'].unique(), sns.color_palette("tab20", len(df['WHO'].unique())))))
legend6 = list(Patch(facecolor=dict(zip(df['WHO_detail'].unique(), sns.color_palette("tab20", len(df['WHO_detail'].unique()))))[name]) for name in dict(zip(df['WHO_detail'].unique(), sns.color_palette("tab20", len(df['WHO_detail'].unique())))))
sns.clustermap(
    expr,
    method="ward",
    metric="euclidean",
    z_score=None,
    standard_scale=0,
    row_cluster=True,
    col_cluster=True,
    xticklabels=True,
    yticklabels=False,
    col_colors=col_colors,
    cmap=cmap,
)
plt.legend(legend1, leiden_color, title='DMR_leiden', bbox_to_anchor=(1,1), bbox_transform=plt.gcf().transFigure, loc='upper right')
plt.legend(legend5, dict(zip(df['WHO'].unique(), sns.color_palette("tab20", len(df['WHO'].unique())))), title='WHO', bbox_to_anchor=(1,1), bbox_transform=plt.gcf().transFigure, loc='upper right')
plt.legend(legend6, dict(zip(df['WHO_detail'].unique(), sns.color_palette("tab20", len(df['WHO_detail'].unique())))), title='WHO_detail', bbox_to_anchor=(1,1), bbox_transform=plt.gcf().transFigure, loc='upper right')
plt.tight_layout()
plt.savefig("HC_STAD_SNUH_Tumor_leiden_vst_DEG_protein.pdf")
plt.clf()

# dfdfdf
import seaborn as sns
sns.boxplot(data=dmr_t.obs[['mean_met_normal', 'mean_met_tumor']], palette={'mean_met_normal': '#00008b','mean_met_tumor': '#8b0000'})
sns.swarmplot(data=dmr_t.obs[['mean_met_normal', 'mean_met_tumor']], color=".2").set(ylabel='Mean Methylation Level')

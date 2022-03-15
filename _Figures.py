####################################################################################
# Boxplot : Mean methylation (Tumor vs Normal)
import pandas as pd
import seaborn as sns

df = dmr_t.obs[['mean_met_normal', 'mean_met_tumor']]
df.rename(columns={'mean_met_normal': 'Normal', 'mean_met_tumor': 'Tumor'}, inplace=True)

sns.boxplot(data=df, palette={'Normal': '#00008b','Tumor': '#8b0000'})
sns.swarmplot(data=df, color="k").set(ylabel='Mean Methylation Level')

new_df = df['Tumor'] - df['Normal']
sns.boxplot(data=new_df, color='lightsalmon')
sns.swarmplot(data=new_df, color='k').set(ylabel='(Tumor - Normal) Mean methylation')

stats.ttest_rel(df['Normal'], df['Tumor'])

####################################################################################
# Violinplot : PMD counts (Tumor vs Normal)
df = pd.read_table("/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/03.WGBS/PMD/PMD_counts.txt", sep='\t')
sns.violinplot(data=df, palette={'Normal':'navy', 'Tumor':'darkred'}, cut=0, scale="count").set_xticklabels(['Normal', 'Tumor'])
sns.stripplot(data=df, color="k").set_title('PMD counts')
